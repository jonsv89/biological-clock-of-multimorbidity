#### This code has been developed by Jon Sanchez-Valle at BSC ####
# Usage: Rscript 02_analyze_AB.R --disease_a E11 --disease_b I10 --data_dir ./CaseControlStudy/ --pool_dir ./pools/ --out_dir ./matched/

## Load libraries ##
library(data.table)
library(fst)

## Read arguments ##
args<-commandArgs(trailingOnly = TRUE)
disease_a<-args[which(args == "--disease_a") + 1]
disease_b<-args[which(args == "--disease_b") + 1]
data_dir<-args[which(args == "--data_dir") + 1]
pool_dir<-args[which(args == "--pool_dir") + 1]
out_dir<-args[which(args == "--out_dir") + 1]

## Stablish thresholds (number of controls...) ##
N_CONTROLS <-5L
MAX_CTRL_COD_FRAC<-0.30
CHUNK_SIZE <-500000L
mem_log<-function(msg){mem_mb<-sum(gc()[, 2])*8/1024 ; cat(sprintf("[MEM %6.0f MB] %s\n", mem_mb, msg))}
mem_log("inicio")

## Load the complete A pool ##
full_file<-file.path(pool_dir, sprintf("pool_full_%s.rds", disease_a))
if(!file.exists(full_file)){cat(sprintf("SKIP: pool_full does not exist for %s.\n", disease_a)); quit(status = 0)}
pool<-readRDS(full_file)
cases<-pool$cases
controls_pool<-pool$controls_pool
rm(pool); gc()
mem_log("pool loaded")

## Detect format: fst (XL) or csv (normal) ##
fst_file<-file.path(pool_dir, sprintf("candidates_%s.fst", disease_a))
csv_file<-file.path(pool_dir, sprintf("candidates_%s.csv", disease_a))
index_file<-file.path(pool_dir, sprintf("candidates_%s_index.rds", disease_a))
use_fst<-file.exists(fst_file) && file.exists(index_file)
use_csv<-file.exists(csv_file)

if(!use_fst && !use_csv){cat(sprintf("SKIP: no existe fichero de candidatos para %s.", disease_a)); quit(status = 0)}
cat(sprintf("Candidates' format: %s\n", ifelse(use_fst, "fst (XL)", "csv")))

## Load B diagnoses ##
cohort<-readRDS(file.path(data_dir, "cohort.rds"))
diag_B<-cohort[cod == disease_b, .(idp, dat_B = dat)]
rm(cohort); gc()
mem_log(sprintf("B diagnoses: %d patients", nrow(diag_B)))

## Filter cases ##
cases_B<-diag_B[idp %in% cases$idp]
cases2<-merge(cases, cases_B, by = "idp", all.x = TRUE)
cases2<-cases2[is.na(dat_B) | dat_B > index_date]
cases2[, dat_B := NULL]
cat(sprintf("Elegible cases for B: %d / %d\n", nrow(cases2), nrow(cases)))

if(nrow(cases2) == 0){cat("SKIP: no elegible cases after filtering previous B"); quit(status = 0)}

## Identify non-elegible controls for B ##
bad_index_B<-controls_pool[index_cod == disease_b, idp]
ctrl_B<-merge(controls_pool[, .(idp, index_date)], diag_B[idp %in% controls_pool$idp], by = "idp", all.x = TRUE)
bad_prior_B<-ctrl_B[!is.na(dat_B) & dat_B <= index_date, idp]
ineligible_ctrls<-union(bad_index_B, bad_prior_B)
cat(sprintf("Ineligible controls for B: %d (index=B: %d, previous B: %d)\n", length(ineligible_ctrls), length(bad_index_B), length(bad_prior_B)))
eligible_case_idps<-cases2$idp

## Read candidates ##
t_read<-proc.time()["elapsed"]
if(use_fst){
  cat("Reading candidates through FST index...")
  index_dt<-readRDS(index_file)
  eligible_idx<-index_dt[case_idp %in% eligible_case_idps]
  cat(sprintf("Cases in index: %d / %d elegible\n", nrow(eligible_idx), length(eligible_case_idps)))
  if(nrow(eligible_idx) == 0){cat("SKIP: no elegible case has candidates"); quit(status = 0)}
  setorder(eligible_idx, row_start)
  ## Merge contiguous ranges to reduce the number of reads ##
  eligible_idx[, gap := row_start - shift(row_end, fill = 0) - 1]
  eligible_idx[, group := cumsum(gap > 100)]
  ranges<-eligible_idx[, .(row_from = first(row_start), row_to = last(row_end), cases = list(case_idp)), by = group]
  cat(sprintf("Reading %d row ranges (from %d cases)...\n", nrow(ranges), nrow(eligible_idx)))
  chunks<-vector("list", nrow(ranges))
  for(i in seq_len(nrow(ranges))){
    chunk<-as.data.table(read_fst(fst_file, from = ranges$row_from[i], to = ranges$row_to[i]))
    ## Filter only elegible case_idp ##
    chunks[[i]]<-chunk[case_idp %in% eligible_case_idps]
  }
  candidates<-rbindlist(chunks, use.names = TRUE)
  rm(chunks); gc()
}else{
  cat("Reading candidates through CSV by chunks...")
  col_names<-names(fread(csv_file, nrows = 0))
  n_total<-nrow(fread(csv_file, select = 1L))
  n_chunks <-ceiling(n_total / CHUNK_SIZE)
  cat(sprintf("Total rows: %d | Chunks: %d\n", n_total, n_chunks))
  filtered_list<-vector("list", n_chunks)
  for(ch in seq_len(n_chunks)){
    chunk<-fread(csv_file, skip = (ch - 1) * CHUNK_SIZE + 1, nrows = CHUNK_SIZE, header = FALSE, col.names = col_names)
    filtered_list[[ch]]<-chunk[case_idp %in% eligible_case_idps]
    if(ch %% 10 == 0){cat(sprintf("Chunk %d / %d processed\n", ch, n_chunks))}
  }
  candidates<-rbindlist(filtered_list, use.names = TRUE)
  rm(filtered_list); gc()
}
cat(sprintf("Read in %.1f min\n", (proc.time()["elapsed"] - t_read) / 60))

## Filter ineligible controls ##
candidates<-candidates[!(ctrl_idp %in% ineligible_ctrls)]
mem_log(sprintf("Candidates after filtering: %d rows", nrow(candidates)))

# Ensure Date types
date_cols<-c("case_index_date", "case_followup_end", "ctrl_index_date", "ctrl_followup_end")
candidates[, (date_cols) := lapply(.SD, as.Date), .SDcols = date_cols]
if(nrow(candidates) == 0){cat("SKIP: without candidates after B filters"); quit(status = 0)}

## Control assignment without reuse ##
n_cands_per_case<-candidates[, .(n_avail = .N), by = case_idp]
case_order<-n_cands_per_case[order(n_avail), case_idp]
case_order<-c(case_order, setdiff(cases2$idp, n_cands_per_case$case_idp))
used_set<-new.env(hash = TRUE, size = nrow(candidates), parent = emptyenv())
setkey(candidates, case_idp)
set.seed(1)
matched_list<-vector("list", length(case_order))
cat(sprintf("Assigning controls for %d cases...\n", length(case_order)))
progress_every<-5000L
t_start<-proc.time()["elapsed"]
for(i in seq_along(case_order)){
  if(i %% progress_every == 0){
    t_elapsed<-proc.time()["elapsed"] - t_start
    t_remaining<-t_elapsed/i*(length(case_order) - i)
    cat(sprintf("[%d / %d] %.1f%% - %.0f seg - ETA: %.0f min\n", i, length(case_order), i/length(case_order)*100, t_elapsed, t_remaining/60))
  }
  cid<-case_order[i]
  available<-candidates[.(cid)]
  if(nrow(available) == 0){next}
  is_used<-vapply(available$ctrl_idp, exists, logical(1), envir = used_set, inherits = FALSE)
  available<-available[!is_used]
  if(nrow(available) == 0){next}
  available<-available[sample(.N)]
  selected<-data.table()
  for(j in seq_len(nrow(available))){
    if(nrow(selected) >= N_CONTROLS){break}
    cand<-available[j]
    n_so_far<-nrow(selected)
    n_this_cod<-if(n_so_far > 0) sum(selected$ctrl_index_cod == cand$ctrl_index_cod) else 0L
    max_allowed<-max(1L, floor((n_so_far + 1L) * MAX_CTRL_COD_FRAC))
    if(n_this_cod < max_allowed){selected<-rbind(selected, cand)}
  }
  if(nrow(selected) == 0){next}
  for(uid in selected$ctrl_idp){assign(uid, TRUE, envir = used_set)}
  matched_list[[i]]<-selected
}

matched<-rbindlist(matched_list, use.names = TRUE)
rm(matched_list, candidates, used_set); gc()
mem_log(sprintf("matching completed: %d pares", nrow(matched)))
if (nrow(matched) == 0) {cat("SKIP: without pairs after matching"); quit(status = 0)}

## Add B diagnosis dates ##
case_B_dates<-diag_B[idp %in% matched$case_idp]
setnames(case_B_dates, c("idp", "dat_B"), c("case_idp", "case_dat_B"))
matched<-merge(matched, case_B_dates, by = "case_idp", all.x = TRUE)

ctrl_B_dates<-diag_B[idp %in% matched$ctrl_idp]
setnames(ctrl_B_dates, c("idp", "dat_B"), c("ctrl_idp", "ctrl_dat_B"))
matched<-merge(matched, ctrl_B_dates, by = "ctrl_idp", all.x = TRUE)
mem_log("B diagnosis dates added")

## Separate complete and incomplete matches ##
n_per_case<-matched[, .(n_ctrl = .N), by = case_idp]
cases_complete<-n_per_case[n_ctrl == N_CONTROLS, case_idp]
cases_incomplete<-n_per_case[n_ctrl <  N_CONTROLS, case_idp]
matched_complete<-matched[case_idp %in% cases_complete]
matched_incomplete<-matched[case_idp %in% cases_incomplete]

## Coverage diagnostics file ##
dist_table<-data.table(n_ctrl = 0:N_CONTROLS)
n_zero<-nrow(cases2) - matched[, uniqueN(case_idp)]
dist_zero<-data.table(n_ctrl = 0L, n_cases = n_zero)
dist_matched<-n_per_case[, .N, by = n_ctrl]
setnames(dist_matched, "N", "n_cases")
dist_full<-merge(dist_table, rbind(dist_zero, dist_matched), by = "n_ctrl", all.x = TRUE)
dist_full[is.na(n_cases), n_cases := 0L]
setorder(dist_full, n_ctrl)

dist_full[, pct := round(n_cases / nrow(cases2) * 100, 2)]
dist_full[, pct_cum := round(cumsum(n_cases) / nrow(cases2) * 100, 2)]
dist_full[, disease_a := disease_a]
dist_full[, disease_b := disease_b]
dist_full[, n_cases_eligible := nrow(cases2)]
dist_full[, n_cases_complete := length(cases_complete)]
dist_full[, n_cases_incomplete := length(cases_incomplete)]
dist_full[, n_cases_zero := n_zero]
dist_full[, pct_complete := round(length(cases_complete) / nrow(cases2) * 100, 2)]

setcolorder(dist_full, c("disease_a", "disease_b", "n_cases_eligible", "n_cases_complete", "pct_complete", "n_cases_incomplete", "n_cases_zero", "n_ctrl", "n_cases", "pct", "pct_cum"))

## Save ##
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_complete<-file.path(out_dir, sprintf("matched_complete_%s_%s.rds", disease_a, disease_b))
saveRDS(list(matched = matched_complete, disease_a = disease_a, disease_b = disease_b), out_complete)

out_all<-file.path(out_dir, sprintf("matched_all_%s_%s.rds", disease_a, disease_b))
saveRDS(list(matched = rbind(matched_complete, matched_incomplete), disease_a = disease_a, disease_b = disease_b), out_all)

out_diag<-file.path(out_dir, sprintf("coverage_%s_%s.csv", disease_a, disease_b))
fwrite(dist_full, out_diag)

mem_log("Files saved")

cat(sprintf("\nSummary %s -> %s:\n", disease_a, disease_b))
cat(sprintf("Elegible cases: %d\n", nrow(cases2)))
cat(sprintf("Cases with 0 controls: %d (%.1f%%)\n", n_zero, n_zero / nrow(cases2) * 100))
cat(sprintf("Cases with 1-4 controls: %d (%.1f%%)\n", length(cases_incomplete), length(cases_incomplete) / nrow(cases2) * 100))
cat(sprintf("Cases with 5 controls: %d (%.1f%%)\n", length(cases_complete), length(cases_complete) / nrow(cases2) * 100))
cat(sprintf("\n  Distribucion detallada:\n"))
print(dist_full[, .(n_ctrl, n_cases, pct, pct_cum)])
cat(sprintf("Saved:"))
cat(sprintf("Complete (5 ctrl): %s\n", out_complete))
cat(sprintf("All (>=1 ctrl): %s\n", out_all))
cat(sprintf("Coverage CSV: %s\n", out_diag))
