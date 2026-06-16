#### This code has been developed by Jon Sanchez-Valle at BSC ####
# Usage: Rscript 01_build_pool_A_v2.R --disease_a E11 --data_dir ./CaseControlStudy/ --out_dir ./pools/

## Load libraries ##
library(data.table)
library(fst)

## Read arguments ##
args<-commandArgs(trailingOnly = TRUE)
disease_a<-args[which(args == "--disease_a") + 1]
data_dir<-args[which(args == "--data_dir") + 1]
out_dir<-args[which(args == "--out_dir") + 1]

## Stablish thresholds (years, prevalence...)
FOLLOWUP_YEARS<-5
MATCH_WINDOW_DAYS<-61L
MAX_PREVALENCE<-0.20
BATCH_SIZE<-5000L

mem_log<-function(msg) {mem_mb<-sum(gc()[, 2])*8/1024 ; cat(sprintf("[MEM %6.0f MB] %s\n", mem_mb, msg))}
cat(sprintf("Get pool for A = %s", disease_a))
mem_log("inicio")
## Load tables ##
cohort   <-readRDS(file.path(data_dir, "cohort.rds"))
valid_any<-readRDS(file.path(data_dir, "valid_diseases_any.rds"))
valid_A  <-readRDS(file.path(data_dir, "valid_diseases_A.rds"))
mem_log("datos cargados")
if(!(disease_a %in% valid_A)){
  n_casos <-cohort[cod == disease_a, .N]
  n_cohort<-cohort[, uniqueN(idp)]
  cat(sprintf("SKIP: %s prevalence %.1f%% > %.0f%%.\n", disease_a, n_casos / n_cohort * 100, MAX_PREVALENCE * 100))
  quit(status = 0)
}
## Cases ##
cases<-cohort[cod == disease_a & followup_years >= FOLLOWUP_YEARS, .(idp, index_date = dat, followup_end, followup_years, sexe, rangos)]
n_cohort<-cohort[, uniqueN(idp)]
prevalence_A<-nrow(cases)/n_cohort

## MAX_CANDIDATES = 200 for all the diseases with more than 50K cases) ##
## For small and medium keep original values ##
n_cases_A<-nrow(cases)
MAX_CANDIDATES_PER_CASE<-as.integer(
  fifelse(n_cases_A >= 50000, 200L,
  fifelse(prevalence_A < 0.01, 400L,
  fifelse(prevalence_A < 0.02, 600L,
  fifelse(prevalence_A < 0.10, 1000L,
  fifelse(prevalence_A < 0.15, 1200L,
                                1500L))))))
cat(sprintf("Prevalence A: %.2f%% -> MAX_CANDIDATES_PER_CASE = %d\n", prevalence_A * 100, MAX_CANDIDATES_PER_CASE))
mem_log(sprintf("casos: %d", nrow(cases)))
if (nrow(cases) == 0) { cat("SKIP: sin casos.\n"); quit(status = 0) }

## Control candidates ##
patients_with_A<-cohort[cod == disease_a, unique(idp)]
all_first<-cohort[!(idp %in% patients_with_A),
                    {
                      idx<-which.min(dat)
                      .(index_date = dat[idx], index_cod = cod[idx], sexe = sexe[idx], rangos = rangos[idx], followup_end = followup_end[idx], followup_years = followup_years[idx])
                    }, by = idp]
controls_pool<-all_first[followup_years >= FOLLOWUP_YEARS & index_cod %in% valid_any]
rm(all_first, patients_with_A); gc()
mem_log(sprintf("controls_pool: %d candidates", nrow(controls_pool)))
setkey(controls_pool, sexe, rangos, index_date)

## Case-by-case search with batch writing to FST files ##
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_fst<-file.path(out_dir, sprintf("candidates_%s.fst", disease_a))
out_index<-file.path(out_dir, sprintf("candidates_%s_index.rds", disease_a))
n_cases<-nrow(cases)
batch_list<-vector("list", BATCH_SIZE)
batch_i<-0L
n_written<-0L
n_cases_with_candidates<-0L
row_counter<-0L
index_list<-vector("list", n_cases)
index_i<-0L
first_batch<-TRUE
set.seed(1)
cat(sprintf("Looking for candidates for %d casos...\n", n_cases))
progress_every<-10000L
t_start<-proc.time()["elapsed"]

flush_batch<-function(){
  if(batch_i == 0L)return()
  batch_dt<-rbindlist(batch_list[seq_len(batch_i)], use.names = TRUE)
  if(first_batch){
    write_fst(batch_dt, out_fst, compress = 50)
    first_batch<<-FALSE
  }else{
    existing<-read_fst(out_fst, as.data.table = TRUE)
    write_fst(rbind(existing, batch_dt), out_fst, compress = 50)
    rm(existing)
  }
  n_written<<-n_written + nrow(batch_dt)
  batch_list<<-vector("list", BATCH_SIZE)
  batch_i<<-0L
  gc()
}

## NOTE: FST does not support efficient appending for large files. ##
## To preserve compatibility and performance, we accumulate all results in memory and write them once at the end. For XL diseases with MAX = 200,
## this is feasible: 547K cases × 200 candidates × ~100 bytes ≈ 11 GB of RAM. ##
## For diseases with a prevalence below 5% (MAX = 1000), we continue using CSV. ##

if(MAX_CANDIDATES_PER_CASE <= 200){
  cat("XL mode: accumulating data in memory for efficient FST writing...\n")
  all_batches<-vector("list", n_cases)
  for(i in seq_len(n_cases)){
    if(i %% progress_every == 0){
      t_elapsed  <-proc.time()["elapsed"] - t_start
      t_remaining<-t_elapsed/i*(n_cases - i)
      mem_log(sprintf("caso %d / %d - ETA: %.0f min", i, n_cases, t_remaining / 60))
    }
    caso<-cases[i]
    date_lo<-caso$index_date - MATCH_WINDOW_DAYS
    date_hi<-caso$index_date + MATCH_WINDOW_DAYS
    cands<-controls_pool[sexe == caso$sexe & rangos == caso$rangos & index_date >= date_lo & index_date <= date_hi]
    if(nrow(cands) == 0){next}
    if(nrow(cands) > MAX_CANDIDATES_PER_CASE){cands<-cands[sample(.N, MAX_CANDIDATES_PER_CASE)]}
    n_cases_with_candidates<-n_cases_with_candidates + 1L
    row_start<-n_written + 1L
    n_rows<-nrow(cands)
    all_batches[[i]]<-data.table(case_idp = caso$idp, case_index_date = caso$index_date, case_followup_end = caso$followup_end,
      ctrl_idp = cands$idp, ctrl_index_date = cands$index_date, ctrl_index_cod = cands$index_cod, ctrl_followup_end = cands$followup_end,
      sexe = caso$sexe, rangos = caso$rangos)
    index_i<-index_i + 1L
    index_list[[index_i]]<-data.table(case_idp = caso$idp, row_start = row_start, row_end = row_start + n_rows - 1L, n_rows = n_rows)
    n_written<-n_written + n_rows
  }
  mem_log("Search completed. Writing FST file...")
  all_dt<-rbindlist(Filter(Negate(is.null), all_batches))
  rm(all_batches); gc()
  write_fst(all_dt, out_fst, compress = 50)
  rm(all_dt); gc()
}else{
  ## Normal mode (prevalence < 5%): write CSV as in v1 and create the index by reading the CSV at the end ##
  out_csv<-file.path(out_dir, sprintf("candidates_%s.csv", disease_a))
  write_header<-TRUE
  for(i in seq_len(n_cases)){
    if(i %% progress_every == 0){
      t_elapsed  <-proc.time()["elapsed"]-t_start
      t_remaining<-t_elapsed/i*(n_cases - i)
      mem_log(sprintf("case %d / %d - ETA: %.0f min", i, n_cases, t_remaining / 60))
    }
    caso <-cases[i]
    date_lo<-caso$index_date-MATCH_WINDOW_DAYS
    date_hi<-caso$index_date+MATCH_WINDOW_DAYS
    cands<-controls_pool[sexe == caso$sexe & rangos == caso$rangos & index_date >= date_lo & index_date <= date_hi]
    if(nrow(cands) == 0) next
    if(nrow(cands) > MAX_CANDIDATES_PER_CASE){cands<-cands[sample(.N, MAX_CANDIDATES_PER_CASE)]}
    n_cases_with_candidates<-n_cases_with_candidates + 1L
    batch_i<-batch_i + 1L
    batch_list[[batch_i]]<-data.table(case_idp = caso$idp, case_index_date = caso$index_date, case_followup_end = caso$followup_end,
      ctrl_idp = cands$idp, ctrl_index_date = cands$index_date, ctrl_index_cod = cands$index_cod, ctrl_followup_end = cands$followup_end, sexe = caso$sexe, rangos = caso$rangos)
    if (batch_i == BATCH_SIZE){
      batch_dt<-rbindlist(batch_list[seq_len(batch_i)], use.names = TRUE)
      fwrite(batch_dt, out_csv, append = !write_header, col.names = write_header)
      write_header<<-FALSE
      n_written<-n_written + nrow(batch_dt)
      batch_list<-vector("list", BATCH_SIZE)
      batch_i<-0L
      gc()
    }
  }
  ## Final flush ##
  if(batch_i > 0){
    batch_dt<-rbindlist(batch_list[seq_len(batch_i)], use.names = TRUE)
    fwrite(batch_dt, out_csv, append = TRUE, col.names = FALSE)
    n_written<-n_written + nrow(batch_dt)
    gc()
  }
  out_fst<-out_csv
}
mem_log(sprintf("Search completed - %d rows written", n_written))
if(n_written == 0){cat("SKIP: ningun caso tiene controles candidatos.\n") ; quit(status = 0)}

## Create an index (only for XL mode with FST) ##
if(MAX_CANDIDATES_PER_CASE <= 200 && index_i > 0){
  cat("Creating index by case_idp...")
  index_dt<-rbindlist(index_list[seq_len(index_i)])
  saveRDS(index_dt, out_index)
  cat(sprintf("Index saved: %d unique cases\n", nrow(index_dt)))
} else if(MAX_CANDIDATES_PER_CASE > 200){
  cat("Creating index from CSV...")
  out_csv<-file.path(out_dir, sprintf("candidates_%s.csv", disease_a))
  idx<-fread(out_csv, select = "case_idp")
  idx[, row_num := .I]
  index_dt <-idx[, .(row_start = first(row_num), row_end = last(row_num), n_rows = .N), by = case_idp]
  rm(idx); gc()
  saveRDS(index_dt, out_index)
  cat(sprintf("Index saved: %d unique cases\n", nrow(index_dt)))
}

## Candidate diagnosis ##
n_cands_per_case<-index_dt[, .(n_cands = n_rows)]
cat(sprintf("\nCandidates by case:\n"))
cat(sprintf("Median : %d\n", as.integer(median(n_cands_per_case$n_cands))))
cat(sprintf("Min  : %d\n", min(n_cands_per_case$n_cands)))
cat(sprintf("Max  : %d\n", max(n_cands_per_case$n_cands)))
cat(sprintf("Saturated (= MAX): %d (%.1f%%)\n", sum(n_cands_per_case$n_cands == MAX_CANDIDATES_PER_CASE), mean(n_cands_per_case$n_cands == MAX_CANDIDATES_PER_CASE) * 100))

## Save metadata and complete pool ##
out_meta<-file.path(out_dir, sprintf("pool_meta_%s.rds", disease_a))
saveRDS(
  list(disease_a = disease_a, n_cases = nrow(cases), n_candidates = n_written, max_candidates_per_case = MAX_CANDIDATES_PER_CASE, prevalence_A = prevalence_A, format = ifelse(MAX_CANDIDATES_PER_CASE <= 200, "fst", "csv")),
  out_meta
)
out_full<-file.path(out_dir, sprintf("pool_full_%s.rds", disease_a))
saveRDS(list(cases = cases, controls_pool = controls_pool, disease_a = disease_a), out_full)

mem_log("Saved")

gb_size<-file.size(out_fst) / 1e9
cat(sprintf("\nPool summary %s:\n", disease_a))
cat(sprintf("Total cases: %d\n", nrow(cases)))
cat(sprintf("Cases with >= 1 candidate: %d (%.1f%%)\n", n_cases_with_candidates, n_cases_with_candidates / nrow(cases) * 100))
cat(sprintf("Row candidates: %d\n", n_written))
cat(sprintf("MAX_CANDIDATES_PER_CASE: %d\n", MAX_CANDIDATES_PER_CASE))
cat(sprintf("Format: %s\n",
            ifelse(MAX_CANDIDATES_PER_CASE <= 200, "fst", "csv")))
cat(sprintf("File size: %.2f GB\n", gb_size))
cat(sprintf("pool_meta: %s\n", out_meta))
cat(sprintf("pool_full: %s\n", out_full))
cat(sprintf("index: %s\n", out_index))
