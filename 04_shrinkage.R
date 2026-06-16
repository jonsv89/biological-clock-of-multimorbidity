#### This code has been developed by Jon Sanchez-Valle at BSC ####
## Usage 1: Rscript 04_shrinkage.R --results_dir ./Results --method rr_contingency --prescreening events --out_dir ./Results/shrinkage_events
## Usage 2: Rscript 04_shrinkage.R --results_dir ./Results --method rr_contingency --prescreening westergaard --out_dir ./Results/shrinkage_westergaard

## Load libraries ##
library(data.table)
library(ashr)

## Read arguments ##
args<-commandArgs(trailingOnly = TRUE)
results_dir<-args[which(args == "--results_dir")+1]
method<-args[which(args == "--method")+1]
out_dir<-args[which(args == "--out_dir")+1]
prescreening<-if("--prescreening" %in% args) args[which(args == "--prescreening") + 1] else "events"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## Validate arguments
if(!method %in% c("rr_contingency", "clogit")){stop("--method debe ser 'rr_contingency' o 'clogit'")}
if(!prescreening %in% c("westergaard", "events")){stop("--prescreening debe ser 'westergaard' o 'events'")}

## Configuration depending on the method ##
if (method == "rr_contingency") {
  in_prefix<-"RR_contingency"
  out_prefix<-"RR_contingency"
  log_col<-"log_RR"
  se_col<-"SE_RR"
  shrunk_col<-"RR_shrunk"
  ci_low_shrunk<-"CI_low_RR_shrunk"
  ci_high_shrunk<-"CI_high_RR_shrunk"
  log_shrunk_col<-"log_RR_shrunk"
  se_shrunk_col<-"SE_RR_shrunk"
} else {
  in_prefix<-"OR_clogit"
  out_prefix<-"OR_clogit"
  log_col<-"log_OR"
  se_col<-"SE"
  shrunk_col<-"OR_shrunk"
  ci_low_shrunk <-"CI_low_OR_shrunk"
  ci_high_shrunk<-"CI_high_OR_shrunk"
  log_shrunk_col<-"log_OR_shrunk"
  se_shrunk_col <-"SE_OR_shrunk"
}

windows<-c("0_1", "0_2", "0_3", "0_4", "0_5", "1_2", "2_3", "3_4", "4_5")
pops<-c("both", "women", "men")

cat(sprintf("Method: %s\n", method))
cat(sprintf("Prescreening: %s\n", prescreening))
cat(sprintf("Input: %s\n", results_dir))
cat(sprintf("Output: %s\n", out_dir))
cat(sprintf("Groups: %d windows x %d populations = %d\n", length(windows), length(pops), length(windows) * length(pops)))

## Process each window × population combination ##
for(w in windows){
  for(pop in pops){
    in_file<-file.path(results_dir, sprintf("%s_%s_%s.txt", in_prefix, pop, w))
    if(!file.exists(in_file)){cat(sprintf("SKIP: doest not exist %s\n", basename(in_file))); next}
    out_file<-file.path(out_dir, sprintf("%s_%s_%s_shrunk.txt", out_prefix, pop, w))
    if(file.exists(out_file)){cat(sprintf("SKIP: already exists %s\n", basename(out_file))); next}
    cat(sprintf("\nProcessing %s / %s / %s...\n", method, pop, w))
    t_start<-proc.time()["elapsed"]
    dt<-fread(in_file)
    ## Force numeric columns (they may be read as character due to Inf or 0 values) ##
    num_cols<-c("OR", "CI_low", "CI_high", "log_OR", "SE", "p", "RR", "CI_low_RR", "CI_high_RR", "log_RR", "SE_RR", "p_RR", "CI_low_OR", "CI_high_OR", "SE_OR", "p_OR")
    for(col in intersect(num_cols, names(dt))){if(!is.numeric(dt[[col]])){dt[, (col) := suppressWarnings(as.numeric(get(col)))]}}
    cat(sprintf("Total rows: %d\n", nrow(dt)))
    ## Filter based on prescreening ##
    if(method == "rr_contingency"){
      if(prescreening == "westergaard"){
        ## Equivalent to the prescreening in Westergaard et al. 2019 ##
        dt[, valid := is.finite(get(log_col)) & is.finite(get(se_col)) & get(se_col) > 0 & !is.na(CI_low_RR) & CI_low_RR >= 1.01 & cases_event > 0 & ctrl_events > 0]
      }else{
        ## At least one event in cases or controls ##
        dt[, valid := is.finite(get(log_col)) & is.finite(get(se_col)) & get(se_col) > 0 & (cases_event > 0 | ctrl_events > 0)]
      }
    }else{
      ## Clogit: filter only finite values ##
      dt[, valid := is.finite(get(log_col)) & is.finite(get(se_col)) & get(se_col) > 0]
    }
    n_valid <-dt[valid == TRUE, .N]
    n_invalid<-dt[valid == FALSE, .N]
    cat(sprintf("Valid for shrinkage: %d\n", n_valid))
    cat(sprintf("Excluded: %d\n", n_invalid))

    ## Initialize shrinkage columns as NA ##
    dt[, (log_shrunk_col) := NA_real_]
    dt[, (se_shrunk_col) := NA_real_]
    dt[, (shrunk_col) := NA_real_]
    dt[, (ci_low_shrunk) := NA_real_]
    dt[, (ci_high_shrunk) := NA_real_]
    dt[, lfsr := NA_real_]
    if(n_valid < 10){
      cat(sprintf("SKIP: very few valids (%d)\n", n_valid))
      dt[, valid := NULL]
      fwrite(dt, out_file, sep = "\t", quote = FALSE)
      next
    }

    ## Apply ashr ##
    sub<-dt[valid == TRUE]
    tryCatch({
      ash_fit<-ash(sub[[log_col]], sub[[se_col]], mixcompdist = "normal")
      dt[valid == TRUE, (log_shrunk_col) := get_pm(ash_fit)]
      dt[valid == TRUE, (se_shrunk_col)  := get_psd(ash_fit)]
      dt[valid == TRUE, lfsr := get_lfsr(ash_fit)]
      dt[valid == TRUE, (shrunk_col) := exp(get(log_shrunk_col))]
      dt[valid == TRUE, (ci_low_shrunk) := exp(get(log_shrunk_col) - 1.96*get(se_shrunk_col))]
      dt[valid == TRUE, (ci_high_shrunk) := exp(get(log_shrunk_col) + 1.96*get(se_shrunk_col))]
      n_sig<-dt[lfsr < 0.05 & !is.na(lfsr), .N]
      cat(sprintf("Significant (lfsr < 0.05): %d\n", n_sig))
    }, error = function(e) {cat(sprintf("ERROR in ashr: %s\n", e$message))})
    dt[, valid := NULL]
    fwrite(dt, out_file, sep = "\t", quote = FALSE)
    t_elapsed<-(proc.time()["elapsed"] - t_start) / 60
    cat(sprintf("Saved in %.1f min: %s\n", t_elapsed, basename(out_file)))
  }
}

## Final summary ##
cat(sprintf("\n%s\n", strrep("=", 60)))
cat(sprintf("SUMMARY SHRINKAGE - method: %s, prescreening: %s\n", method, prescreening))
cat(sprintf("%s\n", strrep("=", 60)))
shrunk_files<-list.files(out_dir, pattern = sprintf("^%s_.*_shrunk\\.txt$", out_prefix))
cat(sprintf("Files generated: %d\n", length(shrunk_files)))

cat(sprintf("Significant associations (lfsr < 0.05) [population=both]:\n"))
for(w in windows){
  f<-file.path(out_dir, sprintf("%s_both_%s_shrunk.txt", out_prefix, w))
  if(!file.exists(f)){next}
  dt<-fread(f)
  if(!"lfsr" %in% names(dt)){next}
  n_pos<-dt[lfsr < 0.05 & !is.na(get(shrunk_col)) & get(shrunk_col) > 1, .N]
  n_neg<-dt[lfsr < 0.05 & !is.na(get(shrunk_col)) & get(shrunk_col) < 1, .N]
  cat(sprintf("%s: %d positives, %d negatives\n", w, n_pos, n_neg))
}
