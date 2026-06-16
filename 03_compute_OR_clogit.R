#### This code has been developed by Jon Sanchez-Valle at BSC ####
# Usage: Rscript 03_compute_OR_clogit.R --disease_a E11 --disease_b I10 --matched_dir ./matched/ --out_dir ./or_clogit/

## Load libraries ##
library(data.table)
library(survival)

## Read arguments ##
args<-commandArgs(trailingOnly = TRUE)
disease_a<-args[which(args == "--disease_a") + 1]
disease_b<-args[which(args == "--disease_b") + 1]
matched_dir<-args[which(args == "--matched_dir") + 1]
out_dir<-args[which(args == "--out_dir") + 1]

out_file<-file.path(out_dir, sprintf("or_cl_%s_%s.rds", disease_a, disease_b))
if(file.exists(out_file)){cat(sprintf("SKIP: already exists %s\n", out_file)); quit(status = 0)}

all_file<-file.path(matched_dir, sprintf("matched_all_%s_%s.rds", disease_a, disease_b))
if(!file.exists(all_file)){cat(sprintf("SKIP: doest not exist matched_all for %s -> %s\n", disease_a, disease_b)); quit(status = 0)}

ma<-readRDS(all_file)$matched
ma<-ma[!is.na(case_idp) & !is.na(ctrl_idp) & !is.na(sexe)]

## Windows ##
windows_cumul<-list(
  "0_1" = list(type = "cumulative", t_start = 0, t_end = 1),
  "0_2" = list(type = "cumulative", t_start = 0, t_end = 2),
  "0_3" = list(type = "cumulative", t_start = 0, t_end = 3),
  "0_4" = list(type = "cumulative", t_start = 0, t_end = 4),
  "0_5" = list(type = "cumulative", t_start = 0, t_end = 5)
)
windows_cont<-list(
  "1_2" = list(type = "continuous_conditional", t_start = 1, t_end = 2),
  "2_3" = list(type = "continuous_conditional", t_start = 2, t_end = 3),
  "3_4" = list(type = "continuous_conditional", t_start = 3, t_end = 4),
  "4_5" = list(type = "continuous_conditional", t_start = 4, t_end = 5)
)
all_windows<-c(windows_cumul, windows_cont)

## Functions ##
make_long<-function(data, t_start, t_end, wtype) {
  days_start<-t_start*365.25
  days_end  <-t_end*365.25
  if(wtype == "cumulative"){
    cases_u<-data[, .(index_date = first(case_index_date), fu_end = first(case_followup_end), dat_B = first(case_dat_B), sexe = first(sexe), pair_id = first(case_idp), exposure = 1L), by = case_idp]
    has_fu<-as.numeric(cases_u$fu_end - cases_u$index_date) >= days_end
    cases_u[, event := as.integer(!is.na(dat_B) & as.numeric(dat_B - index_date) > 0 & as.numeric(dat_B - index_date) <= days_end)]
    cases_u<-cases_u[has_fu]
    ctrls<-data[case_idp %in% cases_u$case_idp, .(index_date = ctrl_index_date, fu_end = ctrl_followup_end, dat_B = ctrl_dat_B, sexe = sexe, pair_id = case_idp, exposure = 0L)]
    has_fu_c<-as.numeric(ctrls$fu_end - ctrls$index_date) >= days_end
    ctrls[, event := as.integer(!is.na(dat_B) & as.numeric(dat_B - index_date) > 0 & as.numeric(dat_B - index_date) <= days_end)]
    ctrls<-ctrls[has_fu_c]
    long<-rbind(cases_u[, .(index_date, sexe, pair_id, exposure, event)], ctrls[, .(index_date, sexe, pair_id, exposure, event)])
  }else{
    cases_u<-data[, .(index_date = first(case_index_date), fu_end = first(case_followup_end), dat_B = first(case_dat_B), sexe = first(sexe), pair_id = first(case_idp), exposure = 1L), by = case_idp]
    had_B<-!is.na(cases_u$dat_B) & as.numeric(cases_u$dat_B - cases_u$index_date) <= days_start
    has_fu<-as.numeric(cases_u$fu_end - cases_u$index_date) >= days_end
    cases_u<-cases_u[!had_B & has_fu]
    cases_u[, event := as.integer(!is.na(dat_B) & as.numeric(dat_B - index_date) > days_start & as.numeric(dat_B - index_date) <= days_end)]
    ctrls<-data[case_idp %in% cases_u$case_idp,.(index_date = ctrl_index_date, fu_end = ctrl_followup_end, dat_B = ctrl_dat_B, sexe = sexe, pair_id = case_idp, exposure = 0L)]
    had_B_c<-!is.na(ctrls$dat_B) & as.numeric(ctrls$dat_B - ctrls$index_date) <= days_start
    has_fu_c<-as.numeric(ctrls$fu_end - ctrls$index_date) >= days_end
    ctrls<-ctrls[!had_B_c & has_fu_c]
    ctrls[, event := as.integer( !is.na(dat_B) & as.numeric(dat_B - index_date) > days_start & as.numeric(dat_B - index_date) <= days_end)]
    valid_pairs<-intersect(cases_u$case_idp, ctrls[, unique(pair_id)])
    cases_u<-cases_u[case_idp %in% valid_pairs]
    ctrls<-ctrls[pair_id %in% valid_pairs]
    long<-rbind(cases_u[, .(index_date, sexe, pair_id, exposure, event)], ctrls[, .(index_date, sexe, pair_id, exposure, event)])
  }
  long[, age_exact := as.numeric(index_date - as.Date("2008-01-01"))/365.25]
  long[, index_year := as.integer(format(index_date, "%Y"))]
  long
}

fit_clogit<-function(long_data){
  tryCatch({
    strat_counts<-long_data[, .(n_cases = sum(exposure == 1), n_ctrl  = sum(exposure == 0)), by = pair_id]
    valid<-strat_counts[n_cases >= 1 & n_ctrl >= 1, pair_id]
    ld<-long_data[pair_id %in% valid]
    if(nrow(ld) == 0 || ld[, uniqueN(pair_id)] < 10) return(NULL)
    fit<-clogit(event ~ exposure + age_exact + index_year + strata(pair_id), data = ld, method = "efron")
    coef_exp<-coef(fit)["exposure"]
    se_exp<-sqrt(diag(vcov(fit)))["exposure"]
    p_exp<-summary(fit)$coefficients["exposure", "Pr(>|z|)"]
    list(OR = exp(coef_exp), CI_low = exp(coef_exp - 1.96 * se_exp), CI_high = exp(coef_exp + 1.96 * se_exp), log_OR = coef_exp, SE = se_exp, p = p_exp, n_pairs = ld[, uniqueN(pair_id)])
  }, error = function(e) NULL)
}

## Calculate OR ##
populations<-list(total = list(label = "total", filter = NULL), mujeres = list(label = "mujeres", filter = "D"), hombres = list(label = "hombres", filter = "H"))
results<-list()
idx <-1L
cat(sprintf("Calculating clogit for %s -> %s...\n", disease_a, disease_b))

for(wname in names(all_windows)){
  w<-all_windows[[wname]]
  for(pop in populations){
    ma_pop<-if(!is.null(pop$filter)) ma[sexe == pop$filter] else ma
    long<-make_long(ma_pop, w$t_start, w$t_end, w$type)
    or_cl<-fit_clogit(long)
    rm(long); gc()

    results[[idx]]<-data.table(
      disease_a = disease_a, disease_b = disease_b, window = wname, window_type = w$type, t_start = w$t_start, t_end = w$t_end, population = pop$label, method = "clogit",
      OR = if (!is.null(or_cl)) or_cl$OR else NA_real_,
      CI_low = if (!is.null(or_cl)) or_cl$CI_low else NA_real_,
      CI_high = if (!is.null(or_cl)) or_cl$CI_high else NA_real_,
      log_OR = if (!is.null(or_cl)) or_cl$log_OR else NA_real_,
      SE = if (!is.null(or_cl)) or_cl$SE else NA_real_,
      p  = if (!is.null(or_cl)) or_cl$p else NA_real_,
      n_cases = if (!is.null(or_cl)) or_cl$n_pairs else NA_integer_,
      n_ctrl = NA_integer_,
      cases_event = NA_integer_,
      ctrl_events = NA_integer_
    )
    idx<-idx + 1L
  }
}
or_results<-rbindlist(results)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(or_results, out_file)
cat(sprintf("Saved: %s (%d rows)\n", out_file, nrow(or_results)))