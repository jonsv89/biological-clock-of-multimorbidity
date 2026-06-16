#### This code has been developed by Jon Sanchez-Valle at BSC ####
# Usage: Rscript 03_compute_RR_contingency.R --disease_a E11 --disease_b I10 --matched_dir ./matched/ --counts_dir ./counts/ --out_dir ./rr_ct/

## Load libraries ##
library(data.table)

## Read arguments ##
args<-commandArgs(trailingOnly = TRUE)
disease_a<-args[which(args == "--disease_a") + 1]
disease_b<-args[which(args == "--disease_b") + 1]
matched_dir<-args[which(args == "--matched_dir") + 1]
counts_dir<-args[which(args == "--counts_dir") + 1]
out_dir<-args[which(args == "--out_dir") + 1]

out_file<-file.path(out_dir, sprintf("rr_ct_%s_%s.rds", disease_a, disease_b))
if(file.exists(out_file)){cat(sprintf("SKIP: already exists %s\n", out_file)); quit(status = 0)}

complete_file<-file.path(matched_dir, sprintf("matched_complete_%s_%s.rds", disease_a, disease_b))
if(!file.exists(complete_file)){cat(sprintf("SKIP: no matched_complete for %s -> %s\n", disease_a, disease_b)); quit(status = 0)}

mc<-readRDS(complete_file)$matched
mc<-mc[!is.na(case_idp) & !is.na(ctrl_idp) & !is.na(sexe)]

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
## RR - Morris & Gardner (1988) ##
## N_A_B = cases with event B (cases_event) ##
## N_A   = total cases in the denominator ##
## N_B   = controls with event B (ctrl_events) ##
## N_0   = total controls in the denonimator ##
## Return NA if no events in cases or controls ##
rr_contingency<-function(N_A_B, N_A, N_B, N_0){
  if(N_A_B == 0 | N_B == 0 | N_A == 0 | N_0 == 0){return(list(RR = NA_real_, CI_low_RR = NA_real_, CI_high_RR = NA_real_, log_RR = NA_real_, SE_RR = NA_real_, p_RR = NA_real_))}
  RR<-(N_A_B/N_A)/(N_B/N_0)
  log_RR<-log(RR)
  SE_RR<-sqrt(1/N_A_B-1/N_A+1/N_B-1/N_0)
  list(RR = RR, CI_low_RR = exp(log_RR - 1.96*SE_RR), CI_high_RR = exp(log_RR + 1.96*SE_RR), log_RR = log_RR, SE_RR = SE_RR, p_RR = 2*pnorm(-abs(log_RR/SE_RR)))
}

## Odds ratio with continuity correction (useful for pairs with zero events) ##
or_contingency<-function(a, b, c, d, cc = 0.5){
  a<-a + cc; b<-b + cc; c<-c + cc; d<-d + cc
  OR<-(a*d)/(b*c)
  lOR<-log(OR)
  SE<-sqrt(1/a + 1/b + 1/c + 1/d)
  list(OR = OR, CI_low = exp(lOR - 1.96*SE), CI_high = exp(lOR + 1.96*SE), log_OR = lOR, SE_OR = SE, p_OR = 2*pnorm(-abs(lOR/SE)))
}

get_abcd<-function(data, t_start, t_end, wtype){
  days_start<-t_start*365.25
  days_end<-t_end*365.25
  if(wtype == "cumulative"){
    cases_u<-data[, .(index_date = first(case_index_date), fu_end = first(case_followup_end), dat_B = first(case_dat_B)), by = case_idp]
    has_fu<-as.numeric(cases_u$fu_end - cases_u$index_date) >= days_end
    event<-!is.na(cases_u$dat_B) & as.numeric(cases_u$dat_B - cases_u$index_date) > 0 & as.numeric(cases_u$dat_B - cases_u$index_date) <= days_end
    a<-sum(event & has_fu, na.rm = TRUE)
    b<-sum(!event & has_fu, na.rm = TRUE)
    has_fu_c<-as.numeric(data$ctrl_followup_end - data$ctrl_index_date) >= days_end
    event_c<-!is.na(data$ctrl_dat_B) & as.numeric(data$ctrl_dat_B - data$ctrl_index_date) > 0 & as.numeric(data$ctrl_dat_B - data$ctrl_index_date) <= days_end
    c<-sum(event_c & has_fu_c, na.rm = TRUE)
    d<-sum(!event_c & has_fu_c, na.rm = TRUE)
  }else{
    cases_u<-data[, .(index_date = first(case_index_date), fu_end = first(case_followup_end), dat_B = first(case_dat_B)), by = case_idp]
    had_B<-!is.na(cases_u$dat_B) & as.numeric(cases_u$dat_B - cases_u$index_date) <= days_start
    has_fu<-as.numeric(cases_u$fu_end - cases_u$index_date) >= days_end
    in_denom<-!had_B & has_fu
    event <-!is.na(cases_u$dat_B) & as.numeric(cases_u$dat_B - cases_u$index_date) > days_start & as.numeric(cases_u$dat_B - cases_u$index_date) <= days_end
    a<-sum(event & in_denom, na.rm = TRUE)
    b<-sum(!event & in_denom, na.rm = TRUE)
    had_B_c<-!is.na(data$ctrl_dat_B) & as.numeric(data$ctrl_dat_B - data$ctrl_index_date) <= days_start
    has_fu_c<-as.numeric(data$ctrl_followup_end - data$ctrl_index_date) >= days_end
    in_denom_c<-!had_B_c & has_fu_c
    event_c<-!is.na(data$ctrl_dat_B) & as.numeric(data$ctrl_dat_B - data$ctrl_index_date) > days_start & as.numeric(data$ctrl_dat_B - data$ctrl_index_date) <= days_end
    c<-sum(event_c & in_denom_c,  na.rm = TRUE)
    d<-sum(!event_c & in_denom_c, na.rm = TRUE)
  }
  list(a = a, b = b, c = c, d = d)
}

## Calculate RR and OR ##
populations<-list(total = list(label = "total", filter = NULL), mujeres = list(label = "mujeres", filter = "D"), hombres = list(label = "hombres", filter = "H"))
results<-list()
idx <-1L
for(wname in names(all_windows)){
  w<-all_windows[[wname]]
  for(pop in populations){
    mc_pop<-if (!is.null(pop$filter)) mc[sexe == pop$filter] else mc
    abcd<-get_abcd(mc_pop, w$t_start, w$t_end, w$type)
    rr<-rr_contingency(abcd$a, abcd$a + abcd$b, abcd$c, abcd$c + abcd$d)
    or<-or_contingency(abcd$a, abcd$b, abcd$c, abcd$d)
    results[[idx]]<-data.table(
      disease_a = disease_a, disease_b = disease_b, window = wname, window_type = w$type, t_start = w$t_start, t_end = w$t_end, population = pop$label, method = "contingency",
      ## RR ##
      RR = rr$RR, CI_low_RR = rr$CI_low_RR, CI_high_RR = rr$CI_high_RR, log_RR = rr$log_RR, SE_RR = rr$SE_RR, p_RR = rr$p_RR,
      ## OR ##
      OR = or$OR, CI_low_OR = or$CI_low, CI_high_OR = or$CI_high, log_OR = or$log_OR, SE_OR = or$SE_OR, p_OR = or$p_OR,
      ## Counts ##
      n_cases = abcd$a + abcd$b, n_ctrl = abcd$c + abcd$d, cases_event = abcd$a, ctrl_events = abcd$c
    )
    idx<-idx + 1L
  }
}

rr_results<-rbindlist(results)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(rr_results, out_file)
cat(sprintf("Saved: %s (%d rows)\n", out_file, nrow(rr_results)))
