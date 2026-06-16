## This code has been developed by Jon Sanchez-Valle at the Barcelona Supercomputing Center ##
## The main focus of the code is to study comorbidities and Catalonia and compare them to the ones in Denmark ##

## load libraries ##
suppressPackageStartupMessages({
  library("data.table")
  require(dplyr)
  library("igraph")
  library("UpSetR")
  library("dendextend")
  library("ggplot2")
  library("gplots")
  library("forestplot")
  library("EbayesThresh")
  library("ashr")
  library("VennDiagram")
  library("gridExtra")
  library("grid")
  library("plotly")
  library("ggalluvial")
  library("MASS")
  library("patchwork")
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(scales)
  library(ggnewscale)
  library(ggrepel)
  library(flextable)
  library(officer)
})

## Create needed directories ##
if("ManuscriptFiles"%in%list.files()==FALSE){dir.create("ManuscriptFiles")}
if("Plots"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/Plots")}
if("Results"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/Results")}
if("IntermediateFiles"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/IntermediateFiles")}
if("Disease_prevalence_number_comorbidities_plot"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/Disease_prevalence_number_comorbidities_plot")}
if("Data"%in%list.files()==FALSE){dir.create("Data")}
if("Tables_for_plots"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/Tables_for_plots")}


## Provide the arguments ##
args<-commandArgs(trailingOnly = TRUE)

## Load functions ##
overlapmetrics<-function(A, B) {
  inter<-intersect(A, B)
  union<-union(A, B)
  data.frame(jaccard = length(inter) / length(union), C_in_D = length(inter) / length(A), D_in_C = length(inter) / length(B))
}

#### Prepare comorbidity networks: remove dagger-asterisk and compute directionality ####
## This must be run BEFORE any analysis section ##
## Follows Westergaard et al. 2019 and Jensen et al. 2014 directionality method ##
if(args[1]=="prepare_networks"){
  source("icd10_dagger_asterisk.R")  # loads ASTERISK_CODES_3DIG and filter_dagger_asterisk()
  if("networks"%in%list.files("Results")==FALSE){dir.create("Results/networks")}
  windows<-c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
  pops<-c("both","women","men")
  for(pop in pops){
    for(w in windows){
      in_file <-sprintf("Results/shrinkage_rr_events/RR_contingency_%s_%s_shrunk.txt", pop, w)
      out_file<-sprintf("Results/networks/RR_net_%s_%s.txt", pop, w)
      if(file.exists(out_file)){cat(sprintf("SKIP: ya existe %s\n", basename(out_file))); next}
      if(!file.exists(in_file)){cat(sprintf("SKIP: no existe %s\n", basename(in_file))); next}
      cat(sprintf("Processing %s / %s...", pop, w))
      dt<-fread(in_file)
      cat(sprintf("Total pairs: %d\n", nrow(dt)))
      ## Keep only pairs with at least 1 event in cases or controls
      dt<-dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk) & is.finite(SE_RR_shrunk)]
      ## Remove dagger-asterisk pairs (specific pairs only, following Westergaard)
      dt<-filter_dagger_asterisk(dt, col_a = "disease_a", col_b = "disease_b")
      ## Compute directionality following Jensen et al. 2014 / Westergaard et al. 2019:
      ## For each pair A->B, find its reverse B->A in the same table.
      ## directionality theta(A->B) = N(A->B) / (N(A->B) + N(B->A))
      ## where N = cases_event (number of cases who develop B after A)
      ## Significance: binomial test (p = 0.5, one-sided: A->B more frequent than B->A)
      cat("Calculating directionality...\n")
      ## Create reverse key for efficient lookup
      dt[, pair_key := paste(disease_a, disease_b, sep = "__")]
      dt[, pair_key_rev := paste(disease_b, disease_a, sep = "__")]
      ## Lookup N(B->A) for each A->B pair
      n_forward<-dt[, .(pair_key, n_forward = cases_event)]
      n_reverse<-dt[, .(pair_key_rev = pair_key, n_reverse = cases_event)]
      dt<-merge(dt, n_reverse, by = "pair_key_rev", all.x = TRUE)
      dt[is.na(n_reverse), n_reverse := 0L]
      ## If N(A->B) + N(B->A) == 0, theta = NA
      dt[, n_total_dir := cases_event + n_reverse]
      dt[, theta := fifelse(n_total_dir > 0, cases_event/n_total_dir, NA_real_)]
      ## Binomial test for preferred direction (one-sided: theta > 0.5)
      ## i.e. A->B occurs more than B->A
      dt[, directionality_pval := {
        mapply(function(n_ab, n_tot){
          if(is.na(n_tot) || n_tot == 0){return(NA_real_)}
          binom.test(n_ab, n_tot, p = 0.5, alternative = "greater")$p.value
        }, cases_event, n_total_dir)
      }]
      ## FDR correction for directionality within this window x pop
      dt[!is.na(directionality_pval), directionality_fdr := p.adjust(directionality_pval, method = "BH")]
      dt[is.na(directionality_pval), directionality_fdr := NA_real_]
      ## Preferred direction: theta > 0.5 & FDR < 0.05
      dt[, preferred_direction := !is.na(directionality_fdr) & directionality_fdr < 0.05 & !is.na(theta) & theta > 0.5]
      ## Clean up temporary columns
      dt[, c("pair_key","pair_key_rev","n_reverse","n_total_dir") := NULL]
      fwrite(dt, out_file, sep = "\t", quote = FALSE)
      cat(sprintf("Saved: %s (%d pares)\n", basename(out_file), nrow(dt)))
    }
  }
  cat("Networks available at: Results/networks/")
}

#### Analyse gender-related differences on disease prevalence and age of diagnoses, including comparison with Denmark ####
## Plot gender-associated disease prevalence differences in Catalonia and Denmark and look for differences between them ##
if(args[1]=="compare_prevalences"){
  #### Catalonia ####
  ## @ @ @@ @@ @ @ ##
  dt<-fread("Diagnoses_20080101_20181231_first_diagnoses_of_each_disease.txt",stringsAsFactors = F,sep="|")
  prevalences<-dcast(dt,cod ~ sexe,value.var = "idp",fun.aggregate = uniqueN)
  prevalences<-cbind(prevalences[[1]],prevalences[[2]]+prevalences[[3]],prevalences[,2:3],3096908,2960950)
  colnames(prevalences)<-c("disease","both","women","men","nwomen","nmen")
  write.table(prevalences,"Data/ICD10_prevalence_Catalonia.txt",quote=F,sep="\t",row.names=F)
  ## Convert disease diagnoses into disease categories ##
  code<-c(paste("A0",0:9,sep=""),paste("A",10:99,sep=""),paste("B0",0:9,sep=""),paste("B",10:99,sep=""),paste("C0",0:9,sep=""),paste("C",10:99,sep=""),
          paste("D0",0:9,sep=""),paste("D",10:48,sep=""),paste("D",50:89,sep=""),paste("E0",0:9,sep=""),paste("E",10:99,sep=""),paste("F0",0:9,sep=""),
          paste("F",10:99,sep=""),paste("G0",0:9,sep=""),paste("G",10:99,sep=""),paste("H0",0:9,sep=""),paste("H",10:59,sep=""),paste("H",60:95,sep=""),
          paste("I0",0:9,sep=""),paste("I",10:99,sep=""),paste("J0",0:9,sep=""),paste("J",10:99,sep=""),paste("K0",0:9,sep=""),paste("K",10:93,sep=""),
          paste("L0",0:9,sep=""),paste("L",10:99,sep=""),paste("M0",0:9,sep=""),paste("M",10:99,sep=""),paste("N0",0:9,sep=""),paste("N",10:99,sep=""),
          paste("O0",0:9,sep=""),paste("O",10:99,sep=""),paste("P0",0:9,sep=""),paste("P",10:96,sep=""),paste("Q0",0:9,sep=""),paste("Q",10:99,sep=""),
          paste("R0",0:9,sep=""),paste("R",10:99,sep=""),paste("S0",0:9,sep=""),paste("S",10:99,sep=""),paste("T0",0:9,sep=""),paste("T",10:98,sep=""))
  cate<-c(rep("I",200),rep("II",149),rep("III",40),rep("IV",100),rep("V",100),rep("VI",100),rep("VII",60),rep("VIII",36),rep("IX",100),rep("X",100),
          rep("XI",94),rep("XII",100),rep("XIII",100),rep("XIV",100),rep("XV",100),rep("XVI",97),rep("XVII",100),rep("XVIII",100),rep("XIX",199))
  catname<-c("Infectious and parasitic","Neoplasms","Blood and blood-forming organs (immune)",
             "Endocrine, nutritional and metabolic","Mental and behavioural",
             "Nervous system","Eye and adnexa","Ear and mastoid process",
             "Circulatory system","Respiratory system","Digestive system",
             "Skin and subcutaneous tissue","Musculoskeletal system and connective tissue",
             "Genitourinary system","Pregnancy, childbirth and the puerperium",
             "Certain conditions originating in the perinatal period","Congenital malformations and chromosomal abnormalities",
             "Symptoms, signs and abnormal laboratory findings","Injury, poisoning")
  names(catname)<-unique(cate)
  names(cate)<-code
  ## Disease category ##
  catename<-as.character(catname[cate[code]])
  names(catename)<-code
  colcod<-c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C","#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54","#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
  names(colcod)<-unique(cate)
  catcol<-as.character(colcod) ; names(catcol)<-catname
  codcol<-names(colcod) ; names(codcol)<-as.character(colcod)
  distocol<-as.character(colcod[as.character(cate)]) ; names(distocol)<-names(cate)
  if(args[2]=="agedistribution"){
    dt[, diseasecategory := as.character(catename[cod])]
    dt[, color := as.character(distocol[cod])]
    ## First diagnosis per patient, disease and sex ##
    dtf<-dt[!is.na(diseasecategory), .SD[1], by = .(idp, cod, sexe, diseasecategory, color)]
    ## Mean age at first diagnosis by gender ##
    dx_mean_age<-dtf[, .(mean_age = mean(edad, na.rm = TRUE)),by = .(cod, sexe, diseasecategory, color)]
    rm(dtf); gc()
    if(length(which(is.na(dx_mean_age$diseasecategory)))>0){dx_mean_age<-dx_mean_age[-which(is.na(dx_mean_age$diseasecategory))]}
    ## Color map ##
    col_map<-unique(dx_mean_age[, .(diseasecategory, color)])
    fill_vals<-setNames(col_map$color, col_map$diseasecategory)
    ## Labels ## 
    sex_labs<-c("D" = "Women", "H" = "Men")
    ## Plot all diseases ##
    p<-ggplot(dx_mean_age,aes(x = mean_age, fill = diseasecategory)) +
      geom_density(position = "stack",alpha = 0.95,linewidth = 0.2,adjust = 1) +
      facet_grid(sexe ~ ., labeller = as_labeller(sex_labs)) +
      scale_fill_manual(values = fill_vals, name = "Disease Category") +
      scale_colour_manual(values = fill_vals, guide = "none") +
      labs(x = "Mean age at first diagnosis",y = "Density") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),legend.position = "right")
    ggsave("ManuscriptFiles/Plots/MeanAgeDiagnosisAll_Catalonia.pdf", p, width = 12, height = 8, useDingbats = FALSE)
    
    ## Plot removing Pregnancy ##
    dx_mean_age<-dx_mean_age[-which(dx_mean_age$diseasecategory=="Pregnancy, childbirth and the puerperium")]
    ## Color map ##
    col_map<-unique(dx_mean_age[, .(diseasecategory, color)])
    fill_vals<-setNames(col_map$color, col_map$diseasecategory)
    p<-ggplot(dx_mean_age,aes(x = mean_age, fill = diseasecategory)) +
      geom_density(position = "stack",alpha = 0.95,linewidth = 0.2,adjust = 1) +
      facet_grid(sexe ~ ., labeller = as_labeller(sex_labs)) +
      scale_fill_manual(values = fill_vals, name = "Disease Category") +
      scale_colour_manual(values = fill_vals, guide = "none") +
      labs(x = "Mean age at first diagnosis",y = "Density") +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),legend.position = "right")
    ggsave("ManuscriptFiles/Plots/MeanAgeDiagnosisAllExceptPregnancy_Catalonia.pdf", p, width = 12, height = 8, useDingbats = FALSE)
   }
  
  prevalences<-cbind(prevalences,as.character(catname[cate[prevalences$disease]]),as.character(distocol[prevalences$disease]))
  colnames(prevalences)[7:8]<-c("diseasecategory","color")
  if(length(which(is.na(prevalences$diseasecategory)))>0){prevalences<-prevalences[-which(is.na(prevalences$diseasecategory))]}
  prevalences<-prevalences[,c(1,7,3:6,8)]
  ## Look for differences in gender prevalence disease by disease ##
  fisher_sex_test<-function(n_w, n_m, N_w, N_m){
    mat<-matrix(c(n_w, N_w - n_w,n_m, N_m - n_m),nrow = 2,byrow = TRUE)
    test<-fisher.test(mat)
    return(c(p_value = test$p.value,OR = unname(test$estimate)))
  }
  res<-t(mapply(fisher_sex_test,prevalences$women,prevalences$men,prevalences$nwomen,prevalences$nmen))
  prevalences$p_value<-res[, "p_value"]
  prevalences$OR<-res[, "OR"]
  # Correct for multiple testing
  prevalences$p_adj<-p.adjust(prevalences$p_value, method = "BH")
  ## Remove pregnancy related diseases ##
  prevalences<-prevalences[-which(prevalences$diseasecategory=="Pregnancy, childbirth and the puerperium")]
  length(intersect(which(prevalences$p_adj<=0.05),which(prevalences$OR>=1.3)))
  length(intersect(which(prevalences$p_adj<=0.05),which(prevalences$OR<=(1/1.3))))
  # Define thresholds
  or_threshold<-1.3
  p_threshold <-0.05
  prevalences[, sex_bias := fifelse(
    p_adj <= p_threshold & OR >= or_threshold, "women",
    fifelse(p_adj <= p_threshold & OR <= 1/or_threshold, "men", "none")
  )]
  write.table(prevalences,"ManuscriptFiles/Results/SupplementaryData_gender_related_diseaseprevalence_differences_Catalonia.txt",quote=F,sep="\t",row.names=F)
  
  #### Gender bias by disease ####
  prevalences2<-prevalences[-which(prevalences$sex_bias=="none")]
  prevalences2[, total_cases := women + men]
  prevalences2$OR<-log(prevalences2$OR)
  prevalences2$total_cases<-log(prevalences2$total_cases)
  prevalences2$OR[which(prevalences2$OR==Inf)]<-11
  prevalences2$OR[which(prevalences2$OR==-Inf)]<--11
  ## Check if OR and total cases are valid ## 
  prevalences_plot<-prevalences2[!is.na(OR) & !is.na(total_cases)]
  ## Stablish the order of the disease categories ##
  disease_order<-unique(prevalences_plot$diseasecategory)
  prevalences_plot[, diseasecategory := factor(diseasecategory, levels = disease_order)]
  ## Plot the association
  pdf("ManuscriptFiles/Plots/Gender_bias_by_disease.pdf",width = 10,height = 10)
  ggplot( prevalences_plot, aes(x = total_cases,y = OR,color = color)) +
    geom_point(size = 3,alpha = 0.6) +
    scale_color_identity(name = "disease category", guide = "legend", breaks = prevalences_plot$color[match(levels(prevalences_plot$diseasecategory), prevalences_plot$diseasecategory)], labels = levels(prevalences_plot$diseasecategory)) +
    xlab("Log (total cases, women + men)") +
    ylab("Log (OR)") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  dev.off()
  
  #### Make the plot interactive ####
  dig3dis<-fread("Data/ICD10_three_digits_names.txt",stringsAsFactors = F)
  tdig3dis<-dig3dis$Name ; names(tdig3dis)<-dig3dis$Code
  prevalences2<-cbind(prevalences2,as.character(tdig3dis[prevalences2$disease]))
  colnames(prevalences2)[13]<-c("diseasename")
  prevalences_plot<-prevalences2[!is.na(OR) & !is.na(total_cases)]
  disease_order<-unique(prevalences_plot$diseasecategory)
  prevalences_plot[, diseasecategory := factor(diseasecategory, levels = disease_order)]
  p<-ggplot(prevalences_plot,
    aes(x = total_cases, y = OR, color = color,
      text = paste0("Disease name: ", diseasename, "<br>", "Disease: ", disease, "<br>", "Category: ", diseasecategory, "<br>", "OR: ", round(exp(OR), 0), "<br>", "Women: ", women,"<br>", "Men: ", men,"<br>"))) +
    geom_point(size = 3, alpha = 0.6) +
    scale_color_identity(name = "disease category", guide = "legend", breaks = prevalences_plot$color[match(levels(prevalences_plot$diseasecategory), prevalences_plot$diseasecategory)], labels = levels(prevalences_plot$diseasecategory)) +
    xlab("Total cases (women + men)") +
    ylab("Odds Ratio (OR)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  ## Make it interactive ## 
  ggplotly(p, tooltip = "text")
  
  #### Enrichment analysis ####
  ## Number of analyzed diseases ##
  n_all<-length(prevalences[, unique(disease)])
  enrichment_by_category<-function(dt, category_col, bias_value){
    biased<-dt[sex_bias == bias_value, unique(disease)]
    if(length(biased) == 0){return(NULL)}
    results<-dt[, {
      diseases_in_cat<-unique(disease); a<-sum(diseases_in_cat %in% biased); b<-length(biased) - a
      c<-length(diseases_in_cat) - a; d<-n_all - a - b - c
      mat<-matrix(c(a, b, c, d), nrow = 2,dimnames = list(Bias = c("InCategory", "OutCategory"),Status = c("Biased", "NotBiased")))
      ft<-fisher.test(mat)
      list(n_biased_in_category = a,n_diseases_in_category = length(diseases_in_cat),OR_enrichment = unname(ft$estimate),p_value = ft$p.value)
    }, by = category_col]
    results[, sex_bias := bias_value]
    return(results)
  }
  # Enrichment for each gender separately ##
  enrichment_results<-rbindlist(list(enrichment_by_category(prevalences,category_col = "diseasecategory",bias_value = "women"), enrichment_by_category(prevalences,category_col = "diseasecategory",bias_value = "men")), use.names = TRUE)
  # Adjust for multiple testing by gender ##
  enrichment_results[, p_adj := p.adjust(p_value, method = "BH"),by = sex_bias]
  ## Plot gender differences
  plot_data<-enrichment_results[, .(diseasecategory, sex_bias, OR_enrichment, log_OR = log(OR_enrichment), n_biased_in_category, p_adj)]
  ## Add color palette ##
  plot_data<-cbind(plot_data,as.character(catcol[plot_data$diseasecategory]))
  colnames(plot_data)[7]<-"color"
  ## Find significant enrichments ##
  plot_data[, significant := p_adj <= 0.05]
  plot_data$significant[which(plot_data$significant==TRUE)]<-"Yes"
  plot_data$significant[which(plot_data$significant==FALSE)]<-"No"
  plot_data[, diseasecategory := factor(diseasecategory,levels = rev(unique(diseasecategory)))]
  categories<-plot_data %>%
    distinct(diseasecategory, .keep_all = TRUE) %>%
    mutate(disease_label = paste0("<span style='color:", color, "'>\u25CF</span> ", diseasecategory))
  categories<-categories %>%
    mutate(disease_label = factor(as.character(disease_label), levels = disease_label))
  
  pdf("ManuscriptFiles/Plots/Enrichment_gender_bias_period_prevalence_by_disease_categories.pdf",width = 10,height = 10)
    ggplot(plot_data, aes(x = log_OR, y = diseasecategory, color = sex_bias, size = n_biased_in_category, alpha = significant)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_alpha_manual(values = c("Yes" = 1, "No" = 0.3)) +
      scale_color_manual(values = c("women" = "#A94C40", "men" = "#5474A5"), labels = c("Men-biased", "Women-biased")) +
      scale_size_continuous(name = "Number of diseases") +
      labs(
        x = "Log odds ratio of enrichment", y = "Disease category (ICD-10)", color = "Gender bias",
        title = "Enrichment of disease categories among diseases\nwith significant gender-biased diagnosis rates"
      ) +
      theme_minimal() +
      theme(legend.position = "right", axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold", hjust=0.5))
  dev.off()
  
  pdf("ManuscriptFiles/Plots/LabelColors.pdf",width = 10,height = 12)
    plot(rep(1,18),1:18,col=plot_data$color[18:1],pch=16,cex=2)
  dev.off()
  
  #### Plot Denmark ####
  denprevalences<-fread("Epidemiology/41467_2019_8475_MOESM4_ESM.txt",stringsAsFactors = F,sep="\t")
  disnames<-denprevalences$Name ; names(disnames)<-denprevalences$Code
  denprevalences<-cbind(denprevalences[[1]],denprevalences[[3]]+denprevalences[[4]],denprevalences[,c(4,3)])
  denprevalences<-cbind(denprevalences,3330464,3579212)
  colnames(denprevalences)<-c("disease","both","women","men","nwomen","nmen")
  ## Convert disease diagnoses into disease categories ##
  code<-c(paste("A0",0:9,sep=""),paste("A",10:99,sep=""),paste("B0",0:9,sep=""),paste("B",10:99,sep=""),paste("C0",0:9,sep=""),paste("C",10:99,sep=""),
          paste("D0",0:9,sep=""),paste("D",10:48,sep=""),paste("D",50:89,sep=""),paste("E0",0:9,sep=""),paste("E",10:99,sep=""),paste("F0",0:9,sep=""),
          paste("F",10:99,sep=""),paste("G0",0:9,sep=""),paste("G",10:99,sep=""),paste("H0",0:9,sep=""),paste("H",10:59,sep=""),paste("H",60:95,sep=""),
          paste("I0",0:9,sep=""),paste("I",10:99,sep=""),paste("J0",0:9,sep=""),paste("J",10:99,sep=""),paste("K0",0:9,sep=""),paste("K",10:93,sep=""),
          paste("L0",0:9,sep=""),paste("L",10:99,sep=""),paste("M0",0:9,sep=""),paste("M",10:99,sep=""),paste("N0",0:9,sep=""),paste("N",10:99,sep=""),
          paste("O0",0:9,sep=""),paste("O",10:99,sep=""),paste("P0",0:9,sep=""),paste("P",10:96,sep=""),paste("Q0",0:9,sep=""),paste("Q",10:99,sep=""),
          paste("R0",0:9,sep=""),paste("R",10:99,sep=""),paste("S0",0:9,sep=""),paste("S",10:99,sep=""),paste("T0",0:9,sep=""),paste("T",10:98,sep=""))
  cate<-c(rep("I",200),rep("II",149),rep("III",40),rep("IV",100),rep("V",100),rep("VI",100),rep("VII",60),rep("VIII",36),rep("IX",100),rep("X",100),
          rep("XI",94),rep("XII",100),rep("XIII",100),rep("XIV",100),rep("XV",100),rep("XVI",97),rep("XVII",100),rep("XVIII",100),rep("XIX",199))
  catname<-c("Infectious and parasitic","Neoplasms","Blood and blood-forming organs (immune)",
             "Endocrine, nutritional and metabolic","Mental and behavioural",
             "Nervous system","Eye and adnexa","Ear and mastoid process",
             "Circulatory system","Respiratory system","Digestive system",
             "Skin and subcutaneous tissue","Musculoskeletal system and connective tissue",
             "Genitourinary system","Pregnancy, childbirth and the puerperium",
             "Certain conditions originating in the perinatal period","Congenital malformations and chromosomal abnormalities",
             "Symptoms, signs and abnormal laboratory findings","Injury, poisoning")
  names(catname)<-unique(cate)
  names(cate)<-code
  ## Disease category ##
  catename<-as.character(catname[cate[code]])
  names(catename)<-code
  colcod<-c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C","#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54","#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
  names(colcod)<-unique(cate)
  catcol<-as.character(colcod) ; names(catcol)<-catname
  codcol<-names(colcod) ; names(codcol)<-as.character(colcod)
  distocol<-as.character(colcod[as.character(cate)]) ; names(distocol)<-names(cate)
  denprevalences<-cbind(denprevalences,as.character(catname[cate[denprevalences$disease]]),as.character(distocol[denprevalences$disease]))
  colnames(denprevalences)[7:8]<-c("diseasecategory","color")
  if(length(which(is.na(denprevalences$diseasecategory)))>0){denprevalences<-denprevalences[-which(is.na(denprevalences$diseasecategory))]}
  denprevalences<-denprevalences[,c(1,7,3:6,8)]
  ## Look for differences in gender prevalence disease by disease ##
  fisher_sex_test<-function(n_w, n_m, N_w, N_m){
    mat<-matrix(c(n_w, N_w - n_w,n_m, N_m - n_m),nrow = 2,byrow = TRUE)
    test<-fisher.test(mat)
    return(c(p_value = test$p.value,OR = unname(test$estimate)))
  }
  res<-t(mapply(fisher_sex_test,denprevalences$women,denprevalences$men,denprevalences$nwomen,denprevalences$nmen))
  denprevalences$p_value<-res[, "p_value"]
  denprevalences$OR<-res[, "OR"]
  # Correct for multiple testing
  denprevalences$p_adj<-p.adjust(denprevalences$p_value, method = "BH")
  write.table(denprevalences,"ManuscriptFiles/Results/SupplementaryData_gender_related_diseaseprevalence_differences_Denmark.txt",quote=F,sep="\t",row.names=F)
  length(intersect(which(denprevalences$p_adj<=0.05),which(denprevalences$OR>=1.3)))
  length(intersect(which(denprevalences$p_adj<=0.05),which(denprevalences$OR<=(1/1.3))))
  # Define thresholds
  or_threshold<-1.3
  p_threshold <-0.05
  denprevalences[, sex_bias := fifelse(
    p_adj <= p_threshold & OR >= or_threshold, "women",
    fifelse(p_adj <= p_threshold & OR <= 1/or_threshold, "men", "none")
  )]
  #### Gender bias by disease ####
  denprevalences2<-denprevalences[-which(denprevalences$sex_bias=="none")]
  denprevalences2[, total_cases := women + men]
  denprevalences2$OR<-log(denprevalences2$OR)
  denprevalences2$total_cases<-log(denprevalences2$total_cases)
  ## Check if OR and total cases are valid ## 
  denprevalences_plot<-denprevalences2[!is.na(OR) & !is.na(total_cases)]
  ## Stablish the order of the disease categories ##
  disease_order<-unique(denprevalences_plot$diseasecategory)
  denprevalences_plot[, diseasecategory := factor(diseasecategory, levels = disease_order)]
  ## Plot the association
  pdf("ManuscriptFiles/Plots/Gender_bias_by_disease_Denmark.pdf",width = 10,height = 10)
  ggplot(denprevalences_plot,
    aes(x = total_cases,y = OR,color = color)) +
    geom_point(size = 3,alpha = 0.6) +
    scale_color_identity(name = "disease category", guide = "legend", breaks = denprevalences_plot$color[match(levels(denprevalences_plot$diseasecategory),denprevalences_plot$diseasecategory)], labels = levels(denprevalences_plot$diseasecategory)) +
    xlab("Log (total cases, women + men)") +
    ylab("Log (OR)") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = "right")
  dev.off()
  #### Enrichment analysis ####
  ## Number of analyzed diseases ##
  n_all<-length(denprevalences[, unique(disease)])
  enrichment_by_category<-function(dt, category_col, bias_value) {
    biased<-dt[sex_bias == bias_value, unique(disease)]
    if(length(biased) == 0){return(NULL)}
    results<-dt[, {
      diseases_in_cat<-unique(disease)
      a<-sum(diseases_in_cat %in% biased); b<-length(biased) - a; c<-length(diseases_in_cat) - a; d<-n_all - a - b - c
      mat<-matrix(c(a, b, c, d), nrow = 2,dimnames = list(Bias = c("InCategory", "OutCategory"),Status = c("Biased", "NotBiased")))
      ft<-fisher.test(mat)
      list(n_biased_in_category = a,n_diseases_in_category = length(diseases_in_cat),OR_enrichment = unname(ft$estimate),p_value = ft$p.value)
    }, by = category_col]
    results[, sex_bias := bias_value]
    return(results)
  }
  ## Enrichment for each gender separately ##
  denenrichment_results<-rbindlist(list(enrichment_by_category(denprevalences,category_col = "diseasecategory",bias_value = "women"), enrichment_by_category(denprevalences,category_col = "diseasecategory",bias_value = "men")), use.names = TRUE)
  ## Adjust for multiple testing by gender ##
  denenrichment_results[, p_adj := p.adjust(p_value, method = "BH"),by = sex_bias]
  ## Plot gender differences
  plot_data_den<-denenrichment_results[, .(diseasecategory, sex_bias, OR_enrichment, log_OR = log(OR_enrichment), n_biased_in_category, p_adj)]
  ## Add color palette ##
  plot_data_den<-cbind(plot_data_den,as.character(catcol[plot_data_den$diseasecategory]))
  colnames(plot_data_den)[7]<-"color"
  ## Find significant enrichments ##
  plot_data_den[, significant := p_adj <= 0.05]
  plot_data_den$significant[which(plot_data_den$significant==TRUE)]<-"Yes"
  plot_data_den$significant[which(plot_data_den$significant==FALSE)]<-"No"
  plot_data_den[, diseasecategory := factor(diseasecategory,levels = rev(unique(diseasecategory)))]
  categories<-plot_data_den %>%
    distinct(diseasecategory, .keep_all = TRUE) %>%
    mutate(disease_label = paste0("<span style='color:", color, "'>\u25CF</span> ", diseasecategory))
  categories<-categories %>%
    mutate(disease_label = factor(as.character(disease_label), levels = disease_label))
  pdf("ManuscriptFiles/Plots/Enrichment_gender_bias_period_prevalence_by_disease_categories_Denmark.pdf",width = 10,height = 10)
  ggplot(plot_data_den, aes(x = log_OR, y = diseasecategory, color = sex_bias, size = n_biased_in_category, alpha = significant)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_alpha_manual(values = c("Yes" = 1, "No" = 0.3)) +
    scale_color_manual(values = c("women" = "#A94C40", "men" = "#5474A5"), labels = c("Men-biased", "Women-biased")) +
    scale_size_continuous(name = "Number of diseases") +
    labs(x = "Log odds ratio of enrichment", y = "Disease category (ICD-10)", color = "Gender bias", title = "Enrichment of disease categories among diseases\nwith significant gender-biased diagnosis rates") +
    theme_minimal() +
    theme(legend.position = "right", axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold", hjust=0.5))
  dev.off()
  pdf("ManuscriptFiles/Plots/LabelColors_Denmark.pdf",width = 10,height = 12)
    plot(rep(1,18),1:18,col=plot_data_den$color[18:1],pch=16,cex=2)
  dev.off()
  
  #### Compare with Denmark ####
  commondis<-intersect(prevalences$disease,denprevalences$disease)
  denwom<-intersect(commondis,denprevalences$disease[which(denprevalences$sex_bias=="women")])
  catwom<-intersect(commondis,prevalences$disease[which(prevalences$sex_bias=="women")])
  denmen<-intersect(commondis,denprevalences$disease[which(denprevalences$sex_bias=="men")])
  catmen<-intersect(commondis,prevalences$disease[which(prevalences$sex_bias=="men")])
  dennone<-intersect(commondis,denprevalences$disease[which(denprevalences$sex_bias=="none")])
  catnone<-intersect(commondis,prevalences$disease[which(prevalences$sex_bias=="none")])
  biasprevalences<-list("Denmark - Women preference"=denwom,"Denmark - Men preference"=denmen, "Catalonia - Women preference"=catwom,"Catalonia - Men preference"=catmen)
  biasprevalences2<-list("Denmark - Women preference"=denwom,"Denmark - Men preference"=denmen, "Catalonia - Women preference"=catwom,"Catalonia - Men preference"=catmen, "Denmark - No preference"=dennone,"Catalonia - No preference"=catnone)
  ## Upset plot ##
  ## Gender associated ##
  upset_data<-fromList(biasprevalences)
  pdf("ManuscriptFiles/Plots/Disease_prevalence_biases_in_Catalonia_and_Denmark.pdf",width = 10,height = 6)
    upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = "")
  dev.off()
  ## Discordant disease prevalences
  setkey(denprevalences,"disease") ; setkey(prevalences,"disease")
  denwom_catmen<-cbind(as.character(disnames[intersect(denwom,catmen)]),denprevalences[intersect(denwom,catmen)][,c(2:4,8,9)],prevalences[intersect(denwom,catmen)][,c(3:4,8,9)])
  denmen_catwom<-cbind(as.character(disnames[intersect(catwom,denmen)]),denprevalences[intersect(catwom,denmen)][,c(2:4,8,9)],prevalences[intersect(catwom,denmen)][,c(3:4,8,9)])
  colnames(denwom_catmen)<-c("disease","diseasecategory","den_wom","den_men","den_pval","den_OR","cat_wom","cat_men","cat_pval","cat_OR")
  colnames(denmen_catwom)<-c("disease","diseasecategory","den_wom","den_men","den_pval","den_OR","cat_wom","cat_men","cat_pval","cat_OR")
  tabla<-rbind(denwom_catmen[,c(1,3,4,6:8,10)],denmen_catwom[,c(1,3,4,6:8,10)])
  write.table(tabla,"ManuscriptFiles/Results/Diseases_with_discordant_sex_related_risk_CatDen.txt",quote=F,sep="\t",row.names=F)
}
#### Plot gender-associated disease age of diagnoses differences in Catalonia and Denmark and look for differences between them ####
if(args[1]=="compare_age_of_diagnoses"){
  #### Catalonia ####
  ## Load diagnoses table ##
  dt<-fread("Diagnoses_20080101_20181231_first_diagnoses_of_each_disease.txt",stringsAsFactors = F,sep="|")
  ## Load the function ##
  ttest_cod<-function(dat) {
    ## Check that both genders are present ##
    if(length(unique(dat$sexe)) != 2) {return(list(diff_mean = NA_real_, p_value = NA_real_))}
    dat_W<-dat[sexe == "D", edad]; dat_M<-dat[sexe == "H", edad]
    ## At least two individuals by gender
    if(length(dat_W) < 2 | length(dat_M) < 2) {return(list(diff_mean = NA_real_, p_value = NA_real_))}
    ## Test t Welch
    t_test<-t.test(dat_W, dat_M, var.equal = FALSE)
    list(diff_mean = t_test$estimate[1] - t_test$estimate[2],p_value   = t_test$p.value)
  }
  results<-dt[, ttest_cod(.SD), by = cod]
  ## Add relevant columns ##
  summary_cod<-dt[, .(nMen    = sum(sexe == "H"),nWomen    = sum(sexe == "D"),
    MenMean = mean(edad[sexe == "H"], na.rm = TRUE),WomenMean = mean(edad[sexe == "D"], na.rm = TRUE),
    MenSD   = sd(edad[sexe == "H"], na.rm = TRUE),WomenSD   = sd(edad[sexe == "D"], na.rm = TRUE)
  ), by = cod][nWomen >= 2 & nMen >= 2]
  ## Put both tables together and correct for multiple testing ##
  catalonia_age<-merge(summary_cod, results, by = "cod", all.x = TRUE)
  catalonia_age[, fdr := p.adjust(p_value, method = "fdr")]
  colnames(catalonia_age)[c(1,8:10)]<-c("Code","DiffMean","pval","FDR")
  catalonia_age<-cbind(catalonia_age[,1:8],catalonia_age[,10])
  ## Save the table ##
  write.table(catalonia_age,"ManuscriptFiles/Results/SupplementaryData_gender_related_ageofdiagnosis_differences_Catalonia.txt",quote=F,sep="\t",row.names=F)
  ## Which diseases are diagnosed later in women? ##
  latterwomen<-catalonia_age[intersect(which(catalonia_age$FDR<=0.05),which(catalonia_age$DiffMean>0))] ; latterwomen<-latterwomen[order(abs(latterwomen$DiffMean),decreasing=T)]
  lattermen<-catalonia_age[intersect(which(catalonia_age$FDR<=0.05),which(catalonia_age$DiffMean<0))] ; lattermen<-lattermen[order(abs(lattermen$DiffMean),decreasing=T)]
  #### Boxplot age ####
  ## Add disease category ##
  dt<-merge(catalonia_age,prevalences[, .(disease, diseasecategory)],by.x = "Code",by.y = "disease",all.x = TRUE)
  ## Add disease category color ##
  dt[, color := distocol[Code]]
  dt[, diseasecategory := factor(diseasecategory)]
  if(length(which(is.na(dt$color)))>0){dt<-dt[-which(is.na(dt$color))]}
  ## Remove pregnancy related disorders for plotting ##
  if(length(grep("O",dt$Code))>0){dt<-dt[-grep("O",dt$Code)]}
  ## For each category, the proportion of diseases with later diagnosis en women is larger than the proportion of later diagnosis in men? ##
  dt2<-copy(dt)
  ## Mean difference (years)
  dt2[, yi := WomenMean - MenMean]
  ## Standard error of the mean difference (Welch formulation)
  dt2[, sei := sqrt((WomenSD^2 / nWomen) + (MenSD^2 / nMen))]
  dt2[, vi := sei^2]
  ## Basic quality filters
  dt2<-dt2[is.finite(yi) & is.finite(sei) & sei > 0 & nMen >= 5 &  nWomen >= 5 & !is.na(diseasecategory)]
  ## Weighted least squares (inverse-variance weighted) ##
  dt2[, diseasecategory := as.factor(diseasecategory)]
  ## Estimate mean difference (in years) per category
  fit_wls<-lm(yi ~ 0 + diseasecategory, data = dt2, weights = 1/vi)
  summ<-summary(fit_wls)
  res_category<-as.data.table(coef(summ), keep.rownames = "term")
  setnames(res_category, c("term","Estimate","Std. Error","t value","Pr(>|t|)"), c("term","beta_years","se","t","p_value"))
  ## Extract clean category name
  res_category[, diseasecategory := sub("^diseasecategory", "", term)]
  ## 95% CI
  res_category[, `:=`(ci_low  = beta_years - 1.96 * se, ci_high = beta_years + 1.96 * se)]
  ## FDR correction across categories
  res_category[, FDR_cat := p.adjust(p_value, method = "BH")]
  ## Add number of diseases per category
  res_category<-merge(res_category, dt2[, .(n_diseases = .N), by = diseasecategory], by = "diseasecategory", all.x = TRUE)[order(FDR_cat)]
  ## Final results table
  res_category
  ## Directional bias among significant diseases (binomial test)
  dt_sig<-dt[FDR <= 0.05]
  dt_sig[, direction := fifelse(
    WomenMean - MenMean > 0, "Later_in_Women",
    fifelse(WomenMean - MenMean < 0, "Later_in_Men", NA_character_)
  )]
  dt_sig<-dt_sig[!is.na(direction)]
  res_binom<-dt_sig[
    ,
    {
      n_total<-.N
      n_women<-sum(direction == "Later_in_Women")
      bt_women<-binom.test(n_women, n_total, p = 0.5, alternative = "greater")
      bt_men  <-binom.test(n_women, n_total, p = 0.5, alternative = "less")
      .(n_sig = n_total, n_later_women = n_women, n_later_men = n_total - n_women, prop_later_women = n_women / n_total, p_women = bt_women$p.value, p_men = bt_men$p.value)
    },
    by = diseasecategory
  ]
  res_binom[, `:=`(FDR_women = p.adjust(p_women, method = "BH"),FDR_men   = p.adjust(p_men, method = "BH"))]
  ## Integrated category-level results
  res_category<-as.data.table(res_category)
  res_binom   <-as.data.table(res_binom)
  ## Clean category names (remove possible leading spaces)
  res_category[, diseasecategory := trimws(diseasecategory)]
  res_binom[, diseasecategory := trimws(as.character(diseasecategory))]
  ## Merge both tables
  supp_table3<-merge(
    res_category[, .(diseasecategory, n_diseases, beta_years, se, ci_low, ci_high, p_value_category = p_value, FDR_category = FDR_cat)],
    res_binom[, .(diseasecategory, n_sig, n_later_women, n_later_men, prop_later_women, p_women, p_men, FDR_women, FDR_men)],
    by = "diseasecategory", all.x = TRUE
  )
  ## Order by FDR (category-level model)
  setorder(supp_table3, FDR_category)
  supp_table3[, `:=`(beta_years = round(beta_years, 2), se = round(se, 2), ci_low = round(ci_low, 2), ci_high = round(ci_high, 2), prop_later_women = round(prop_later_women, 3))]
  ## Define direction label for easier interpretation
  supp_table3[, dominant_direction := fifelse(
    FDR_category <= 0.05 & beta_years > 0, "Later in women",
    fifelse(FDR_category <= 0.05 & beta_years < 0, "Later in men", "Not significant")
  )]
  ## Reorder columns for clarity
  setcolorder(supp_table3, c("diseasecategory", "n_diseases", "beta_years", "ci_low", "ci_high", "p_value_category", "FDR_category", "n_sig", "n_later_women", "n_later_men", "prop_later_women", "FDR_women", "FDR_men", "dominant_direction"))
  ## Write table
  fwrite(supp_table3, "ManuscriptFiles/Results/Supplementary_Table_AgeDiagnosisDifference_categorylevel.csv")
  ## Plot with all the diseases ##
  pdf(file="ManuscriptFiles/Plots/Gender_differences_agediagnosis_all.pdf",width = 12,height = 8)
  ggplot(dt, aes(x = DiffMean, y = diseasecategory, fill = diseasecategory)) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.1, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = diseasecategory), size = 2, alpha = 0.5, width = 0.2) +
    scale_fill_manual(values = unique(dt$color[order(dt$diseasecategory)])) +
    scale_color_manual(values = unique(dt$color[order(dt$diseasecategory)])) +
    labs(x = "DiffMean (Women − Men)", y = "Disease category", fill = "Category", color = "Category") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90"))
  dev.off()
  ## Only FDR<=0.05 ##
  dt_sig<-dt[FDR <= 0.05]
  pdf(file="ManuscriptFiles/Plots/Gender_differences_agediagnosis_significant.pdf",width = 12,height = 8)
  ggplot(dt_sig, aes(x = DiffMean, y = diseasecategory, fill = diseasecategory)) +
    geom_violin(alpha = 0.4, color = NA) +
    geom_boxplot(width = 0.1, alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = diseasecategory), size = 2, alpha = 0.5, width = 0.2) +
    scale_fill_manual(values = unique(dt_sig$color[order(dt_sig$diseasecategory)])) +
    scale_color_manual(values = unique(dt_sig$color[order(dt_sig$diseasecategory)])) +
    labs(x = "DiffMean (Women − Men)", y = "Disease category", fill = "Category", color = "Category") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90"))
  dev.off()
  
  #### Age vs. Prevalence (sum of the number of patients) ####
  dt<-merge(catalonia_age,prevalences[, .(disease, diseasecategory)],by.x = "Code",by.y = "disease",all.x = TRUE)
  ## Total prevalence ## 
  dt[, prevalence := nMen + nWomen]
  ## Add the color ##
  dt[, color := distocol[Code]]
  ## Identify significances ##
  dt[, signif := ifelse(FDR <= 0.05, "FDR ≤ 0.05", "Not significant")]
  dt[, signif := factor(signif, levels = c("Not significant", "FDR ≤ 0.05"))]
  ## Remove pregnancy related disorders for plotting ##
  if(length(grep("O",dt$Code))>0){dt<-dt[-grep("O",dt$Code)]}
  write.table(dt,"ManuscriptFiles/IntermediateFiles/AgeDifferencesCompleteListCatalonia.txt",quote=F,sep="\t",row.names=F)
  
  ## Plot ##
  pdf(file="ManuscriptFiles/Plots/Gender_differences_prevalence_number_vs_age.pdf",width = 12,height = 8)
  ggplot(dt, aes(x = prevalence, y = DiffMean)) +
    geom_point(aes(color = diseasecategory,shape = signif),size = 3.5,alpha = 0.65) +
    scale_x_log10() +
    scale_color_manual(values = unique(dt$color[order(dt$diseasecategory)])) +
    scale_shape_manual(values = c("Not significant" = 16, "FDR ≤ 0.05" = 17)) +
    labs(x = "Prevalence (log scale)",y = "DiffMean (Women − Men)",color = "Disease category",shape = "Significance") +
    theme_bw()
  dev.off()
  
  #### Age vs. Prevalence differences (OR) ####
  ## Merge by disease ##
  dt<-merge(catalonia_age,prevalences,by.x = "Code",by.y = "disease",all = FALSE)
  ## Log(OR) ##
  dt[, logOR := log(OR)]
  ## Add disease category name ##
  dt[, catename := diseasecategory]
  ## Add the color of each disease ##
  dt[, color := distocol[Code]]
  write.table(dt,"ManuscriptFiles/IntermediateFiles/AgeAndPrevalenceDifferencesCompleteListCatalonia.txt",quote=F,sep="\t",row.names=F)
  ## Remove NAs ##
  xlims<-range(dt$DiffMean, na.rm = TRUE)
  ylims<-range(dt$logOR, na.rm = TRUE)
  ## Create a color vector ##
  cat_colors<-unique(dt[, .(catename = diseasecategory, color)])
  cat_colors_vec<-setNames(cat_colors$color, cat_colors$catename)
  ## Plot ##
  p<-ggplot(dt, aes(x = DiffMean, y = logOR)) +
    geom_point(aes(color = catename),size = 3.5,alpha = 0.65) +
    scale_color_manual(values = cat_colors_vec) +
    labs(x = "DiffMean (Women − Men)",y = "log(OR) (Women vs. Men)",color = "Disease category") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(color = "grey90"))
  
  pdf(file="ManuscriptFiles/Plots/Gender_differences_prevalence_vs_age.pdf",width = 12,height = 8)
    print(p)
  dev.off()
  
  #### Age vs. Prevalence differences (OR) - Interactive ####
  dig3dis<-fread("Data/ICD10_three_digits_names.txt",stringsAsFactors = F)
  tdig3dis<-dig3dis$Name ; names(tdig3dis)<-dig3dis$Code
  dt[, disease_name := tdig3dis[Code]]
  p<-ggplot(dt, aes(x = DiffMean, y = logOR)) +
    geom_point(aes(color = catename, text = paste0("Code: ", Code, "<br>", "Disease: ", disease_name, "<br>", "Category: ", catename, "<br>", "DiffMean: ", round(DiffMean, 2), "<br>", "log(OR): ", round(logOR, 2))), size = 3.5, alpha = 0.65) +
    scale_color_manual(values = cat_colors_vec) +
    labs(x = "DiffMean (Women − Men)", y = "log(OR) (Women vs. Men)", color = "Disease category") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey90"))
  ggplotly(p, tooltip = "text")
  
  #### Compare with Denmark ####
  denmark_age<-fread("Epidemiology/41467_2019_8475_MOESM5_ESM.txt",stringsAsFactors = F) ; setkey(denmark_age,"Code")
  ## Make the difference be women - men ##
  denmark_age$DiffMean<-denmark_age$DiffMean*(-1)
  setkey(catalonia_age,"Code")
  commondis<-intersect(denmark_age$Code,catalonia_age$Code)
  cdenmark_age<-denmark_age[commondis] ; ccatalonia_age<-catalonia_age[commondis]
  ## Ensure both tables present the same ordering ##
  stopifnot(all(cdenmark_age$Code == ccatalonia_age$Code))
  ## Combine both tables ##
  dt<-data.table(Code = cdenmark_age$Code,Diff_Denmark   = cdenmark_age$DiffMean,Diff_Catalonia = ccatalonia_age$DiffMean,FDR_Denmark    = cdenmark_age$FDR,FDR_Catalonia  = ccatalonia_age$FDR)
  ## Get significant categories ##
  dt[, signif := fifelse(
    FDR_Denmark <= 0.05 & FDR_Catalonia <= 0.05, "Both",
    fifelse(FDR_Denmark <= 0.05, "Denmark",
            fifelse(FDR_Catalonia <= 0.05, "Catalonia", "None"))
  )]
  ## Convert into factor for ordering labels ##
  dt[, signif := factor(signif,levels = c("None", "Denmark", "Catalonia", "Both"))]
  ## Add color information ##
  dt[, color := distocol[Code]]
  dt[, catename := as.character(catename[dt$Code])]
  write.table(dt,"ManuscriptFiles/IntermediateFiles/AgeDifferencesCompleteListCataloniavsDenmark.txt",quote=F,sep="\t",row.names=F)
  ## Plot ##
  pdf(file="ManuscriptFiles/Plots/Age_difference_correlation_countries.pdf",width = 10,height = 8)
    ggplot(dt, aes(x = Diff_Denmark, y = Diff_Catalonia)) +
      geom_point(aes(shape = signif, color = color),alpha = 0.6,size = 2.5) +
      scale_color_identity() +
      scale_shape_manual(values = c("None" = 16,"Denmark" = 17,"Catalonia" = 15,"Both" = 18)) +
      labs(x = "Age difference (Denmark)",y = "Age difference (Catalonia)",shape = "Significance (FDR ≤ 0.05)") +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(color = "grey90"))
  dev.off()
  ## Upset plot ##
  lista<-list("Later Women (Catalonia)"=ccatalonia_age$Code[intersect(which(ccatalonia_age$FDR<=0.05),which(ccatalonia_age$DiffMean>0))],
              "Later Men (Catalonia)"=ccatalonia_age$Code[intersect(which(ccatalonia_age$FDR<=0.05),which(ccatalonia_age$DiffMean<0))],
              "Later Women (Denmark)"=cdenmark_age$Code[intersect(which(cdenmark_age$FDR<=0.05),which(cdenmark_age$DiffMean>0))],
              "Later Men (Denmark)"=cdenmark_age$Code[intersect(which(cdenmark_age$FDR<=0.05),which(cdenmark_age$DiffMean<0))])
  upset_data<-fromList(lista)
  pdf("ManuscriptFiles/Plots/Disease_age_gender_biases_in_Catalonia_and_Denmark.pdf",width = 10,height = 6)
    upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = "")
  dev.off()
  ## Discordant disease age differences
  setkey(denprevalences,"disease") ; setkey(prevalences,"disease")
  tabla<-dt[c(intersect(intersect(which(dt$Diff_Catalonia>0),which(dt$FDR_Catalonia<=0.05)),intersect(which(dt$Diff_Denmark<0),which(dt$FDR_Denmark<=0.05))),
    intersect(intersect(which(dt$Diff_Catalonia<0),which(dt$FDR_Catalonia<=0.05)),intersect(which(dt$Diff_Denmark>0),which(dt$FDR_Denmark<=0.05))))]
  tabla<-tabla[,c(1:5)]
  tabla$Diff_Denmark<-round(tabla$Diff_Denmark,2)
  tabla$Diff_Catalonia<-round(tabla$Diff_Catalonia,2)
  tabla$Code<-as.character(disnames[tabla$Code])
  colnames(tabla)[1]<-"Disease"
  write.table(tabla,"ManuscriptFiles/Results/Diseases_with_discordant_sex_related_ageDiagnosis_CatDen.txt",quote=F,sep="\t",row.names=F)
}

#### Calculate correlations between disease prevalence and number of comorbidities (separately for source and sink nodes) ####
if(args[1]=="calculate_correlations_between_diseaseprevalence_and_numberofcomorbidities"){
  #### Correlations between number of comorbidities and disease prevalence ####
  prevalences<-fread("Data/ICD10_prevalence_Catalonia.txt",stringsAsFactors = F,sep="\t")
  colnames(prevalences)[1]<-"diseases" ; setkey(prevalences,"diseases")
  ficheros<-list.files("Results/shrinkage_rr_events/")
  comparisons<-c("both","women","men")
  ## Incremental_time_windows ##
  ficheros1<-ficheros[grep("RR_contingency_",ficheros)]
  for(z in comparisons){
    fichs<-ficheros1[grep(paste("_",z,sep=""),ficheros1)]
    lesult<-list()
    for(a in fichs){
      ## Load the three tables ##
      bayesall<-fread(paste("Results/shrinkage_rr_events/",a,sep=""),stringsAsFactors = F,sep="\t")
      ## Get the significant associations ##
      sbayesall<-bayesall[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
      ## correlations ##
      if(gsub("_.+","",gsub("RR_contingency_","",a))=="both"){patsbydisease<-prevalences[,1:2] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub("_.+","",gsub("RR_contingency_","",a))=="women"){patsbydisease<-prevalences[,c(1,3)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub("_.+","",gsub("RR_contingency_","",a))=="men"){patsbydisease<-prevalences[,c(1,4)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      ## significant bayes all ##
      ## @@ @@ @@ @ @ @@ @@ @@ ##
      ## Counts in Disease1 ##
      count_d1<-sbayesall %>% dplyr::count(disease_a, name = "n_Disease1") %>% rename(diseases = disease_a)
      ## Counts in Disease2 ##
      count_d2<-sbayesall %>% dplyr::count(disease_b, name = "n_Disease2") %>% rename(diseases = disease_b)
      ## Add to patsbydisease
      resulta<-patsbydisease %>% left_join(count_d1, by = "diseases") %>% left_join(count_d2, by = "diseases") %>%
        mutate(n_Disease1 = coalesce(n_Disease1, 0L),n_Disease2 = coalesce(n_Disease2, 0L))
      if(length(intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)))>0){resulta<-resulta[-intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)),]}
      lesult[[gsub("_shrunk.+","",gsub(paste0("RR_contingency_",z,"_"),"",a))]]<-resulta
    }
    ## Plot ##
    pdf(file=paste("ManuscriptFiles/Disease_prevalence_number_comorbidities_plot/",z,"_incremental_timewindows.pdf",sep=""),width = 18,height = 8)
      par(mfrow = c(2, 5))
      ## 1 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 2 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 3 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 4 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 5 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 2 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 3 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 4 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 10 ##
      plot(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 5 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease2)$p.value,6),sep=""),cex=0.7)
    dev.off()
    print(paste(z,"finished!"))
  }
  ## Continuous_time_windows ##
  ficheros2<-c(ficheros[grep("_0_1_",ficheros)],ficheros[grep("_1_2_",ficheros)],ficheros[grep("_2_3_",ficheros)],ficheros[grep("_3_4_",ficheros)],ficheros[grep("_4_5_",ficheros)])
  for(z in comparisons){
    fichs<-ficheros2[grep(paste("_",z,sep=""),ficheros2)]
    lesult<-list()
    for(a in fichs){
      ## Load the three tables ##
      bayesall<-fread(paste("Results/shrinkage_rr_events/",a,sep=""),stringsAsFactors = F,sep="\t")
      ## Get the significant associations ##
      sbayesall<-bayesall[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
      ## Correlations ##
      if(gsub("_.+","",gsub("RR_contingency_","",a))=="both"){patsbydisease<-prevalences[,1:2] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub("_.+","",gsub("RR_contingency_","",a))=="women"){patsbydisease<-prevalences[,c(1,3)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub("_.+","",gsub("RR_contingency_","",a))=="men"){patsbydisease<-prevalences[,c(1,4)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      ## Counts in Disease1 ##
      count_d1<-sbayesall %>% dplyr::count(disease_a, name = "n_Disease1") %>% rename(diseases = disease_a)
      ## Counts in Disease2 ##
      count_d2<-sbayesall %>% dplyr::count(disease_b, name = "n_Disease2") %>% rename(diseases = disease_b)
      ## Add to patsbydisease
      resulta<-patsbydisease %>% left_join(count_d1, by = "diseases") %>% left_join(count_d2, by = "diseases") %>%
        mutate(n_Disease1 = coalesce(n_Disease1, 0L),n_Disease2 = coalesce(n_Disease2, 0L))
      if(length(intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)))>0){resulta<-resulta[-intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)),]}
      lesult[[gsub("_shrunk.+","",gsub(paste0("RR_contingency_",z,"_"),"",a))]]<-resulta
    }
    ## Plot ##
    pdf(file=paste("ManuscriptFiles/Disease_prevalence_number_comorbidities_plot/",z,"_consecutive_timewindows.pdf",sep=""),width = 18,height = 8)
      par(mfrow = c(2, 5))  # 2 filas, 3 columnas
      ## 1 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 1 - 2 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 2 - 3 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 3 - 4 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 4 - 5 years\nsource disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease1)$estimate,3), "\npval = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 1 - 2 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 2 - 3 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 3 - 4 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 10 ##
      plot(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 4 - 5 years\nsink disease", xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease2)$estimate,3), "\npval = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease2)$p.value,6),sep=""),cex=0.7)
    dev.off()
    print(paste(z,"finished!"))
  }
}

#### Get the biological clock ####
if(args[1] == "biological_clock") {
  ## ------------------------------------------------------------------ ##
  ## ICD-10 category mapping (shared with rest of manuscript_analyses)   ##
  ## ------------------------------------------------------------------ ##
  code<-c(paste0("A0",0:9), paste0("A",10:99),
            paste0("B0",0:9), paste0("B",10:99),
            paste0("C0",0:9), paste0("C",10:99),
            paste0("D0",0:9), paste0("D",10:48), paste0("D",50:89),
            paste0("E0",0:9), paste0("E",10:99),
            paste0("F0",0:9), paste0("F",10:99),
            paste0("G0",0:9), paste0("G",10:99),
            paste0("H0",0:9), paste0("H",10:59), paste0("H",60:95),
            paste0("I0",0:9), paste0("I",10:99),
            paste0("J0",0:9), paste0("J",10:99),
            paste0("K0",0:9), paste0("K",10:93),
            paste0("L0",0:9), paste0("L",10:99),
            paste0("M0",0:9), paste0("M",10:99),
            paste0("N0",0:9), paste0("N",10:99),
            paste0("O0",0:9), paste0("O",10:99),
            paste0("P0",0:9), paste0("P",10:96),
            paste0("Q0",0:9), paste0("Q",10:99),
            paste0("R0",0:9), paste0("R",10:99),
            paste0("S0",0:9), paste0("S",10:99),
            paste0("T0",0:9), paste0("T",10:98))
  cate<-c(rep("I",200), rep("II",149), rep("III",40),
            rep("IV",100), rep("V",100), rep("VI",100),
            rep("VII",60), rep("VIII",36), rep("IX",100),
            rep("X",100), rep("XI",94), rep("XII",100),
            rep("XIII",100), rep("XIV",100), rep("XV",100),
            rep("XVI",97), rep("XVII",100), rep("XVIII",100),
            rep("XIX",199))
  catname<-c("Infectious","Neoplasms","Blood/Immune",
               "Endocrine/Metabolic","Mental","Nervous",
               "Eye","Ear","Circulatory","Respiratory",
               "Digestive","Skin","Musculoskeletal",
               "Genitourinary","Pregnancy","Perinatal",
               "Congenital","Symptoms","Injury")
  names(catname)<-unique(cate)
  names(cate)   <-code
  SIG_FILTER<-quote(lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)
  
  #### Load cumulative windows and identify persistent pairs ####
  windows_cumul<-c("0_1","0_2","0_3","0_4","0_5")
  all_cumul<-rbindlist(lapply(windows_cumul, function(w){
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) { cat("  Missing:", f, "\n"); return(NULL) }
    dt<-fread(f)
    dt<-dt[eval(SIG_FILTER)]
    dt[, `:=`(window = w, window_num = as.integer(gsub("0_","",w)), pair = paste(disease_a, disease_b, sep = "_"))]
    dt
  }))
  pairs_per_window<-all_cumul[, .(n_windows = uniqueN(window)), by = pair]
  persistent_pairs<-pairs_per_window[n_windows == 5, pair]
  cat(sprintf("Pairs with data in >=1 window: %d\n", nrow(pairs_per_window)))
  cat(sprintf("Persistent pairs (all 5 windows): %d\n", length(persistent_pairs)))
  
  #### Compute weighted slope for each persistent pair ####
  rr_traj<-all_cumul[pair %in% persistent_pairs, .(pair, window_num, RR_shrunk, SE_RR_shrunk, disease_a, disease_b)]
  rr_traj[, catA := catname[cate[disease_a]]]
  rr_traj[, catB := catname[cate[disease_b]]]
  slopes<-rr_traj[, {
    x <-window_num
    y <-RR_shrunk
    se<-SE_RR_shrunk
    if(any(is.na(se)) || any(se <= 0) || length(x) < 3){
      .(slope = NA_real_, slope_p = NA_real_)
    }else{
      w_fit<-1/se^2
      fit <-lm(y ~ x, weights = w_fit)
      sm <-summary(fit)$coefficients
      .(slope   = sm["x","Estimate"], slope_p = sm["x","Pr(>|t|)"])
    }
  }, by = .(pair, catA, catB)]
  slopes[!is.na(slope_p), slope_fdr := p.adjust(slope_p, method = "BH")]
  slopes[is.na(slope_p),  slope_fdr := NA_real_]
  
  #### Category-level slope summary + statistical tests ####
  traj_by_cat<-slopes[!is.na(catA) & !is.na(slope), .(
    n_total = .N,
    med_slope = median(slope),
    pct_decr = mean(slope < -0.05 & !is.na(slope_fdr) & slope_fdr < 0.05) * 100,
    pct_stable = mean(abs(slope) <= 0.05 | is.na(slope_fdr) | slope_fdr >= 0.05) * 100,
    pct_incr = mean(slope >  0.05 & !is.na(slope_fdr) & slope_fdr < 0.05) * 100
  ), by = catA][order(med_slope)]
  ## Kruskal-Wallis global
  kw_global<-kruskal.test(slope ~ catA, data = slopes[!is.na(catA) & !is.na(slope)])
  cat(sprintf("Kruskal-Wallis across categories: chi2=%.1f, p=%.2e\n", kw_global$statistic, kw_global$p.value))
  ## Post-hoc: each category vs Endocrine/Metabolic (most stable)
  ref_slopes<-slopes[catA == "Endocrine/Metabolic", slope]
  posthoc<-slopes[!is.na(catA) & !is.na(slope) & catA != "Endocrine/Metabolic", {
    wt<-wilcox.test(slope, ref_slopes, exact = FALSE)
    .(p_vs_ref = wt$p.value, med_slope = median(slope))
  }, by = catA]
  posthoc[, p_adj := p.adjust(p_vs_ref, method = "BH")]
  ## Assign biological clock groups
  traj_by_cat[, group := fcase(med_slope < -0.045, "Fast attenuation\n(episodic/reactive)", med_slope < -0.032, "Moderate attenuation\n(mixed)", rep(TRUE, .N), "Slow attenuation\n(chronic/progressive)")]
  ## Save
  fwrite(traj_by_cat, "ManuscriptFiles/Results/Biological_clock_slopes_by_category.txt", sep = "\t", quote = FALSE)
  
  #### Conditional windows - persistent and late-emerging risk ####
  windows_cont<-c("1_2","2_3","3_4","4_5")
  all_cont<-rbindlist(lapply(windows_cont, function(w) {
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) return(NULL)
    dt<-fread(f)
    dt<-dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
    dt[, `:=`(window = w, window_num = as.integer(gsub("_.+","",w)), pair = paste(disease_a, disease_b, sep = "_"), sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
    dt[, catA := catname[cate[disease_a]]]
    dt[, catB := catname[cate[disease_b]]]
    dt
  }))
  ## Keep only pairs present in all 4 conditional windows
  pairs_cont4<-all_cont[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]
  cont_wide<-data.table::dcast(all_cont[pair %in% pairs_cont4, .(pair, window, sig, RR_shrunk, catA, catB)], pair + catA + catB ~ window, value.var = c("sig", "RR_shrunk"))
  cont_wide[, pattern := fcase(
    (sig_1_2 == TRUE | sig_2_3 == TRUE) & (sig_3_4 == FALSE & sig_4_5 == FALSE), "Early only (risk fades)",
    (sig_1_2 == TRUE | sig_2_3 == TRUE) & (sig_3_4 == TRUE  | sig_4_5 == TRUE), "Persistent conditional risk",
    (sig_1_2 == FALSE & sig_2_3 == FALSE) & (sig_3_4 == TRUE  | sig_4_5 == TRUE), "Late emerging risk",
    rep(TRUE, .N), "Not significant"
  )]
  print(cont_wide[, .N, by = pattern][order(-N)])
  ## Get the percentage of comorbidities for each index category (Skin, Eye, Mental...) in each group ##
  ## The sum of the pct in each category should be 100 ##
  pat_by_cat<-cont_wide[!is.na(catA) & !is.na(pattern), .(
    n_total = .N,
    pct_early = mean(pattern == "Early only (risk fades)")*100,
    pct_persist = mean(pattern == "Persistent conditional risk")*100,
    pct_late = mean(pattern == "Late emerging risk")*100,
    pct_notsig = mean(pattern == "Not significant")*100
  ), by = catA][order(-pct_persist)]
  fwrite(cont_wide[, .(pair, catA, catB, pattern, RR_shrunk_1_2, RR_shrunk_2_3, RR_shrunk_3_4, RR_shrunk_4_5)], "ManuscriptFiles/Results/Biological_clock_conditional_patterns.txt", sep = "\t", quote = FALSE)
  ## Top late-emerging pairs in biphasic categories
  top_late<-cont_wide[pattern == "Late emerging risk" & catA %in% c("Circulatory","Endocrine/Metabolic","Skin","Infectious","Injury") & !is.na(RR_shrunk_4_5), .(pair, catA, catB, RR_3_4 = RR_shrunk_3_4, RR_4_5 = RR_shrunk_4_5)][order(-RR_4_5)][1:30]
  fwrite(top_late, "ManuscriptFiles/Results/Biological_clock_late_emerging_top.txt", sep = "\t", quote = FALSE)
  
  #### Combine dimensions — 2D clock space ####
  clock_2d<-merge(traj_by_cat[, .(catA, med_slope, group, n_total)], pat_by_cat[, .(catA, pct_persist, pct_late)], by = "catA")
  ## Get the median of slope and pct_persistent to create the quadrants ##
  x_mid<-median(clock_2d$med_slope)
  y_mid<-median(clock_2d$pct_persist)
  ## Get the classification on each quadrant ##
  clock_2d[, quadrant_test := fcase(
    med_slope <= x_mid & pct_persist >= y_mid, "Biphasic",
    med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive",
    med_slope <= x_mid & pct_persist <  y_mid, "Purely episodic",
    rep(TRUE, .N), "Chronic stable"
  )]
  print(clock_2d[, .(catA, med_slope, pct_persist, quadrant_test)])
  ## Summary table for the paper
  quadrant_cats<-list(
    "Purely episodic" = clock_2d$catA[which(clock_2d$quadrant_test=="Purely episodic")],
    "Chronic stable" = clock_2d$catA[which(clock_2d$quadrant_test=="Chronic stable")],
    "Chronic progressive" = clock_2d$catA[which(clock_2d$quadrant_test=="Chronic progressive")],
    "Biphasic" = clock_2d$catA[which(clock_2d$quadrant_test=="Biphasic")]
  )
  
  #### Find examples that are representative and interesting for each quadrant ####
  ## Load data needed
  rr_01<-fread("Results/shrinkage_rr_events/RR_contingency_both_0_1_shrunk.txt")
  rr_05<-fread("Results/shrinkage_rr_events/RR_contingency_both_0_5_shrunk.txt")
  rr_01[, pair := paste(disease_a, disease_b, sep="_")]
  rr_05[, pair := paste(disease_a, disease_b, sep="_")]
  ## Add quadrant info to slopes and cont_wide
  slopes[, quadrant := clock_2d$quadrant_test[match(catA, clock_2d$catA)]]
  cont_wide[, quadrant := clock_2d$quadrant_test[match(catA, clock_2d$catA)]]
  ## Purely episodic ##
  ## High RR in 0-1, strong decay, low conditional persistence
  episodic_cands<-slopes[catA %in% quadrant_cats[["Purely episodic"]] & !is.na(slope) & slope < 0, .(pair, catA, slope)][rr_01[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100,.(pair, RR_01 = RR_shrunk)], on = "pair", nomatch = 0][cont_wide[pattern == "Not significant", .(pair)], on = "pair", nomatch = 0][order(-RR_01)][1:10, .(pair, catA, RR_01, slope)]
  print(episodic_cands)
  
  ## Chronic stable ##
  ## High RR in 0-5, flat slope, some conditional persistence
  stable_cands<-slopes[catA %in% quadrant_cats[["Chronic stable"]] & !is.na(slope), .(pair, catA, slope)][rr_05[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100, .(pair, RR_05 = RR_shrunk)], on = "pair", nomatch = 0][order(-RR_05)][1:10, .(pair, catA, RR_05, slope)]
  print(stable_cands)
  
  ## Chronic progressive ##
  ## Significant in all 4 conditional windows, stable or increasing cumulative RR
  progressive_cands<-cont_wide[catA %in% quadrant_cats[["Chronic progressive"]] & sig_1_2 == TRUE & sig_2_3 == TRUE & sig_3_4 == TRUE & sig_4_5 == TRUE, .(pair, catA, RR_1_2 = RR_shrunk_1_2, RR_4_5 = RR_shrunk_4_5)][slopes[, .(pair, slope)], on = "pair", nomatch = 0][order(slope)][1:10, .(pair, catA, RR_1_2, RR_4_5, slope)]
  print(progressive_cands)
  
  ## Biphasic ##
  ## Non-sig or low RR in 1-2 conditional, high RR in 4-5 conditional
  biphasic_cands<-cont_wide[catA %in% quadrant_cats[["Biphasic"]] & sig_1_2 == FALSE & sig_2_3 == FALSE & sig_4_5 == TRUE & !is.na(RR_shrunk_4_5), .(pair, catA, RR_1_2 = RR_shrunk_1_2, RR_4_5 = RR_shrunk_4_5)][order(-RR_4_5)][1:10, .(pair, catA, RR_1_2, RR_4_5)]
  print(biphasic_cands)
  
  #### Find more extreme examples ####
  ## Purely episodic ##
  episodic_best<-slopes[catA %in% quadrant_cats[["Purely episodic"]] & !is.na(slope) & slope < 0 & slope_fdr < 0.05, .(pair, catA, slope)][rr_01[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100, .(pair, RR_01 = RR_shrunk)], on = "pair", nomatch = 0][cont_wide[pattern == "Not significant", .(pair)], on = "pair", nomatch = 0][order(-RR_01)][1:15, .(pair, catA, RR_01, slope)]
  ## Chronic stable ##
  stable_best<-cont_wide[catA %in% quadrant_cats[["Chronic stable"]] & pattern == "Persistent conditional risk" & sig_1_2 == TRUE & sig_2_3 == TRUE & sig_3_4 == TRUE & sig_4_5 == TRUE, .(pair, catA, RR_1_2 = RR_shrunk_1_2, RR_4_5 = RR_shrunk_4_5, ratio = RR_shrunk_4_5 / RR_shrunk_1_2)][slopes[slope_fdr >= 0.05, .(pair, slope)], on = "pair", nomatch = 0][order(-ratio)][1:15, .(pair, catA, RR_1_2, RR_4_5, ratio, slope)]
  ## Chronic progressive ##
  progressive_best<-cont_wide[catA %in% quadrant_cats[["Chronic progressive"]] & sig_1_2 == TRUE & sig_2_3 == TRUE & sig_3_4 == TRUE & sig_4_5 == TRUE, .(pair, catA, RR_1_2 = RR_shrunk_1_2, RR_4_5 = RR_shrunk_4_5, ratio = RR_shrunk_4_5 / RR_shrunk_1_2)][slopes[slope_fdr < 0.05 & slope < 0, .(pair, slope)], on = "pair", nomatch = 0][order(-ratio)][1:15, .(pair, catA, RR_1_2, RR_4_5, ratio, slope)]
  ## Biphasic ##
  biphasic_best<-cont_wide[catA %in% quadrant_cats[["Biphasic"]] & sig_1_2 == FALSE & sig_2_3 == FALSE & sig_4_5 == TRUE & !is.na(RR_shrunk_1_2) & RR_shrunk_1_2 > 0, .(pair, catA, RR_1_2 = RR_shrunk_1_2, RR_4_5 = RR_shrunk_4_5, ratio = RR_shrunk_4_5 / RR_shrunk_1_2)][order(-ratio)][1:15, .(pair, catA, RR_1_2, RR_4_5, ratio)]
  
  cat("PURELY EPISODIC")
  print(episodic_best)
  cat("CHRONIC STABLE")
  print(stable_best)
  cat("CHRONIC PROGRESSIVE")
  print(progressive_best)
  cat("BIPHASIC")
  print(biphasic_best)
  
  examples<-c(
    "Purely episodic" = "A63\u2192A64 (Anogenital warts\u2192Other venereal diseases)",
    "Chronic stable" = "F10\u2192K70 (Alcohol disorder\u2192Liver cirrhosis)",
    "Chronic progressive" = "J01\u2192J32 (Acute sinusitis\u2192Chronic sinusitis)",
    "Biphasic" = "R06\u2192J96 (Dyspnoea\u2192Chronic Respiratory failure)"
  )
  summary_table<-rbindlist(lapply(names(quadrant_cats), function(q) {
    cats<-quadrant_cats[[q]]
    data.table(
      Pattern = q,
      Categories = paste(cats, collapse = ", "),
      Example = examples[q],
      N_pairs_early = cont_wide[catA %in% cats & pattern == "Early only (risk fades)", .N],
      N_pairs_persistent = cont_wide[catA %in% cats & pattern == "Persistent conditional risk", .N],
      N_pairs_late = cont_wide[catA %in% cats & pattern == "Late emerging risk", .N]
    )
  }))
  fwrite(summary_table, "ManuscriptFiles/Results/Biological_clock_summary_table.txt", sep = "\t", quote = FALSE)
  cat("Saved: Biological_clock_summary_table.txt\n")
  print(summary_table)
  
  #### Get the RR of each of the examples ####
  ## Accumulated ##
  print(all_cumul[pair == "A63_A64", .(window_num, RR = round(RR_shrunk,2), CI_low = round(CI_low_RR_shrunk,2), cases_event)])
  ## Conditional ##
  lapply(c("1_2","2_3","3_4","4_5"), function(w) {
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt<-fread(f)
    dt[disease_a == "A63" & disease_b == "A64",.(window=w, RR=round(RR_shrunk,2), cases_event, sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  }) |> rbindlist()
  
  ## Accumulated ##
  print(all_cumul[pair == "F10_K70", .(window_num, RR = round(RR_shrunk,2), CI_low = round(CI_low_RR_shrunk,2), cases_event)])
  ## Conditional ##
  lapply(c("1_2","2_3","3_4","4_5"), function(w) {
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt<-fread(f)
    dt[disease_a == "F10" & disease_b == "K70", .(window=w, RR=round(RR_shrunk,2), cases_event, sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  }) |> rbindlist()
  
  ## Accumulated ##
  print(all_cumul[pair == "J01_J32", .(window_num, RR = round(RR_shrunk,2), CI_low = round(CI_low_RR_shrunk,2), cases_event)])
  ## Conditional ##
  lapply(c("1_2","2_3","3_4","4_5"), function(w) {
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    dt<-fread(f)
    dt[disease_a == "J01" & disease_b == "J32", .(window=w, RR=round(RR_shrunk,2), cases_event, sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  }) |> rbindlist()
  
  #### Sex-stratified slopes ####
  all_cumul_sex<-rbindlist(lapply(c("both","women","men"), function(pop) {
    rbindlist(lapply(windows_cumul, function(w) {
      f<-sprintf("Results/shrinkage_rr_events/RR_contingency_%s_%s_shrunk.txt", pop, w)
      if(!file.exists(f)){return(NULL)}
      dt<-fread(f)
      dt<-dt[eval(SIG_FILTER)]
      dt[, `:=`(window = w, window_num = as.integer(gsub("0_","",w)), population = pop, pair = paste(disease_a, disease_b, sep = "_"))]
      dt
    }))
  }))
  pers_by_pop<-all_cumul_sex[, .(n_w = uniqueN(window)), by = .(pair, population)]
  pers_by_pop<-pers_by_pop[n_w == 5]
  slopes_sex<-all_cumul_sex[pers_by_pop, on = .(pair, population)][, {
    x <-window_num
    y <-RR_shrunk
    se<-SE_RR_shrunk
    if (any(is.na(se)) || any(se <= 0) || length(x) < 3) {
      .(slope = NA_real_)
    } else {
      w_fit<-1 / se^2
      fit <-lm(y ~ x, weights = w_fit)
      .(slope = coef(fit)["x"])
    }
  }, by = .(pair, population)]
  slopes_sex[, catA := catname[cate[gsub("_.+","",pair)]]]
  sex_gradient<-slopes_sex[!is.na(catA) & !is.na(slope), .(med_slope = median(slope)), by = .(catA, population)]
  sex_gradient_wide<-dcast(sex_gradient, catA ~ population, value.var = "med_slope")
  if ("women" %in% names(sex_gradient_wide) & "men"   %in% names(sex_gradient_wide)) {sex_gradient_wide[, diff_wm := women - men]}
  ## Sex conditional patterns ##
  pat_by_cat_sex<-rbindlist(lapply(c("women","men"), function(pop) {
    all_cont_sex<-rbindlist(lapply(windows_cont, function(w) {
      f<-sprintf("Results/shrinkage_rr_events/RR_contingency_%s_%s_shrunk.txt", pop, w)
      if (!file.exists(f)) return(NULL)
      dt<-fread(f)
      dt<-dt[!is.na(RR_shrunk) & is.finite(log_RR_shrunk)]
      dt[, `:=`(window = w, window_num = as.integer(gsub("_.+","",w)), pair = paste(disease_a, disease_b, sep = "_"), sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
      dt[, catA := catname[cate[disease_a]]]
      dt
    }))
    p4<-all_cont_sex[, .(n_w = uniqueN(window)), by = pair][n_w == 4, pair]
    cw<-dcast(all_cont_sex[pair %in% p4, .(pair, window, sig, RR_shrunk, catA)], pair + catA ~ window, value.var = c("sig","RR_shrunk"))
    cw[, pattern := fcase(
      (sig_1_2 == TRUE | sig_2_3 == TRUE) & (sig_3_4 == FALSE & sig_4_5 == FALSE), "Early only",
      (sig_1_2 == TRUE | sig_2_3 == TRUE) & (sig_3_4 == TRUE  | sig_4_5 == TRUE), "Persistent conditional risk",
      (sig_1_2 == FALSE & sig_2_3 == FALSE) & (sig_3_4 == TRUE  | sig_4_5 == TRUE), "Late emerging risk", rep(TRUE, .N), "Not significant")]
    result<-cw[!is.na(catA) & !is.na(pattern), .(pct_persist = mean(pattern == "Persistent conditional risk") * 100), by = catA]
    result[, population := pop]
    result
  }))
  
  #### STEP 7: Figures ####
  ## Thresholds for quadrant lines
  x_mid<-median(clock_2d$med_slope)
  y_mid<-median(clock_2d$pct_persist)
  clock_2d[, quadrant := fcase(
    med_slope <= x_mid & pct_persist >= y_mid, "Biphasic\n(fast decay + persistent risk)",
    med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive\n(slow decay + persistent risk)",
    med_slope <= x_mid & pct_persist <  y_mid, "Purely episodic\n(fast decay + no persistence)",
    rep(TRUE, .N), "Chronic stable\n(slow decay + no persistence)")]
  ## Short labels for plotting
  clock_2d[, catA_short := fcase(catA == "Endocrine/Metabolic", "Endocrine/\nMetabolic", catA == "Musculoskeletal", "Musculo-\nskeletal", catA == "Blood/Immune", "Blood/\nImmune", catA == "Genitourinary", "Genito-\nurinary", rep(TRUE, .N), catA)]
  traj_by_cat[, catA := factor(catA, levels = traj_by_cat[order(med_slope), catA])]
  
  p1d<-ggplot(traj_by_cat, aes(x = med_slope, y = catA, color = group, size = n_total)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.85) +
    scale_color_manual(values = c("Fast attenuation\n(episodic/reactive)" = "#1B6CA8", "Moderate attenuation\n(mixed)" = "#1B6CA8", "Slow attenuation\n(chronic/progressive)" = "#1B6CA8")) +
    scale_size_continuous(range = c(3, 10), name = "N comorbidities") +
    labs(
      x = "Median RR slope across follow-up windows (RR units/year)", y = NULL, color = "Attenuation pattern",
      title = "The biological clock of multimorbidity",
      subtitle = paste0( "Rate of change in comorbidity strength (0-1 to 0-5 years) by ICD-10 category.\n", "Negative slope = association weakens over time.")
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "right", panel.grid.minor = element_blank(), panel.grid.major.y = element_line(color = "grey92"), plot.title = element_text(face = "bold"), axis.text.y = element_text(size = 11))
  ggsave("ManuscriptFiles/Plots/Biological_clock_1D.pdf", p1d, width = 11, height = 8, useDingbats = FALSE)
  cat("  Saved: Biological_clock_1D.pdf\n")
  
  p2d<-ggplot(clock_2d, aes(x = med_slope, y = pct_persist, color = quadrant, size = n_total)) +
    annotate("rect", xmin=-Inf, xmax=x_mid, ymin=y_mid, ymax=Inf, fill="#FFF3CD", alpha=0.4) +
    annotate("rect", xmin=x_mid, xmax=Inf, ymin=y_mid, ymax=Inf, fill="#D4EDDA", alpha=0.4) +
    annotate("rect", xmin=-Inf, xmax=x_mid, ymin=-Inf, ymax=y_mid, fill="#F8D7DA", alpha=0.4) +
    annotate("rect", xmin=x_mid, xmax=Inf, ymin=-Inf, ymax=y_mid, fill="#CCE5FF", alpha=0.4) +
    geom_vline(xintercept=x_mid, linetype="dashed", color="grey50", linewidth=0.5) +
    geom_hline(yintercept=y_mid, linetype="dashed", color="grey50", linewidth=0.5) +
    geom_point(alpha=0.85) +
    geom_text_repel(aes(label=catA_short), size=6, color="grey20", max.overlaps=20, box.padding=0.4) +
    scale_color_manual(values = c("Biphasic\n(fast decay + persistent risk)"  = "#E8861A", "Chronic progressive\n(slow decay + persistent risk)" = "#5A9E6F",
      "Purely episodic\n(fast decay + no persistence)" = "#B30000", "Chronic stable\n(slow decay + no persistence)" = "#1B6CA8"
    )) +
    scale_size_continuous(range=c(4,12), name="N comorbidities") +
    annotate("text", x=x_mid-0.003, y=Inf, label="", hjust=1, vjust=1.5, color="grey40", size=4.5, fontface="italic") +
    annotate("text", x=x_mid+0.001, y=Inf, label="", hjust=0, vjust=1.5, color="grey40", size=4.5, fontface="italic") +
    annotate("text", x=-Inf, y=y_mid+0.3, label="", hjust=-0.1, vjust=0, color="grey40", size=4.5, fontface="italic") +
    annotate("text", x=-Inf, y=y_mid-0.3, label="", hjust=-0.1, vjust=1, color="grey40", size=4.5, fontface="italic") +
    labs(
      x = "Rate of RR attenuation (median slope, cumulative windows)", y = "% of comorbidities with persistent conditional risk (1-2 to 4-5 yr)", color = "Biological clock pattern",
      title = "The two-dimensional biological clock of multimorbidity",
      subtitle = paste0("X-axis: how fast comorbidity strength decays over cumulative follow-up.\n", "Y-axis: proportion of disease pairs with persistent risk ", "even if not detected in early windows.\n", "Each point = one ICD-10 category. Size = number of comorbidity pairs.")
    ) +
    theme_minimal(base_size=15) +
    theme(legend.position = "right", panel.grid.minor = element_blank(),
      plot.title = element_text(face="bold", size=15), plot.subtitle = element_text(color="grey40", size=10),
      axis.title = element_text(size=13), axis.text = element_text(size=12),
      legend.text = element_text(size=11), legend.title = element_text(size=12)
    )
  ggsave("ManuscriptFiles/Plots/Biological_clock_2D2_B.pdf", p2d, width=13, height=9, useDingbats=FALSE)
  
  sex_2d<-merge(sex_gradient[population != "both", .(catA, population, med_slope)], pat_by_cat_sex, by = c("catA","population"))
  both_2d<-clock_2d[, .(catA, med_slope_both = med_slope, pct_persist_both = pct_persist)]
  sex_2d<-merge(sex_2d, both_2d, by = "catA")
  sex_wide<-dcast(sex_2d, catA ~ population, value.var = c("med_slope","pct_persist"))
  p_sex<-ggplot() +
    geom_vline(xintercept=x_mid, linetype="dashed", color="grey60", linewidth=0.4) +
    geom_hline(yintercept=y_mid, linetype="dashed", color="grey60", linewidth=0.4) +
    geom_point(data=clock_2d, aes(x=med_slope, y=pct_persist), color="grey75", size=3, alpha=0.5) +
    geom_segment(data=sex_wide, aes(x=med_slope_women, y=pct_persist_women, xend=med_slope_men, yend=pct_persist_men), arrow=arrow(length=unit(0.2,"cm"), type="closed"), color="grey30", linewidth=0.6, alpha=0.7, na.rm=TRUE) +
    geom_point(data=sex_2d, aes(x=med_slope, y=pct_persist, color=population), size=3, alpha=0.85) +
    geom_text_repel(data=sex_2d[population=="women"], aes(x=med_slope, y=pct_persist, label=catA), size=6, color="grey20", max.overlaps=15) +
    scale_color_manual(values=c("women"="#A94C40","men"="#5474A5"), labels=c("women"="Women","men"="Men")) +
    labs(
      x = "Rate of RR attenuation (median slope)", y = "% comorbidities with persistent conditional risk", color = "Sex",
      title = "Sex differences in the biological clock of multimorbidity",
      subtitle = "Arrows point from women to men. Grey = combined population."
    ) +
    theme_minimal(base_size=15) +
    theme(plot.title=element_text(face="bold", hjust=0.5, size=15), plot.subtitle = element_text(hjust=0.5, size=10),
          axis.title = element_text(size=13), axis.text = element_text(size=12),
          legend.text = element_text(size=12), legend.title = element_text(size=13))
  ggsave("ManuscriptFiles/Plots/Biological_clock_2D_sex_B.pdf", p_sex, width=12, height=9, useDingbats=FALSE)
  
  plot_comorbidity_clock<-function(pair_id, title_label) {
    da<-gsub("_.+","",pair_id)
    db<-gsub(".+_","",pair_id)
    cumul_ex<-rbindlist(lapply(windows_cumul, function(w) {
      f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
      if(!file.exists(f)){return(NULL)}
      dt<-fread(f)
      dt<-dt[disease_a == da & disease_b == db]
      if(nrow(dt) == 0){return(NULL)}
      dt[, `:=`(window_num = as.integer(gsub("0_","",w)), type = "Cumulative (0 to k years)")]
      dt
    }))
    cond_ex<-rbindlist(lapply(windows_cont, function(w) {
      f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
      if(!file.exists(f)){return(NULL)}
      dt<-fread(f)
      dt<-dt[disease_a == da & disease_b == db]
      if (nrow(dt) == 0){return(NULL)}
      dt[, `:=`(window_num = as.integer(gsub("_.+","",w)) + 1L, type = "Conditional (year k-1 to k)")]
      dt
    }))
    plot_dt<-rbindlist(list(cumul_ex, cond_ex), fill=TRUE)
    if(nrow(plot_dt) == 0){return(NULL)}
    ## Significance: full criteria including cases_event >= 100 ##
    plot_dt[, sig := lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    ## Separate data for lines (always drawn) and points (shaped by sig) ##
    ggplot(plot_dt, aes(x = window_num, y = RR_shrunk, color = type, fill = type)) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "grey60") +
      geom_ribbon(aes(ymin = CI_low_RR_shrunk, ymax = CI_high_RR_shrunk), alpha = 0.15, color = NA, na.rm = TRUE) +
      geom_line(aes(linetype = type), linewidth = 1.2, na.rm = TRUE) +
      geom_point(aes(shape = sig), size = 3, na.rm = TRUE) +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), labels = c("TRUE" = "Significant", "FALSE" = "Not significant"), name = "") +
      scale_linetype_manual(values = c("Cumulative (0 to k years)" = "solid", "Conditional (year k-1 to k)"  = "dashed"), guide = "none") +
      scale_color_manual(values = c("Cumulative (0 to k years)" = "#1B6CA8", "Conditional (year k-1 to k)" = "#B30000"), name = "Window type") +
      scale_fill_manual(values = c("Cumulative (0 to k years)" = "#1B6CA8", "Conditional (year k-1 to k)" = "#B30000"), guide = "none") +
      scale_x_continuous(breaks = 1:5, labels = c("Year 1", "Year 2", "Year 3", "Year 4", "Year 5")) +
      labs( x = "Follow-up endpoint (years since index diagnosis)", y = "Relative Risk (shrunk)", title = title_label, subtitle = sprintf("%s \u2192 %s", da, db)) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face = "bold",hjust=0.5), plot.subtitle = element_text(hjust=0.5))
  }
  p_ex1<-plot_comorbidity_clock("A63_A64", "Purely episodic\n(Anogenital warts \u2192 Other venereal diseases)")
  p_ex2<-plot_comorbidity_clock("F10_K70", "Chronic stable\n(Alcohol disorder \u2192 Liver cirrhosis)")
  p_ex3<-plot_comorbidity_clock("J01_J32", "Chronic progressive\n(Acute sinusitis \u2192 Chronic sinusitis)")
  p_ex4<-plot_comorbidity_clock("R06_J96", "Late emerging\n(Dyspnoea \u2192 Chronic respiratory failure)")
  p_examples<-(p_ex1 | p_ex2) / (p_ex3 | p_ex4) +
    plot_annotation(
      title = "Examples of the four biological clock patterns",
      subtitle = paste0(
        "Blue = cumulative RR (0 to endpoint year, solid line).  ",
        "Red = conditional RR (year k-1 to k given no prior B, dashed line).  ",
        "Filled circles = significant (lfsr < 0.05, CI\u2093 \u2265 1.01, \u2265100 cases).\n",
        "Note: the biphasic pattern is a category-level property; ",
        "individual pairs may show late-emerging conditional risk ",
        "rather than a strictly biphasic trajectory."
      ),
      theme = theme(plot.title = element_text(face = "bold", size = 18), plot.subtitle = element_text(size = 13, color = "grey40"))
    )
  ggsave("ManuscriptFiles/Plots/Biological_clock_examples_B.pdf", p_examples, width = 14, height = 10, useDingbats = FALSE)
  
  #### Directionality within the biological clock ####
  ## Load networks with directionality already calculated ##
  load_net<-function(pop, w){
    f<-sprintf("Results/networks/RR_net_%s_%s.txt", pop, w)
    if (!file.exists(f)) return(NULL)
    dt<-fread(f)
    dt<-dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    dt[, `:=`(population = pop, window = w, pair = paste(disease_a, disease_b, sep = "_"), pair_rev = paste(disease_b, disease_a, sep = "_"))]
    dt[, catA := catname[cate[disease_a]]]
    dt[, catB := catname[cate[disease_b]]]
    dt
  }
  net_01<-load_net("both", "0_1")
  net_05<-load_net("both", "0_5")
  ## For each A->B pair in 0-1 with preferred direction A->B,
  ## check whether in 0-5 the B->A pair has preferred direction B->A
  dir_01<-net_01[preferred_direction == TRUE, .(pair, pair_rev, catA, catB, theta_01 = theta, n_AB_01 = cases_event)]
  dir_05_rev<-net_05[preferred_direction == TRUE, .(pair_rev = pair, theta_05_rev = theta, n_BA_05 = cases_event)]
  ## Merge: A->B preferred in 0–1 + B->A preferred in 0–5
  reversals_window<-merge(dir_01, dir_05_rev, by = "pair_rev")
  cat(sprintf("  Pairs with directional reversal (0-1 -> 0-5): %d\n", nrow(reversals_window)))
  ## Add the RR of both pairs in 0-5 to keep context ##
  rr_05<-net_05[, .(pair, RR_05 = RR_shrunk)]
  reversals_window<-merge(reversals_window, rr_05, by = "pair", all.x = TRUE)
  reversals_window[, catA := catname[cate[gsub("_.+","",pair)]]]
  reversals_window[, catB := catname[cate[gsub(".+_","",pair)]]]
  ## Enrichment by category ##
  reversal_by_cat<-reversals_window[!is.na(catA), .(n_reversals = .N), by = catA][order(-n_reversals)]
  ## As a percentage of pairs with preferred direction in 0-1 ##
  dir_01_totals<-dir_01[!is.na(catA), .(n_directed = .N), by = catA]
  reversal_by_cat<-merge(reversal_by_cat, dir_01_totals, by = "catA", all.x = TRUE)
  reversal_by_cat[, pct_reversed := n_reversals / n_directed * 100]
  reversal_by_cat[order(-pct_reversed)]
  fwrite(reversals_window, "ManuscriptFiles/Results/Directional_reversals_window_01_vs_05.txt", sep = "\t", quote = FALSE)
  fwrite(reversal_by_cat, "ManuscriptFiles/Results/Directional_reversals_by_category.txt", sep = "\t", quote = FALSE)
  
  net_05_w<-load_net("women", "0_5")
  net_05_m<-load_net("men",   "0_5")
  if(!is.null(net_05_w) && !is.null(net_05_m)){
    ## Pairs with preferred direction A->B in women ##
    dir_w<-net_05_w[preferred_direction == TRUE, .(pair, pair_rev, theta_women = theta, n_AB_women = cases_event)]
    ## Pairs with preferred direction A->B in men ##
    dir_m_rev<-net_05_m[preferred_direction == TRUE, .(pair_rev = pair, theta_men_rev = theta, n_BA_men = cases_event)]
    ## A->B preferred in women + B->A preferred in men ##
    reversals_sex<-merge(dir_w, dir_m_rev, by = "pair_rev")
    reversals_sex[, catA := catname[cate[gsub("_.+","",pair)]]]
    reversals_sex[, catB := catname[cate[gsub(".+_","",pair)]]]
    cat(sprintf("  Pairs with sex-based directional reversal (W: A->B, M: B->A): %d\n", nrow(reversals_sex)))
    ## Top reversals por fuerza de theta (diferencia entre mujeres y hombres)
    reversals_sex[, theta_diff := theta_women - (1 - theta_men_rev)]
    top_sex_reversals<-reversals_sex[order(-abs(theta_diff))][1:50]
    fwrite(reversals_sex, "ManuscriptFiles/Results/Directional_reversals_sex_women_vs_men_0_5.txt", sep = "\t", quote = FALSE)
    ## Enrichment by category
    reversal_sex_by_cat<-reversals_sex[!is.na(catA), .(n_reversals = .N), by = catA][order(-n_reversals)]
    ## As a percentage of pairs with preferred direction in women ##
    dir_w_totals<-dir_w[, .(pair, catA = catname[cate[gsub("_.+","",pair)]])]
    dir_w_totals<-dir_w_totals[!is.na(catA), .(n_directed_w = .N), by = catA]
    reversal_sex_by_cat<-merge(reversal_sex_by_cat, dir_w_totals, by = "catA", all.x = TRUE)
    reversal_sex_by_cat[, pct_reversed := n_reversals / n_directed_w * 100]
    fwrite(reversal_sex_by_cat, "ManuscriptFiles/Results/Directional_reversals_sex_by_category.txt", sep = "\t", quote = FALSE)
    
    ## Add reversal information to clock_2d ##
    reversal_sex_clock<-merge(clock_2d[, .(catA, med_slope, pct_persist, quadrant)], reversal_sex_by_cat[, .(catA, n_reversals, pct_reversed, n_directed_w)], by = "catA", all.x = TRUE)
    ## Distinguish between "0 reversals but directed pairs exist" and "no data" ##
    reversal_sex_clock[is.na(n_directed_w), data_available := FALSE]
    reversal_sex_clock[!is.na(n_directed_w), data_available := TRUE]
    reversal_sex_clock[is.na(pct_reversed), pct_reversed := 0]
    p_rev_sex<-ggplot(reversal_sex_clock, aes(x = med_slope, y = pct_persist)) +
      annotate("rect", xmin=-Inf, xmax=x_mid, ymin=y_mid, ymax=Inf, fill="#FFF3CD", alpha=0.3) +
      annotate("rect", xmin=x_mid, xmax=Inf, ymin=y_mid, ymax=Inf, fill="#D4EDDA", alpha=0.3) +
      annotate("rect", xmin=-Inf, xmax=x_mid, ymin=-Inf, ymax=y_mid, fill="#F8D7DA", alpha=0.3) +
      annotate("rect", xmin=x_mid, xmax=Inf, ymin=-Inf, ymax=y_mid, fill="#CCE5FF", alpha=0.3) +
      geom_vline(xintercept=x_mid, linetype="dashed", color="grey50", linewidth=0.4) +
      geom_hline(yintercept=y_mid, linetype="dashed", color="grey50", linewidth=0.4) +
      geom_point(data = reversal_sex_clock[data_available == TRUE], aes(color = quadrant, size = pct_reversed), alpha = 0.85) +
      geom_point(data = reversal_sex_clock[data_available == FALSE], color = "grey70", size = 3, shape = 1, stroke = 1) +
      ## Labels >1% with percentages ##
      geom_text_repel(data = reversal_sex_clock[data_available == TRUE & pct_reversed > 1], aes(label = paste0(catA, "\n(", round(pct_reversed, 1), "%)")), size = 3.2, color = "grey20", max.overlaps = 20, box.padding = 0.5) +
      geom_text_repel(data = reversal_sex_clock[data_available == TRUE & pct_reversed <= 1], aes(label = catA), size = 2.5, color = "grey50", max.overlaps = 20, box.padding = 0.3) +
      geom_text_repel(data = reversal_sex_clock[data_available == FALSE], aes(label = catA), size = 2.5, color = "grey70", max.overlaps = 20, fontface = "italic", box.padding = 0.3) +
      scale_color_manual(values = c(
        "Biphasic\n(fast decay + persistent risk)" = "#E8861A",
        "Chronic progressive\n(slow decay + persistent risk)"= "#5A9E6F",
        "Purely episodic\n(fast decay + no persistence)" = "#B30000",
        "Chronic stable\n(slow decay + no persistence)" = "#1B6CA8"
      )) +
      scale_size_continuous(range = c(2, 14), name = "% pairs with\nsex reversal", breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5)) +
      labs(
        x = "Rate of RR attenuation (median slope)", y = "% comorbidities with persistent conditional risk", color = "Biological clock pattern",
        title = "Sex-based directional reversals within the biological clock",
        subtitle = paste0("Point size = % of directed pairs where A->B is preferred in women ", "but B->A is preferred in men (window 0-5).\n", "Labels show % reversal for categories >1%. ", "Grey open circles = no directed pairs in women.")
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
    ggsave("ManuscriptFiles/Plots/Biological_clock_directional_reversals_sex.pdf", p_rev_sex, width = 13, height = 9, useDingbats = FALSE)
  }
  
  #### Biological clock example table - Supplementary Table ####
  clock_examples<-list(
    "A63→A64\n(Episodic)\nAnogenital warts→Other venereal diseases" = "A63_A64",
    "F10→K70\n(Chronic stable)\nAlcohol disorder→Cirrhosis" = "F10_K70",
    "J01→J32\n(Chronic progressive)\nAcute→Chronic sinusitis" = "J01_J32",
    "R06→J96\n(Biphasic)\nDyspnoea→Chronic resp. failure" = "R06_J96"
  )
  ## Window labels ##
  window_labels<-c("0_1" = "Cumulative 0–1 yr", "0_2" = "Cumulative 0–2 yr", "0_3" = "Cumulative 0–3 yr", "0_4" = "Cumulative 0–4 yr", "0_5" = "Cumulative 0–5 yr",
    "1_2" = "Conditional 1–2 yr", "2_3" = "Conditional 2–3 yr", "3_4" = "Conditional 3–4 yr", "4_5" = "Conditional 4–5 yr"
  )
  all_windows<-c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
  ## Load data for all pairs and windows ##
  supp_table<-rbindlist(lapply(all_windows, function(w){
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if(!file.exists(f)){return(NULL)}
    dt<-fread(f)
    dt[, pair := paste(disease_a, disease_b, sep="_")]
    dt[pair %in% unlist(clock_examples), .(pair, window = w, RR = RR_shrunk, CI_low = CI_low_RR_shrunk, CI_high = CI_high_RR_shrunk, cases = cases_event, sig = lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)]
  }))
  ## Format RR string: "X.XX (X.XX–X.XX)*" where * = significant ##
  supp_table[, RR_fmt := sprintf("%s%.2f (%.2f–%.2f)",
    ifelse(sig, "", ""),
    RR, CI_low, CI_high
  )]
  supp_table[, cell := sprintf("%s\n[n=%d]%s", RR_fmt, cases,
                               ifelse(sig, "*", ""))]
  ## Add window label and order ##
  supp_table[, window_label := window_labels[window]]
  supp_table[, window_order := match(window, all_windows)]
  ## Add pair label ##
  pair_to_label<-setNames(names(clock_examples), unlist(clock_examples))
  supp_table[, pair_label := pair_to_label[pair]]
  ## Wide format: rows = windows, cols = pairs ##
  supp_wide<-as.data.table(dcast(supp_table[order(window_order)], window_order + window_label ~ pair_label, value.var = "cell"))
  supp_wide[, window_order := NULL]
  setnames(supp_wide, "window_label", "Window")
  ## Reorder columns ##
  col_order<-c("Window", names(clock_examples))
  setcolorder(supp_wide, col_order)
  ## Add section separator rows ##
  cumul_rows<-supp_wide[1:5]
  cond_rows <-supp_wide[6:9]
  ## Print to check ##
  print(supp_wide)
  
  ## Save as formatted docx table ##
  # Clean column names for display
  col_names_clean<-c("Window", "A63→A64\n(Episodic)\nAnogenital warts→Other venereal diseases", "F10→K70\n(Chronic stable)\nAlcohol disorder→Cirrhosis", "J01→J32\n(Chronic progressive)\nAcute→Chronic sinusitis", "R06→J96\n(Biphasic)\nDyspnoea→Chronic resp. failure")
  ft<-flextable(supp_wide) %>%
    set_header_labels(Window = "Window", .list = setNames(as.list(col_names_clean[-1]), names(clock_examples))) %>%
    add_header_row(values = c("", "Purely episodic", "Chronic stable", "Chronic progressive", "Biphasic"), colwidths = c(1, 1, 1, 1, 1)) %>%
    bg(i = 1:5, bg = "#F0F4F8") %>%
    bg(i = 6:9, bg = "#FFF8F0") %>%
    bold(i = c(1, 6), bold = FALSE) %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "body") %>%
    width(j = 1, width = 1.5) %>%
    width(j = 2:5, width = 1.8) %>%
    border_outer(part = "all", border = fp_border(color = "grey40", width = 1)) %>%
    border_inner_h(part = "body", border = fp_border(color = "grey70", width = 0.5)) %>%
    border_inner_v(part = "all", border = fp_border(color = "grey70", width = 0.5)) %>%
    # Add footnote
    add_footer_lines("RR = posterior shrunk relative risk; 95% CI = credible interval. 
     * Significant: lfsr < 0.05, lower bound CI ≥ 1.01, ≥100 cases. 
     n = number of cases (index disease patients developing secondary disease).
     Blue shading = cumulative windows (0 to k years); 
     Orange shading = conditional windows (year k-1 to k, given no prior secondary disease)."
    ) %>%
    fontsize(size = 8, part = "footer") %>%
    color(color = "grey40", part = "footer")
  ## Save ##
  doc<-read_docx() %>%
    body_add_par("Supplementary Table S_BiologicalClock_Examples. Temporal trajectories of relative risk for the four paradigmatic examples of the biological clock patterns.", style = "heading 2") %>%
    body_add_flextable(ft)
  print(doc, target = "ManuscriptFiles/Results/Supp_Table_BiologicalClock_Examples.docx")
}

#### Biological clock level 2 ####
if(args[1]=="BiologicalClock_level2"){
  if(!dir.exists("ManuscriptFiles/Plots/Subcategory")){dir.create("ManuscriptFiles/Plots/Subcategory", recursive=TRUE)}
  if(!dir.exists("ManuscriptFiles/Results/Subcategory")){dir.create("ManuscriptFiles/Results/Subcategory", recursive=TRUE)}
  ## ICD-10 category mapping (same as main script) ##
  code<-c(paste0("A0",0:9), paste0("A",10:99),
            paste0("B0",0:9), paste0("B",10:99),
            paste0("C0",0:9), paste0("C",10:99),
            paste0("D0",0:9), paste0("D",10:48), paste0("D",50:89),
            paste0("E0",0:9), paste0("E",10:99),
            paste0("F0",0:9), paste0("F",10:99),
            paste0("G0",0:9), paste0("G",10:99),
            paste0("H0",0:9), paste0("H",10:59), paste0("H",60:95),
            paste0("I0",0:9), paste0("I",10:99),
            paste0("J0",0:9), paste0("J",10:99),
            paste0("K0",0:9), paste0("K",10:93),
            paste0("L0",0:9), paste0("L",10:99),
            paste0("M0",0:9), paste0("M",10:99),
            paste0("N0",0:9), paste0("N",10:99),
            paste0("O0",0:9), paste0("O",10:99),
            paste0("P0",0:9), paste0("P",10:96),
            paste0("Q0",0:9), paste0("Q",10:99),
            paste0("R0",0:9), paste0("R",10:99),
            paste0("S0",0:9), paste0("S",10:99),
            paste0("T0",0:9), paste0("T",10:98))
  cate<-c(rep("I",200), rep("II",149), rep("III",40),
            rep("IV",100), rep("V",100), rep("VI",100),
            rep("VII",60), rep("VIII",36), rep("IX",100),
            rep("X",100), rep("XI",94), rep("XII",100),
            rep("XIII",100), rep("XIV",100), rep("XV",100),
            rep("XVI",97), rep("XVII",100), rep("XVIII",100),
            rep("XIX",199))
  catname<-c("Infectious","Neoplasms","Blood/Immune",
               "Endocrine/Metabolic","Mental","Nervous",
               "Eye","Ear","Circulatory","Respiratory",
               "Digestive","Skin","Musculoskeletal",
               "Genitourinary","Pregnancy","Perinatal",
               "Congenital","Symptoms","Injury")
  names(catname)<-unique(cate)
  names(cate)<-code
  SIG_FILTER<-quote(lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100)
  windows_cumul<-c("0_1","0_2","0_3","0_4","0_5")
  windows_cont <-c("1_2","2_3","3_4","4_5")
  # Load subcategory mapping
  subcat_map<-fread("./icd10_3digitos_categoria_subcategoria_who2016.csv")
  subcat_map<-subcat_map[, .(
    icd10_3 = as.character(icd10_3),
    subcategoria_rango = as.character(subcategoria_rango),
    subcategoria_nombre = as.character(subcategoria_nombre),
    categoria_nombre = as.character(categoria_nombre)
  )]
  # Named vectors for lookup
  subcat_rango<-setNames(subcat_map$subcategoria_rango, subcat_map$icd10_3)
  subcat_nombre<-setNames(subcat_map$subcategoria_nombre, subcat_map$icd10_3)
  cat_nombre<-setNames(subcat_map$categoria_nombre, subcat_map$icd10_3)
  cat(sprintf("Subcategory map loaded: %d codes, %d subcategories, %d categories\n", nrow(subcat_map), uniqueN(subcat_map$subcategoria_rango), uniqueN(subcat_map$categoria_nombre)))
  
  #### Load cumulative windows and identify persistent pairs ####
  all_cumul<-rbindlist(lapply(windows_cumul, function(w) {
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) { cat("  Missing:", f, "\n"); return(NULL) }
    dt<-fread(f)
    dt<-dt[eval(SIG_FILTER)]
    dt[, `:=`(window = w, window_num = as.integer(gsub("0_","",w)), pair = paste(disease_a, disease_b, sep="_"))]
    dt[, subcatA := subcat_rango[disease_a]]
    dt[, subcatA_nombre := subcat_nombre[disease_a]]
    dt[, catA    := catname[cate[disease_a]]]
    dt
  }))
  pairs_per_window<-all_cumul[, .(n_windows=uniqueN(window)), by=pair]
  persistent_pairs<-pairs_per_window[n_windows==5, pair]
  cat(sprintf("Pairs with data in >=1 window: %d\n", nrow(pairs_per_window)))
  cat(sprintf("Persistent pairs (all 5 windows): %d\n", length(persistent_pairs)))
  # Check subcategory coverage
  n_with_subcat<-all_cumul[pair %in% persistent_pairs & !is.na(subcatA), uniqueN(pair)]
  cat(sprintf("Persistent pairs with subcategory: %d (%.1f%%)\n", n_with_subcat, n_with_subcat/length(persistent_pairs)*100))
  
  #### Compute weighted slope per persistent pair ####
  rr_traj<-all_cumul[pair %in% persistent_pairs, .(pair, window_num, RR_shrunk, SE_RR_shrunk, disease_a, disease_b, subcatA, subcatA_nombre, catA)]
  slopes<-rr_traj[, {
    x<-window_num
    y<-RR_shrunk
    se<-SE_RR_shrunk
    if(any(is.na(se)) || any(se<=0) || length(x)<3){
      .(slope=NA_real_, slope_p=NA_real_)
    } else {
      w_fit<-1/se^2
      fit<-lm(y ~ x, weights=w_fit)
      sm<-summary(fit)$coefficients
      .(slope = sm["x","Estimate"], slope_p = sm["x","Pr(>|t|)"])
    }
  }, by=.(pair, subcatA, subcatA_nombre, catA)]
  slopes[!is.na(slope_p), slope_fdr := p.adjust(slope_p, method="BH")]
  slopes[is.na(slope_p),  slope_fdr := NA_real_]
  
  #### Subcategory-level slope summary ####
  traj_by_subcat<-slopes[!is.na(subcatA) & !is.na(slope), .(n_total = .N, med_slope = median(slope), mean_slope = mean(slope), sd_slope = sd(slope), catA = first(catA), subcatA_nombre = first(subcatA_nombre)), by=subcatA][order(med_slope)]
  cat(sprintf("Subcategories with slope calculated: %d\n", nrow(traj_by_subcat)))
  cat(sprintf("Subcategories with >=10 pairs: %d\n", sum(traj_by_subcat$n_total >= 10)))
  ## Filter to subcategories with enough pairs for stable estimates ##
  MIN_PAIRS_SUBCAT<-10
  traj_by_subcat_filt<-traj_by_subcat[n_total >= MIN_PAIRS_SUBCAT]
  cat(sprintf("Subcategories after filtering (>=%d pairs): %d\n", MIN_PAIRS_SUBCAT, nrow(traj_by_subcat_filt)))
  ## Save ##
  fwrite(traj_by_subcat, "ManuscriptFiles/Results/Subcategory/Clock_slopes_by_subcategory.txt", sep="\t", quote=FALSE)
  
  #### Conditional windows at subcategory level ####
  colcod<-c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C","#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54","#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
  names(colcod)<-unique(cate)
  # Map from categoria_nombre (long) to color
  cat_colors<-setNames(as.character(colcod), catname)
  
  ## Also map from short catname used in subcat_map$categoria_nombre ##
  ## Check which names are used in subcat_map ##
  cat("Category names in subcat_map:\n")
  print(unique(subcat_map$categoria_nombre))
  ## Filter out pairs with NA subcatA before dcast ##
  all_cont_filt<-all_cont[pair %in% pairs_cont4 & !is.na(subcatA), .(pair, window, sig, RR_shrunk, subcatA, subcatA_nombre, catA)]
  cat(sprintf("Pairs with subcatA assigned: %d of %d\n", uniqueN(all_cont_filt$pair), length(pairs_cont4)))
  cont_wide<-dcast(all_cont_filt, pair + subcatA + subcatA_nombre + catA ~ window, value.var = c("sig","RR_shrunk"))
  ## Verify column names ##
  cat("Columns after dcast:", paste(names(cont_wide), collapse=", "), "\n")
  ## Check expected columns exist ##
  expected_sig<-c("sig_1_2","sig_2_3","sig_3_4","sig_4_5")
  expected_rr <-c("RR_shrunk_1_2","RR_shrunk_2_3","RR_shrunk_3_4","RR_shrunk_4_5")
  ## Classify pattern ##
  cont_wide[, pattern := fcase(
    (sig_1_2==TRUE | sig_2_3==TRUE) & (sig_3_4==FALSE & sig_4_5==FALSE), "Early only (risk fades)",
    (sig_1_2==TRUE | sig_2_3==TRUE) & (sig_3_4==TRUE  | sig_4_5==TRUE), "Persistent conditional risk",
    (sig_1_2==FALSE & sig_2_3==FALSE) & (sig_3_4==TRUE | sig_4_5==TRUE), "Late emerging risk", rep(TRUE,.N), "Not significant"
  )]
  print(cont_wide[, .N, by=pattern][order(-N)])
  ## Subcategory-level conditional summary ##
  pat_by_subcat<-cont_wide[!is.na(subcatA) & !is.na(pattern), .(
    n_total = .N,
    pct_early = mean(pattern=="Early only (risk fades)")*100,
    pct_persist = mean(pattern=="Persistent conditional risk")*100,
    pct_late = mean(pattern=="Late emerging risk")*100,
    pct_notsig = mean(pattern=="Not significant")*100,
    catA = first(catA),
    subcatA_nombre = first(subcatA_nombre)
  ), by=subcatA][order(-pct_persist)]
  print(head(pat_by_subcat, 10))
  fwrite(pat_by_subcat, "ManuscriptFiles/Results/Subcategory/Clock_conditional_by_subcategory.txt", sep="\t", quote=FALSE)
  #### 2D clock with correct colors ####
  clock_2d_subcat<-merge(traj_by_subcat_filt[, .(subcatA, subcatA_nombre, catA, med_slope, n_total)], pat_by_subcat[, .(subcatA, pct_persist, pct_late)], by="subcatA")
  x_mid_sub<-median(clock_2d_subcat$med_slope)
  y_mid_sub<-median(clock_2d_subcat$pct_persist)
  cat(sprintf("x_mid_sub: %.5f\n", x_mid_sub))
  cat(sprintf("y_mid_sub: %.3f\n",  y_mid_sub))
  clock_2d_subcat[, quadrant := fcase(
    med_slope <= x_mid_sub & pct_persist >= y_mid_sub, "Late-emerging\n(fast decay + persistent risk)",
    med_slope >  x_mid_sub & pct_persist >= y_mid_sub, "Chronic progressive\n(slow decay + persistent risk)",
    med_slope <= x_mid_sub & pct_persist <  y_mid_sub, "Purely episodic\n(fast decay + no persistence)",
    rep(TRUE,.N), "Chronic stable\n(slow decay + no persistence)"
  )]
  ## Check what category names look like in clock_2d_subcat ##
  cat("\nCategory names in clock_2d_subcat:\n")
  print(unique(clock_2d_subcat$catA))
  ## Build color mapping matching actual catA names in data ##
  cat_colors_sub<-setNames(
    as.character(colcod[names(colcod) %in% names(catname)[catname %in% unique(clock_2d_subcat$catA)]]),
    catname[catname %in% unique(clock_2d_subcat$catA)]
  )
  ## Check overlap ##
  cat("\nMatched categories:", length(cat_colors_sub), "\n")
  unmatched<-setdiff(unique(clock_2d_subcat$catA), names(cat_colors_sub))
  if (length(unmatched) > 0) {
    cat("Unmatched categories:", paste(unmatched, collapse=", "), "\n")
    ## Add grey for unmatched ##
    extra_colors<-setNames(rep("grey50", length(unmatched)), unmatched)
    cat_colors_sub<-c(cat_colors_sub, extra_colors)
  }
  
  p2d_sub<-ggplot(clock_2d_subcat, aes(x=med_slope, y=pct_persist, color=catA, size=n_total)) +
    annotate("rect", xmin=-Inf, xmax=x_mid_sub, ymin=y_mid_sub, ymax=Inf,   fill="#FFF3CD", alpha=0.4) +
    annotate("rect", xmin=x_mid_sub, xmax=Inf, ymin=y_mid_sub, ymax=Inf,   fill="#D4EDDA", alpha=0.4) +
    annotate("rect", xmin=-Inf, xmax=x_mid_sub, ymin=-Inf, ymax=y_mid_sub,  fill="#F8D7DA", alpha=0.4) +
    annotate("rect", xmin=x_mid_sub, xmax=Inf, ymin=-Inf, ymax=y_mid_sub,  fill="#CCE5FF", alpha=0.4) +
    geom_vline(xintercept=x_mid_sub, linetype="dashed", color="grey50", linewidth=0.5) +
    geom_hline(yintercept=y_mid_sub, linetype="dashed", color="grey50", linewidth=0.5) +
    geom_point(alpha=0.8) +
    geom_text_repel(aes(label=subcatA), size=2.3, color="grey20", max.overlaps=40, box.padding=0.3, segment.size=0.3, min.segment.length=0.2) +
    scale_color_manual(values=cat_colors_sub, name="ICD-10 category") +
    scale_size_continuous(range=c(2,10), name="N persistent pairs") +
    annotate("text", x=x_mid_sub-0.002, y=Inf, label="Fast decay", hjust=1, vjust=1.5, color="grey40", size=3.5, fontface="italic") +
    annotate("text", x=x_mid_sub+0.001, y=Inf, label="Slow decay", hjust=0, vjust=1.5, color="grey40", size=3.5, fontface="italic") +
    annotate("text", x=-Inf, y=y_mid_sub+0.2, label="High persistent risk", hjust=-0.05, vjust=0, color="grey40", size=3.2, fontface="italic") +
    annotate("text", x=-Inf, y=y_mid_sub-0.2, label="Low persistent risk", hjust=-0.05, vjust=1, color="grey40", size=3.2, fontface="italic") +
    labs(x = "Rate of RR attenuation (median slope, cumulative windows)", y = "% comorbidities with persistent conditional risk",
      title = "Two-dimensional biological clock at ICD-10 subcategory level",
      subtitle = paste0("Each point = one ICD-10 subcategory (>=", MIN_PAIRS_SUBCAT, " persistent pairs). Colour = ICD-10 chapter.")
    ) +
    theme_minimal(base_size=12) +
    theme(legend.position = "right", panel.grid.minor = element_blank(), plot.title = element_text(face="bold"), plot.subtitle = element_text(color="grey40", size=9))
  ggsave("ManuscriptFiles/Plots/Subcategory/Clock_2D_subcategory.pdf", p2d_sub, width=16, height=11, useDingbats=FALSE)
  
  ## Save results ##
  fwrite(clock_2d_subcat[order(catA, med_slope)], "ManuscriptFiles/Results/Subcategory/Clock_2D_subcategory.txt", sep="\t", quote=FALSE)
  ## Within-category heterogeneity ##
  heterogeneity_by_cat<-clock_2d_subcat[, .(
    n_subcategories = .N,
    slope_range = round(max(med_slope) - min(med_slope), 4),
    slope_min = round(min(med_slope), 4),
    slope_max = round(max(med_slope), 4),
    persist_range = round(max(pct_persist) - min(pct_persist), 1),
    persist_min = round(min(pct_persist), 1),
    persist_max = round(max(pct_persist), 1),
    n_quadrants = uniqueN(quadrant),
    quadrants = paste(sort(unique(gsub("\n.*","",quadrant))), collapse=", ")
  ), by=catA][order(-slope_range)]
  print(heterogeneity_by_cat)
  fwrite(heterogeneity_by_cat, "ManuscriptFiles/Results/Subcategory/Clock_within_category_heterogeneity.txt", sep="\t", quote=FALSE)
}

#### Biological clock validation ####
if(args[1]=="Validate_Biological_Clock"){
  #### Pair-level bootstrap ####
  set.seed(1)
  N_BOOT<-1000
  ## Pre-merge slopes and conditional pattern for efficiency ##
  slopes_pat<-merge(slopes[!is.na(slope) & !is.na(catA), .(pair, catA, slope)], cont_wide[!is.na(catA), .(pair, pattern)], by="pair", all.x=TRUE)
  slopes_pat[is.na(pattern), pattern := "Not significant"]
  ## Bootstrap function for one replicate ##
  boot_one<-function(dt){dt[, .SD[sample(.N, .N, replace=TRUE)], by=catA][, .(med_slope = median(slope, na.rm=TRUE), pct_persist = mean(pattern=="Persistent conditional risk", na.rm=TRUE)*100), by=catA]}
  ## Run bootstrap ##
  boot_results<-rbindlist(lapply(seq_len(N_BOOT), function(i){
    res<-boot_one(slopes_pat)
    res[, replicate := i]
    res
  }))
  cat(sprintf("Bootstrap complete: %d replicates x %d categories\n", N_BOOT, uniqueN(boot_results$catA)))
  ## Compute bootstrap CIs and quadrant retention ##
  boot_summary<-boot_results[, {
    ## Assign quadrant in each replicate ##
    quad<-fcase(med_slope <= x_mid & pct_persist >= y_mid, "Late-emerging", med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive", med_slope <= x_mid & pct_persist <  y_mid, "Purely episodic", rep(TRUE, .N), "Chronic stable")
    ## Original quadrant ##
    orig_quad<-clock_2d[catA == .BY$catA, gsub("\n.*","", quadrant)]
    if (length(orig_quad) == 0) orig_quad<-NA_character_
    .(
      slope_mean = mean(med_slope),
      slope_lo = quantile(med_slope, 0.025),
      slope_hi = quantile(med_slope, 0.975),
      persist_mean = mean(pct_persist),
      persist_lo = quantile(pct_persist, 0.025),
      persist_hi = quantile(pct_persist, 0.975),
      quadrant_retention = mean(quad == orig_quad, na.rm=TRUE) * 100,
      original_quadrant  = orig_quad[1]
    )
  }, by=catA]
  ## Print BOOTSTRAP SUMMARY ##
  print(boot_summary[order(-quadrant_retention),.(catA, slope_lo = round(slope_lo, 4), slope_hi = round(slope_hi, 4), persist_lo = round(persist_lo, 1), persist_hi = round(persist_hi, 1), quadrant_retention = round(quadrant_retention, 1), original_quadrant)])
  fwrite(boot_summary, "ManuscriptFiles/Results/Biological_clock_bootstrap_summary.txt", sep="\t", quote=FALSE)
  
  #### Leave-one-index-disease-out sensitivity ####
  ## For each disease in each category, remove all pairs where disease_a == that disease and recompute category position ##
  all_diseases_by_cat<-slopes_pat[!is.na(catA), .(diseases = list(unique(gsub("_.+","",pair)))), by=catA]
  loo_results<-rbindlist(lapply(seq_len(nrow(all_diseases_by_cat)), function(i) {
    cat_i <-all_diseases_by_cat$catA[i]
    dis_i <-all_diseases_by_cat$diseases[[i]]
    rbindlist(lapply(dis_i, function(d) {
      dt_sub<-slopes_pat[catA == cat_i & gsub("_.+","",pair) != d]
      if (nrow(dt_sub) < 5) return(NULL)
      data.table(catA  = cat_i, left_out = d, n_remaining = nrow(dt_sub), med_slope = median(dt_sub$slope, na.rm=TRUE), pct_persist = mean(dt_sub$pattern=="Persistent conditional risk", na.rm=TRUE) * 100)
    }))
  }))
  ## Assign quadrant for each LOO replicate ##
  loo_results[, quadrant_loo := fcase(med_slope <= x_mid & pct_persist >= y_mid, "Late-emerging", med_slope >  x_mid & pct_persist >= y_mid, "Chronic progressive", med_slope <= x_mid & pct_persist <  y_mid, "Purely episodic", rep(TRUE, .N), "Chronic stable")]
  ## Merge with original quadrant ##
  loo_results<-merge(loo_results, clock_2d[, .(catA, original_quadrant = gsub("\n.*","", quadrant))], by="catA")
  ## Stability: % of LOO replicates retaining original quadrant ##
  loo_summary<-loo_results[, .(
    n_diseases = .N, quadrant_retention = mean(quadrant_loo == original_quadrant)*100,
    slope_range = max(med_slope) - min(med_slope), persist_range = max(pct_persist) - min(pct_persist),
    original_quadrant = first(original_quadrant), most_influential = left_out[which.max(abs(med_slope - median(med_slope)))]
  ), by=catA][order(-quadrant_retention)]
  
  ## Print LOO summary ##
  print(loo_summary[, .(catA, n_diseases, quadrant_retention = round(quadrant_retention, 1), slope_range = round(slope_range, 4), persist_range = round(persist_range, 1), original_quadrant, most_influential)])
  fwrite(loo_summary, "ManuscriptFiles/Results/Biological_clock_LOO_summary.txt", sep="\t", quote=FALSE)
  
  #### Alternative quadrant definitions ####
  ## K-means with k=4 on normalized coordinates ##
  clock_norm<-copy(clock_2d)
  clock_norm[, slope_z := scale(med_slope)[,1]]
  clock_norm[, persist_z := scale(pct_persist)[,1]]
  set.seed(1)
  km<-kmeans(clock_norm[, .(slope_z, persist_z)], centers=4, nstart=50, iter.max=100)
  clock_norm[, cluster_kmeans := km$cluster]
  ## Match clusters to quadrant names by centroid position ##
  centers<-as.data.table(km$centers)
  centers[, cluster := 1:.N]
  centers[, quadrant_km := fcase(slope_z <= 0 & persist_z >= 0, "Late-emerging", slope_z >  0 & persist_z >= 0, "Chronic progressive", slope_z <= 0 & persist_z <  0, "Purely episodic",rep(TRUE, .N), "Chronic stable")]
  clock_norm<-merge(clock_norm, centers[, .(cluster, quadrant_km)], by.x="cluster_kmeans", by.y="cluster")
  ## Concordance between median-split and k-means ##
  concordance<-merge(clock_norm[, .(catA, quadrant_km)], clock_2d[, .(catA, quadrant_orig = gsub("\n.*","", quadrant))], by="catA")
  concordance[, same := quadrant_km == quadrant_orig]
  ## Print MEDIAN SPLIT vs K-MEANS CONCORDANCE ##
  print(concordance[, .(catA, quadrant_orig, quadrant_km, same)])
  cat(sprintf("Overall concordance: %d/%d categories (%.1f%%)\n", sum(concordance$same), nrow(concordance), mean(concordance$same)*100))
  fwrite(concordance, "ManuscriptFiles/Results/Biological_clock_quadrant_sensitivity.txt", sep="\t", quote=FALSE)
  
  #### STEP 4: Combined supplementary table ####
  robustness_table<-merge(
    boot_summary[, .(catA, original_quadrant, slope_CI = sprintf("(%.4f, %.4f)", round(slope_lo,4), round(slope_hi,4)), persist_CI  = sprintf("(%.1f, %.1f)", round(persist_lo,1), round(persist_hi,1)), boot_retention = round(quadrant_retention, 1))],
    loo_summary[, .(catA, loo_retention = round(quadrant_retention, 1), most_influential)],
    by="catA"
  )
  robustness_table<-merge(
    robustness_table,
    concordance[, .(catA, quadrant_kmeans = quadrant_km, kmeans_agrees = same)],
    by="catA"
  )
  ## Print full robustness table ##
  print(robustness_table[order(-boot_retention)])
  fwrite(robustness_table, "ManuscriptFiles/Results/Biological_clock_robustness_table.txt", sep="\t", quote=FALSE)
  ## Bootstrap function corrected ##
  boot_one_v2<-function(dt) {
    res<-dt[, .SD[sample(.N, .N, replace=TRUE)], by=catA][, .(med_slope   = median(slope, na.rm=TRUE), pct_persist = mean(pattern=="Persistent conditional risk", na.rm=TRUE)*100), by=catA]
    ## Recalculate thresholds in each replicate ##
    x_mid_b<-median(res$med_slope)
    y_mid_b<-median(res$pct_persist)
    res[, quadrant_b := fcase(med_slope <= x_mid_b & pct_persist >= y_mid_b, "Late-emerging", med_slope >  x_mid_b & pct_persist >= y_mid_b, "Chronic progressive", med_slope <= x_mid_b & pct_persist < y_mid_b, "Purely episodic", rep(TRUE, .N), "Chronic stable")]
    res
  }
  ## Run bootstrap v2 ##
  boot_results_v2<-rbindlist(lapply(seq_len(N_BOOT), function(i) {
    res<-boot_one_v2(slopes_pat)
    res[, replicate := i]
    res
  }))
  
  ## Re-calculate retention with thresholds re-calculated in each replica ##
  boot_summary_v2<-boot_results_v2[, {
    orig_quad<-clock_2d[catA == .BY$catA, gsub("\n.*","", quadrant)]
    if (length(orig_quad) == 0) orig_quad<-NA_character_
    .(slope_lo = quantile(med_slope, 0.025), slope_hi = quantile(med_slope, 0.975), persist_lo = quantile(pct_persist, 0.025), persist_hi = quantile(pct_persist, 0.975), quadrant_retention = mean(quadrant_b == orig_quad, na.rm=TRUE) * 100, original_quadrant = orig_quad[1])
  }, by=catA]
  print(boot_summary_v2[order(-quadrant_retention), .(catA, slope_CI  = sprintf("(%.4f, %.4f)", round(slope_lo,4), round(slope_hi,4)), persist_CI = sprintf("(%.1f, %.1f)", round(persist_lo,1), round(persist_hi,1)), quadrant_retention = round(quadrant_retention, 1), original_quadrant)])
  ## Rename areas ##
  clock_2d[, quadrant_clean := fcase(grepl("Biphasic|Late-emerging", quadrant), "Late-emerging", grepl("Chronic progressive", quadrant), "Chronic progressive", grepl("Chronic stable", quadrant), "Chronic stable", grepl("episodic|Episodic", quadrant), "Episodic")]
  quadrant_colors<-c("Late-emerging" = "#E8861A", "Chronic progressive" = "#5A9E6F", "Episodic" = "#B30000", "Chronic stable" = "#1B6CA8")
  boot_results_v2[, quadrant_clean := gsub("Biphasic", "Late-emerging", gsub("\n.*","", clock_2d$quadrant[match(catA, clock_2d$catA)]))]
  boot_results_v2[, quadrant_clean := fcase(grepl("Biphasic|Late-emerging", clock_2d$quadrant[match(catA, clock_2d$catA)]), "Late-emerging", grepl("Chronic progressive", clock_2d$quadrant[match(catA, clock_2d$catA)]), "Chronic progressive", grepl("Chronic stable", clock_2d$quadrant[match(catA, clock_2d$catA)]), "Chronic stable", rep(TRUE, .N), "Episodic")]
  
  p_boot<-ggplot() +
    annotate("rect", xmin=-Inf, xmax=x_mid, ymin=y_mid, ymax=Inf, fill="#FFF3CD", alpha=0.3) +
    annotate("rect", xmin=x_mid, xmax=Inf,  ymin=y_mid, ymax=Inf, fill="#D4EDDA", alpha=0.3) +
    annotate("rect", xmin=-Inf, xmax=x_mid, ymin=-Inf,  ymax=y_mid, fill="#F8D7DA", alpha=0.3) +
    annotate("rect", xmin=x_mid, xmax=Inf,  ymin=-Inf,  ymax=y_mid, fill="#CCE5FF", alpha=0.3) +
    geom_vline(xintercept=x_mid, linetype="dashed", color="grey50", linewidth=0.4) +
    geom_hline(yintercept=y_mid, linetype="dashed", color="grey50", linewidth=0.4) +
    geom_point(data=boot_results_v2[replicate <= 200], aes(x=med_slope, y=pct_persist, color=quadrant_clean), alpha=0.05, size=0.8) +
    stat_ellipse(data=boot_results_v2, aes(x=med_slope, y=pct_persist, color=quadrant_clean), level=0.95, linewidth=0.8) +
    geom_point(data=merge(clock_2d, boot_summary_v2[, .(catA, quadrant_retention)], by="catA"), aes(x=med_slope, y=pct_persist, color=quadrant_clean), size=4, shape=16) +
    geom_text_repel(
      data=merge(clock_2d, boot_summary_v2[, .(catA, quadrant_retention)], by="catA"),
      aes(x=med_slope, y=pct_persist, label=sprintf("%s\n(%.0f%%)", catA_short, quadrant_retention)),
      size=2.8, color="grey20", max.overlaps=20, box.padding=0.3
    ) +
    scale_color_manual(values=quadrant_colors, name="Biological clock pattern") +
    labs(
      x = "Rate of RR attenuation (median slope, cumulative windows)", y = "% comorbidities with persistent conditional risk",
      title = "Stability of the biological clock under pair-level bootstrap",
      subtitle = paste0("Points = original category estimates. ", "Ellipses = 95% bootstrap CI (1,000 replicates, ", "pair-level resampling within each ICD-10 chapter).\n", "Percentages = quadrant retention rate across bootstrap replicates.")
    ) +
    theme_minimal(base_size=12) +
    theme( legend.position = "bottom", panel.grid.minor = element_blank(), plot.title = element_text(face="bold"), plot.subtitle = element_text(color="grey40", size=9))
  ggsave("ManuscriptFiles/Plots/Supp_Biological_clock_bootstrap.pdf", p_boot, width=13, height=10, useDingbats=FALSE)
  
  ## Regenerar robustness_table ##
  robustness_table_v2<-merge(
    boot_summary_v2[, .(catA, original_quadrant,
                        slope_CI = sprintf("(%.4f, %.4f)", round(slope_lo,4), round(slope_hi,4)),
                        persist_CI  = sprintf("(%.1f, %.1f)", round(persist_lo,1), round(persist_hi,1)),
                        boot_retention = round(quadrant_retention, 1))],
    loo_summary[, .(catA, loo_retention = round(quadrant_retention, 1), most_influential)], by="catA"
  )
  robustness_table_v2<-merge(robustness_table_v2, concordance[, .(catA, quadrant_kmeans = quadrant_km, kmeans_agrees = same)], by="catA")
  robustness_table_v2[, original_quadrant := gsub("Purely episodic", "Episodic", gsub("Biphasic", "Late-emerging", original_quadrant))]
  robustness_table_v2[, quadrant_kmeans := gsub("Purely episodic", "Episodic", gsub("Biphasic", "Late-emerging", quadrant_kmeans))]
  print(robustness_table_v2[order(-boot_retention)])
  fwrite(robustness_table_v2, "ManuscriptFiles/Results/Biological_clock_robustness_table_v2.txt", sep="\t", quote=FALSE)
  fwrite(robustness_table_v2[order(-boot_retention)], "ManuscriptFiles/Results/Biological_clock_robustness_table_v2.txt", sep="\t", quote=FALSE)
}

#### Look for overlaps between women and men networks, changes in their evolution over time, and differences in directionality ####
if(args[1]=="Sex_window_differences"){
  #### Network overlap by sex and time window ####
  overlap_sex_windows<-rbindlist(lapply(c("0_1","0_2","0_3","0_4","0_5"), function(w) {
    f_w<-sprintf("Results/networks/RR_net_women_%s.txt", w)
    f_m<-sprintf("Results/networks/RR_net_men_%s.txt", w)
    if(!file.exists(f_w) || !file.exists(f_m)){return(NULL)}
    net_w<-fread(f_w)[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    net_m<-fread(f_m)[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    pairs_w<-paste(net_w$disease_a, net_w$disease_b, sep="_")
    pairs_m<-paste(net_m$disease_a, net_m$disease_b, sep="_")
    shared<-intersect(pairs_w, pairs_m)
    only_w<-setdiff(pairs_w, pairs_m)
    only_m<-setdiff(pairs_m, pairs_w)
    data.table(window = w, n_women = length(pairs_w), n_men = length(pairs_m), n_shared = length(shared), n_only_w = length(only_w), n_only_m = length(only_m),
      jaccard = round(length(shared)/length(union(pairs_w, pairs_m)), 3),
      pct_w_in_m = round(length(shared)/length(pairs_w) * 100, 1),
      pct_m_in_w = round(length(shared)/length(pairs_m) * 100, 1)
    )
  }))
  ## Network overlap by sex and window ##
  print(overlap_sex_windows)
  
  #### Pairs with preferred directionality by sex and window ####
  dir_sex_windows<-rbindlist(lapply(c("0_1","0_2","0_3","0_4","0_5"), function(w) {
    f_w<-sprintf("Results/networks/RR_net_women_%s.txt", w)
    f_m<-sprintf("Results/networks/RR_net_men_%s.txt", w)
    if(!file.exists(f_w) || !file.exists(f_m)){return(NULL)}
    net_w<-fread(f_w)[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100 & preferred_direction == TRUE]
    net_m<-fread(f_m)[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100 & preferred_direction == TRUE]
    pairs_w<-paste(net_w$disease_a, net_w$disease_b, sep="_")
    pairs_m<-paste(net_m$disease_a, net_m$disease_b, sep="_")
    ## Pairs with concordant directionality ##
    shared_dir<-intersect(pairs_w, pairs_m)
    ## Reversions: A->B in women, B->A in men ##
    pairs_w_rev<-paste(net_w$disease_b, net_w$disease_a, sep="_")
    reversals<-intersect(pairs_w_rev, pairs_m)
    data.table(
      window = w, n_dir_women = length(pairs_w), n_dir_men = length(pairs_m), n_concordant = length(shared_dir), n_reversals = length(reversals),
      pct_concordant = round(length(shared_dir)/length(union(pairs_w, pairs_m)) * 100, 1),
      pct_reversed   = round(length(reversals)/length(union(pairs_w, pairs_m)) * 100, 1)
    )
  }))
  ## Print directionality by sex and window ##
  print(dir_sex_windows)
  
  #### Directionality stability: 0–1 vs. each time window ####
  f_01_w<-"Results/networks/RR_net_women_0_1.txt"
  f_01_m<-"Results/networks/RR_net_men_0_1.txt"
  net_01_w<-fread(f_01_w)[preferred_direction == TRUE & lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  net_01_m<-fread(f_01_m)[preferred_direction == TRUE & lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  stability_dir<-rbindlist(lapply(c("0_2","0_3","0_4","0_5"), function(w) {
    f_w<-sprintf("Results/networks/RR_net_women_%s.txt", w)
    f_m<-sprintf("Results/networks/RR_net_men_%s.txt", w)
    if(!file.exists(f_w) || !file.exists(f_m)){return(NULL)}
    net_w<-fread(f_w)[preferred_direction == TRUE & lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    net_m<-fread(f_m)[preferred_direction == TRUE & lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    ## Women: pairs directed in 0-1 that reverse in window w ##
    pairs_01_w<-paste(net_01_w$disease_a, net_01_w$disease_b, sep="_")
    pairs_w<-paste(net_w$disease_a, net_w$disease_b, sep="_")
    rev_w_name<-paste(net_01_w$disease_b, net_01_w$disease_a, sep="_")
    reversed_w<-sum(rev_w_name %in% pairs_w)
    stable_w<-sum(pairs_01_w %in% pairs_w)
    ## Men ##
    pairs_01_m<-paste(net_01_m$disease_a, net_01_m$disease_b, sep="_")
    pairs_m<-paste(net_m$disease_a, net_m$disease_b, sep="_")
    rev_m_name<-paste(net_01_m$disease_b, net_01_m$disease_a, sep="_")
    reversed_m<-sum(rev_m_name %in% pairs_m)
    stable_m<-sum(pairs_01_m %in% pairs_m)
    
    data.table(window = w, 
      n_dir_01_w = length(pairs_01_w), n_stable_w = stable_w, n_reversed_w = reversed_w, pct_stable_w = round(stable_w/length(pairs_01_w) * 100, 1), pct_reversed_w = round(reversed_w/length(pairs_01_w) * 100, 1),
      n_dir_01_m = length(pairs_01_m), n_stable_m = stable_m, n_reversed_m = reversed_m, pct_stable_m = round(stable_m/length(pairs_01_m)*100, 1), pct_reversed_m = round(reversed_m/length(pairs_01_m)*100, 1)
    )
  }))
  ## Print directional stability 0-1 vs. other windows ##
  print(stability_dir)
  ## Load 0-5 networks by sex ##
  net_w_05<-fread("Results/networks/RR_net_women_0_5.txt")[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  net_m_05<-fread("Results/networks/RR_net_men_0_5.txt")[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  net_w_05[, pair := paste(disease_a, disease_b, sep="_")]
  net_m_05[, pair := paste(disease_a, disease_b, sep="_")]
  ## Add categories ##
  net_w_05[, catA := catname[cate[disease_a]]]
  net_w_05[, catB := catname[cate[disease_b]]]
  net_m_05[, catA := catname[cate[disease_a]]]
  net_m_05[, catB := catname[cate[disease_b]]]
  ## Classify pairs: shared, only_women, only_men ##
  pairs_w<-net_w_05$pair
  pairs_m<-net_m_05$pair
  shared<-intersect(pairs_w, pairs_m)
  net_w_05[, status := ifelse(pair %in% shared, "shared", "only_women")]
  net_m_05[, status := ifelse(pair %in% shared, "shared", "only_men")]
  ## Combine ##
  net_sex<-rbindlist(list(net_w_05[, .(pair, catA, catB, status, sex = "women")], net_m_05[!pair %in% shared, .(pair, catA, catB, status, sex = "men")]))
  cat("Shared:", length(shared), "\n")
  cat("Only women:", sum(net_sex$status == "only_women"), "\n")
  cat("Only men:", sum(net_sex$status == "only_men"), "\n")
  
  #### Figures ####
  ## PANEL A: Enrichment catAxcatB en shared vs only_women+only_men ##
  universe_AB<-net_sex[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
  shared_AB<-net_sex[status == "shared" & !is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
  setnames(shared_AB, "N", "n_shared")
  enrich_shared<-merge(universe_AB, shared_AB, by = c("catA","catB"), all.x = TRUE)
  enrich_shared[is.na(n_shared), n_shared := 0L]
  total_shared<-length(shared)
  total_universe<-nrow(net_sex[!is.na(catA) & !is.na(catB)])
  enrich_shared<-enrich_shared[, {
    a1<-n_shared
    b1<-total_shared - a1
    a0<-N - a1
    b0<-(total_universe - total_shared) - a0
    if (any(c(a1,b1,a0,b0) < 0)) {
      .(odds = NA_real_, p = NA_real_)
    } else {
      ft<-fisher.test(matrix(c(a1,b1,a0,b0), nrow=2, byrow=TRUE))
      .(odds = as.numeric(ft$estimate), p = ft$p.value)
    }
  }, by = .(catA, catB, N, n_shared)]
  enrich_shared[, p_adj := p.adjust(p, "BH")]
  enrich_shared[, log_odds := log(odds)]
  enrich_shared[, sig := !is.na(p_adj) & p_adj <= 0.05]
  enrich_shared[, fill_value := fifelse(sig, log_odds, 0)]
  cat("TOP ENRICHED IN SHARED")
  print(enrich_shared[sig==TRUE & odds>1][order(-odds)][1:10, .(catA, catB, n_shared, N, odds=round(odds,2), p_adj=round(p_adj,4))])
  cat("TOP DEPLETED IN SHARED")
  print(enrich_shared[sig==TRUE & odds<1][order(odds)][1:10, .(catA, catB, n_shared, N, odds=round(odds,2), p_adj=round(p_adj,4))])
  
  ## PANEL B: Enrichment catAxcatB en women vs men ##
  women_AB<-net_w_05[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
  men_AB  <-net_m_05[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB)]
  setnames(women_AB, "N", "n_women")
  setnames(men_AB, "N", "n_men")
  enrich_sex<-merge(women_AB, men_AB, by = c("catA","catB"), all = TRUE)
  enrich_sex[is.na(n_women), n_women := 0L]
  enrich_sex[is.na(n_men), n_men := 0L]
  total_women<-nrow(net_w_05[!is.na(catA) & !is.na(catB)])
  total_men<-nrow(net_m_05[!is.na(catA) & !is.na(catB)])
  enrich_sex<-enrich_sex[, {
    a1<-n_women
    b1<-total_women - a1
    a0<-n_men
    b0<-total_men - a0
    if (any(c(a1,b1,a0,b0) < 0)) {
      .(odds = NA_real_, p = NA_real_)
    } else {
      ft<-fisher.test(matrix(c(a1,b1,a0,b0), nrow=2, byrow=TRUE))
      .(odds = as.numeric(ft$estimate), p = ft$p.value)
    }
  }, by = .(catA, catB, n_women, n_men)]
  enrich_sex[, p_adj := p.adjust(p, "BH")]
  enrich_sex[, log_odds := log(odds)]
  enrich_sex[, sig := !is.na(p_adj) & p_adj <= 0.05]
  enrich_sex[, fill_value := fifelse(sig, log_odds, 0)]
  cat("TOP ENRICHED IN WOMEN")
  print(enrich_sex[sig==TRUE & odds>1][order(-odds)][1:10,.(catA, catB, n_women, n_men, odds=round(odds,2), p_adj=round(p_adj,4))])
  cat("TOP ENRICHED IN MEN")
  print(enrich_sex[sig==TRUE & odds<1][order(odds)][1:10,.(catA, catB, n_women, n_men, odds=round(odds,2), p_adj=round(p_adj,4))])
  
  #### Plot them ####
  abbrev<-c(
    "Infectious" = "Infect.", "Neoplasms" = "Neopl.", "Blood/Immune" = "Blood/Imm.", "Endocrine/Metabolic" = "Endocr.", "Mental" = "Mental",
    "Nervous" = "Nervous", "Eye" = "Eye", "Ear" = "Ear", "Circulatory" = "Circ.", "Respiratory" = "Resp.", "Digestive" = "Digest.", "Skin" = "Skin",
    "Musculoskeletal" = "Muscul.", "Genitourinary" = "Genito.", "Pregnancy" = "Preg.", "Perinatal" = "Perinat.", "Congenital" = "Congenit.", "Symptoms" = "Sympt.","Injury" = "Injury"
  )
  plot_enrich_heatmap<-function(dt, title_label,low_col = "#C53030", high_col = "#2B6CB0",exclude_cats = NULL) {
    dt2<-copy(dt)
    if (!is.null(exclude_cats)){dt2<-dt2[!catA %in% exclude_cats & !catB %in% exclude_cats]}
    dt2[, catA_s := abbrev[catA]]
    dt2[, catB_s := abbrev[catB]]
    catA_levels<-sort(unique(dt2$catA_s))
    catB_levels<-sort(unique(dt2$catB_s))
    grid<-CJ(catA_s = catA_levels, catB_s = catB_levels)
    grid<-merge(grid, dt2[, .(catA_s, catB_s, fill_value, sig)], by = c("catA_s","catB_s"), all.x = TRUE)
    grid[is.na(fill_value), fill_value := NA_real_]
    ggplot(grid, aes(x = catB_s, y = catA_s, fill = fill_value)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(data = grid[!is.na(sig) & sig == TRUE & !is.na(fill_value) & abs(fill_value) > 0.5 & is.finite(fill_value)], aes(label = sprintf("%.1f", fill_value)), size = 2.5, color = "grey20") +
      scale_fill_gradient2(midpoint = 0, low = low_col, mid = "white", high = high_col, na.value = "grey93", name = "log(OR)\n(sig only)", limits = c(-3, 3), oob = scales::squish) +
      labs(x = "Secondary disease category (B)", y = "Index disease category (A)", title = title_label) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 9), panel.grid = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5), legend.position = "right")
  }
  
  p_shared<-plot_enrich_heatmap(enrich_shared, "A) Category-pair enrichment among shared vs sex-specific comorbidities\n(blue = overrepresented in shared; red = underrepresented)", exclude_cats = c("Pregnancy", "Perinatal"))
  p_sex<-plot_enrich_heatmap(enrich_sex, "B) Category-pair enrichment in women vs men\n(blue = overrepresented in women; red = overrepresented in men)", exclude_cats = c("Pregnancy", "Perinatal"))
  p_combined<-p_shared/p_sex +
    plot_annotation(title = "Sex differences in comorbidity network structure (window 0-5 years)", theme = theme(plot.title = element_text(face = "bold", size = 13,hjust = 0.5)))
  ggsave("ManuscriptFiles/Plots/Sex_network_enrichment_heatmaps.pdf", p_combined, width = 14, height = 20, useDingbats = FALSE)
  
  ## Examples of directional reversions women vs. men ##
  net_w_full<-fread("Results/networks/RR_net_women_0_5.txt")[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  net_m_full<-fread("Results/networks/RR_net_men_0_5.txt")[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  net_w_full[, pair := paste(disease_a, disease_b, sep="_")]
  net_w_full[, pair_rev := paste(disease_b, disease_a, sep="_")]
  net_m_full[, pair := paste(disease_a, disease_b, sep="_")]
  ## Pairs with preferred directionality in women ##
  dir_w<-net_w_full[preferred_direction == TRUE, .(pair, pair_rev, disease_a, disease_b, catA = catname[cate[disease_a]], catB = catname[cate[disease_b]], theta_w = theta, RR_w = RR_shrunk, n_AB_w = cases_event)]
  ## Pairs with preferred directionality in men ##
  dir_m_rev<-net_m_full[preferred_direction == TRUE, .(pair_rev = pair, theta_m_rev = theta, RR_m_rev = RR_shrunk, n_BA_m = cases_event)]
  ## Reversions: A->B in women, B->A in men ##
  reversals_sex_examples<-merge(dir_w, dir_m_rev, by = "pair_rev")
  reversals_sex_examples[, theta_diff := theta_w - (1 - theta_m_rev)]
  ## Top examples by RR in women ##
  print(reversals_sex_examples[order(-RR_w)][1:20, .(pair, catA, catB, RR_w = round(RR_w, 2), theta_w = round(theta_w, 2), RR_m_rev = round(RR_m_rev, 2), theta_m = round(1 - theta_m_rev, 2), n_AB_w, n_BA_m)])
  ## Top examples by theta_diff ##
  print(reversals_sex_examples[order(-abs(theta_diff))][1:20, .(pair, catA, catB, RR_w = round(RR_w, 2), theta_w = round(theta_w, 2), RR_m_rev = round(RR_m_rev, 2), theta_m = round(1 - theta_m_rev, 2), theta_diff = round(theta_diff, 2))])
  ## Examples directional reversions 0-1 vs 0-5 ##
  net_01_both<-fread("Results/networks/RR_net_both_0_1.txt")[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100 & preferred_direction == TRUE]
  net_05_both<-fread("Results/networks/RR_net_both_0_5.txt")[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100 & preferred_direction == TRUE]
  net_01_both[, pair := paste(disease_a, disease_b, sep="_")]
  net_01_both[, pair_rev := paste(disease_b, disease_a, sep="_")]
  net_05_both[, pair := paste(disease_a, disease_b, sep="_")]
  net_01_both[, catA := catname[cate[disease_a]]]
  net_01_both[, catB := catname[cate[disease_b]]]
  dir_05_rev<-net_05_both[preferred_direction == TRUE, .(pair_rev = pair, theta_05_rev = theta, RR_05_rev = RR_shrunk)]
  reversals_window_examples<-merge(net_01_both[, .(pair, pair_rev, catA, catB, theta_01 = theta, RR_01 = RR_shrunk)], dir_05_rev, by = "pair_rev")
  ## Top by RR in 0-1 ##
  print(reversals_window_examples[order(-RR_01)][1:20, .(pair, catA, catB, RR_01 = round(RR_01, 2), theta_01 = round(theta_01, 2), RR_05_rev = round(RR_05_rev, 2), theta_05 = round(1 - theta_05_rev, 2))])
  ## Numbers for the manuscript ##
  ## Overlap general 0-5 ##
  cat("GENERAL OVERLAP 0-5")
  cat("Women:", nrow(net_w_05), "\n")
  cat("Men:", nrow(net_m_05), "\n")
  cat("Shared:", length(shared), "\n")
  cat("Only women:", sum(net_w_05$status == "only_women"), "\n")
  cat("Only men:", nrow(net_m_05[!pair %in% shared]), "\n")
  cat("% men in women:", round(length(shared)/nrow(net_m_05)*100,1), "\n")
  cat("% women in men:", round(length(shared)/nrow(net_w_05)*100,1), "\n")
  
  ## Top enriched in shared (Panel A) ##
  print(enrich_shared[sig==TRUE & odds>1 & is.finite(log_odds)][order(-odds)][1:10, .(catA, catB, n_shared, N, odds=round(odds,2), p_adj=round(p_adj,5))])
  print(enrich_shared[sig==TRUE & odds<1 & catA=="Genitourinary"][order(odds)][1:5, .(catA, catB, n_shared, N, odds=round(odds,2))])
  ## Top enriched in women vs men (Panel B) ##
  print(enrich_sex[sig==TRUE & odds>1 & is.finite(log_odds)][order(-odds)][1:10, .(catA, catB, n_women, n_men, odds=round(odds,2), p_adj=round(p_adj,5))])
  print(enrich_sex[sig==TRUE & odds<1 & is.finite(log_odds)][order(odds)][1:10, .(catA, catB, n_women, n_men, odds=round(odds,2), p_adj=round(p_adj,5))])
  ## Top reversions (sex) by RR
  print(reversals_sex_examples[order(-RR_w)][1:10, .(pair, catA, catB, RR_w=round(RR_w,2), theta_w=round(theta_w,2), RR_m_rev=round(RR_m_rev,2), theta_m=round(1-theta_m_rev,2))])
  ## Reversions 0-1 vs 0-5 ##
  rev_win_by_cat<-reversals_window_examples[!is.na(catA), .N, by=catA][order(-N)]
  dir_01_totals<-net_01_both[, .(n_dir=.N), by=.(catA=catname[cate[disease_a]])]
  rev_win_by_cat<-merge(rev_win_by_cat, dir_01_totals, by="catA")
  rev_win_by_cat[, pct := round(N/n_dir*100,1)]
  print(rev_win_by_cat[order(-pct)])
}

#### Check the temporal robustness of comorbidity associations ####
if(args[1]=="Temporal_robustness_of_comorbidities"){
  ## Load all the windows, cumulative and conditional ##
  all_windows_9<-c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
  window_type  <-c(rep("cumulative",5), rep("conditional",4))
  all_9<-rbindlist(lapply(seq_along(all_windows_9), function(i) {
    w<-all_windows_9[i]
    wt<-window_type[i]
    f<-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if(!file.exists(f)) return(NULL)
    dt<-fread(f)
    dt<-dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    dt[, `:=`(window = w, window_type = wt, pair = paste(disease_a, disease_b, sep="_"))]
    dt
  }))
  
  ## For each pair, which window detects it? ##
  pair_windows<-all_9[, .(
    in_0_1 = any(window == "0_1"),
    in_0_2 = any(window == "0_2"),
    in_0_3 = any(window == "0_3"),
    in_0_4 = any(window == "0_4"),
    in_0_5 = any(window == "0_5"),
    in_1_2 = any(window == "1_2"),
    in_2_3 = any(window == "2_3"),
    in_3_4 = any(window == "3_4"),
    in_4_5 = any(window == "4_5"),
    n_cumul = sum(window_type == "cumulative"),
    n_cond  = sum(window_type == "conditional")
  ), by = pair]
  
  dim(pair_windows)
  ## Biological taxonomy ##
  pair_windows[, bio_class := fcase(
    n_cumul == 5 & n_cond == 4, "Omnipresent (all 9 windows)",
    n_cumul == 5 & n_cond > 0, "Persistent cumulative + partial conditional",
    n_cumul == 5 & n_cond == 0, "Persistent cumulative only (no conditional risk)",
    in_0_5 == TRUE & in_0_1 == FALSE & n_cond == 0, "Long-window only (late cumulative, no early signal)",
    in_0_5 == TRUE & in_0_1 == FALSE & n_cond > 0, "Late cumulative + conditional risk",
    n_cumul == 0 & n_cond > 0, "Conditional only (no cumulative signal)",
    in_0_1 == TRUE & in_0_5 == FALSE & n_cond == 0, "Early transient (short window only)",
    in_0_1 == TRUE & in_0_5 == FALSE & n_cond > 0, "Early cumulative + late conditional",
    rep(TRUE, .N), "Mixed/partial"
  )]
  ## Counts by class ##
  class_counts<-pair_windows[, .N, by = bio_class][order(-N)]
  print(class_counts)
  ## Add ICD-10 categories ##
  pair_windows[, catA := catname[cate[gsub("_.+","",pair)]]]
  pair_windows[, catB := catname[cate[gsub(".+_","",pair)]]]
  ## Distribution by category ##
  class_by_cat<-pair_windows[!is.na(catA), .N, by = .(catA, bio_class)]
  class_by_cat[, pct := N / sum(N) * 100, by = catA]
  class_by_cat_wide<-dcast(class_by_cat, catA ~ bio_class, value.var = "pct", fill = 0)
  print(class_by_cat_wide)
  write.table(class_by_cat_wide,"Class_by_type.txt",quote=F,sep="\t",row.names=F)
  
  ## Read the data ##
  class_pct<-fread("Class_by_type.txt")
  ## Rename the columns for the plot ##
  setnames(class_pct, 
           old = c(
             "Omnipresent (all 9 windows)", "Persistent cumulative only (no conditional risk)", "Persistent cumulative + partial conditional",
             "Long-window only (late cumulative, no early signal)", "Late cumulative + conditional risk",
             "Conditional only (no cumulative signal)", "Early transient (short window only)", "Early cumulative + late conditional", "Mixed/partial"),
           new = c("Omnipresent", "Persistent\n(no conditional)", "Persistent +\nconditional", "Long-window\nonly", "Late +\nconditional",
             "Conditional\nonly", "Early\ntransient", "Early +\nlate cond.", "Mixed/\npartial")
  )
  ## Long format for ggplot ##
  class_long<-melt(class_pct,  id.vars = "catA", variable.name = "bio_type", value.name = "pct")
  ## Order categories by omnipresent % ##
  cat_order<-class_pct[order(`Omnipresent`), catA]
  class_long[, catA := factor(catA, levels = cat_order)]
  ## Order classes by importance ##
  type_order<-c("Omnipresent", "Persistent +\nconditional", "Persistent\n(no conditional)", "Long-window\nonly", "Late +\nconditional", "Conditional\nonly", "Early\ntransient", "Early +\nlate cond.", "Mixed/\npartial")
  class_long[, bio_type := factor(bio_type, levels = type_order)]
  ## Plot heatmap ##
  p_heatmap<-ggplot(class_long, aes(x = bio_type, y = catA, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.1f", pct)), size = 2.8, color = "grey20") +
    scale_fill_gradient(low = "white", high = "#1B4521", name = "% of pairs", limits = c(0, 100)) +
    labs(x = "Biological detection pattern", y = "ICD-10 category (index disease)",
      title = "Temporal detection taxonomy of comorbidity pairs",
      subtitle = paste0("Each cell shows the percentage of comorbidity pairs in that category ", "belonging to each detection pattern.\n", "Ordered by % omnipresent (most temporally robust).")
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10), axis.text.y = element_text(size = 11),
      panel.grid = element_blank(), plot.title = element_text(face = "bold"), plot.subtitle = element_text(color = "grey40", size = 9), legend.position = "right")
  
  ggsave("ManuscriptFiles/Plots/Biological_clock_taxonomy_heatmap.pdf", p_heatmap, width = 14, height = 9, useDingbats = FALSE)
  print(p_heatmap)
  
  ## Alternatively, stacked bar ordered by % omnipresent ##
  ## Groups 9 classes into 4 narratives ##
  class_long[, narrative_group := fcase(
    bio_type %in% c("Omnipresent", "Persistent +\nconditional", "Persistent\n(no conditional)"), "Robust\n(detectable early & late)",
    bio_type %in% c("Long-window\nonly", "Late +\nconditional"), "Late-emerging\n(requires long follow-up)",
    bio_type %in% c("Conditional\nonly", "Early\ntransient", "Early +\nlate cond."), "Transient/conditional\n(window-specific)",
    rep(TRUE, .N), "Mixed/partial"
  )]
  ## Sum by narrative group
  narr_summary<-class_long[, .(pct = sum(pct)), by = .(catA, narrative_group)]
  ## Order by robust %
  robust_order<-narr_summary[narrative_group == "Robust\n(detectable early & late)", .(catA, pct_robust = pct)][order(pct_robust), catA]
  narr_summary[, catA := factor(catA, levels = robust_order)]
  narr_summary[, narrative_group := factor(narrative_group, levels = c( "Robust\n(detectable early & late)", "Late-emerging\n(requires long follow-up)", "Transient/conditional\n(window-specific)", "Mixed/partial"))]
  p_stacked<-ggplot(narr_summary, aes(x = pct, y = catA, fill = narrative_group)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = c("Robust\n(detectable early & late)" = "#1B6CA8", "Late-emerging\n(requires long follow-up)" = "#E8861A", "Transient/conditional\n(window-specific)" = "#B30000", "Mixed/partial" = "#CCCCCC"), name = "Detection pattern") +
    geom_vline(xintercept = c(25,50,75), linetype = "dashed", color = "white", linewidth = 0.3) +
    labs(x = "% of comorbidity pairs", y = NULL, title = "Temporal robustness of comorbidities by disease category") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold",hjust = 0.5),
      plot.subtitle = element_text(color = "grey40", size = 10), axis.text.y = element_text(size = 11), legend.position = "bottom")
  ggsave("ManuscriptFiles/Plots/Biological_clock_taxonomy_stacked.pdf", p_stacked, width = 12, height = 9, useDingbats = FALSE)
  print(p_stacked)
}

#### Do the analysis looking for pairs of diseases and not only focusing on the index disease ####
if(args[1] == "Temporal_robustness_of_comorbidities_by_pairs") {
  ## Convert disease diagnoses into disease categories ##
  code<-c(paste("A0",0:9,sep=""),paste("A",10:99,sep=""),paste("B0",0:9,sep=""),paste("B",10:99,sep=""),paste("C0",0:9,sep=""),paste("C",10:99,sep=""),
          paste("D0",0:9,sep=""),paste("D",10:48,sep=""),paste("D",50:89,sep=""),paste("E0",0:9,sep=""),paste("E",10:99,sep=""),paste("F0",0:9,sep=""),
          paste("F",10:99,sep=""),paste("G0",0:9,sep=""),paste("G",10:99,sep=""),paste("H0",0:9,sep=""),paste("H",10:59,sep=""),paste("H",60:95,sep=""),
          paste("I0",0:9,sep=""),paste("I",10:99,sep=""),paste("J0",0:9,sep=""),paste("J",10:99,sep=""),paste("K0",0:9,sep=""),paste("K",10:93,sep=""),
          paste("L0",0:9,sep=""),paste("L",10:99,sep=""),paste("M0",0:9,sep=""),paste("M",10:99,sep=""),paste("N0",0:9,sep=""),paste("N",10:99,sep=""),
          paste("O0",0:9,sep=""),paste("O",10:99,sep=""),paste("P0",0:9,sep=""),paste("P",10:96,sep=""),paste("Q0",0:9,sep=""),paste("Q",10:99,sep=""),
          paste("R0",0:9,sep=""),paste("R",10:99,sep=""),paste("S0",0:9,sep=""),paste("S",10:99,sep=""),paste("T0",0:9,sep=""),paste("T",10:98,sep=""))
  cate<-c(rep("I",200),rep("II",149),rep("III",40),rep("IV",100),rep("V",100),rep("VI",100),rep("VII",60),rep("VIII",36),rep("IX",100),rep("X",100),
          rep("XI",94),rep("XII",100),rep("XIII",100),rep("XIV",100),rep("XV",100),rep("XVI",97),rep("XVII",100),rep("XVIII",100),rep("XIX",199))
  catname<-c("Infectious and parasitic","Neoplasms","Blood and blood-forming organs (immune)",
             "Endocrine, nutritional and metabolic","Mental and behavioural",
             "Nervous system","Eye and adnexa","Ear and mastoid process",
             "Circulatory system","Respiratory system","Digestive system",
             "Skin and subcutaneous tissue","Musculoskeletal system and connective tissue",
             "Genitourinary system","Pregnancy, childbirth and the puerperium",
             "Certain conditions originating in the perinatal period","Congenital malformations and chromosomal abnormalities",
             "Symptoms, signs and abnormal laboratory findings","Injury, poisoning")
  names(catname)<-unique(cate)
  names(cate)<-code
  ## Disease category ##
  catename<-as.character(catname[cate[code]])
  names(catename)<-code
  colcod<-c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C","#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54","#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
  names(colcod)<-unique(cate)
  catcol<-as.character(colcod) ; names(catcol)<-catname
  codcol<-names(colcod) ; names(codcol)<-as.character(colcod)
  distocol<-as.character(colcod[as.character(cate)]) ; names(distocol)<-names(cate)
  discat<-cbind(code,catname[cate[code]],distocol[code])
  colnames(discat)<-c("disease","category","color")
  discat<-as.data.table(discat)
  ## Load all the windows, cumulative and conditional ##
  all_windows_9<-c("0_1","0_2","0_3","0_4","0_5","1_2","2_3","3_4","4_5")
  window_type  <-c(rep("cumulative",5), rep("conditional",4))
  all_9<-rbindlist(lapply(seq_along(all_windows_9), function(i) {
    w <-all_windows_9[i]
    wt<-window_type[i]
    f <-sprintf("Results/shrinkage_rr_events/RR_contingency_both_%s_shrunk.txt", w)
    if (!file.exists(f)) return(NULL)
    dt<-fread(f)
    dt<-dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    dt[, `:=`(window = w, window_type = wt, pair = paste(disease_a, disease_b, sep = "_"))]
    dt
  }))
  ## Detect windows per pair ##
  pair_windows<-all_9[, .(
    in_0_1  = any(window == "0_1"),
    in_0_2  = any(window == "0_2"),
    in_0_3  = any(window == "0_3"),
    in_0_4  = any(window == "0_4"),
    in_0_5  = any(window == "0_5"),
    in_1_2  = any(window == "1_2"),
    in_2_3  = any(window == "2_3"),
    in_3_4  = any(window == "3_4"),
    in_4_5  = any(window == "4_5"),
    n_cumul = sum(window_type == "cumulative"),
    n_cond  = sum(window_type == "conditional")
  ), by = pair]
  ## Biological taxonomy ##
  pair_windows[, bio_class := fcase(
    n_cumul == 5 & n_cond == 4, "Omnipresent (all 9 windows)",
    n_cumul == 5 & n_cond > 0, "Persistent cumulative + partial conditional",
    n_cumul == 5 & n_cond == 0, "Persistent cumulative only (no conditional risk)",
    in_0_5 == TRUE & in_0_1 == FALSE & n_cond == 0, "Long-window only (late cumulative, no early signal)",
    in_0_5 == TRUE & in_0_1 == FALSE & n_cond > 0, "Late cumulative + conditional risk",
    n_cumul == 0 & n_cond > 0, "Conditional only (no cumulative signal)",
    in_0_1 == TRUE & in_0_5 == FALSE & n_cond == 0, "Early transient (short window only)",
    in_0_1 == TRUE & in_0_5 == FALSE & n_cond > 0, "Early cumulative + late conditional",
    rep(TRUE, .N), "Mixed/partial"
  )]
  ## Add ICD-10 categories ##
  pair_windows[, catA := catname[cate[gsub("_.+", "", pair)]]]
  pair_windows[, catB := catname[cate[gsub(".+_", "", pair)]]]
  ## Distribution by catA x catB pair ##
  class_by_catAB<-pair_windows[!is.na(catA) & !is.na(catB), .N, by = .(catA, catB, bio_class)]
  class_by_catAB[, pct := N / sum(N)*100, by = .(catA, catB)]
  class_by_catAB_wide<-dcast(class_by_catAB, catA + catB ~ bio_class, value.var = "pct", fill = 0)
  fwrite(class_by_catAB_wide, "ManuscriptFiles/Results/Biological_clock_taxonomy_catAB.txt", sep = "\t", quote = FALSE)
  ## Figures: heatmap catA x catB for each bio_class ##
  ## Category abbreviations for readable axis labels
  abbrev<-c(
    "Infectious and parasitic" = "Infect.", "Neoplasms" = "Neopl.", "Blood and blood-forming organs (immune)" = "Blood/Imm.", "Endocrine, nutritional and metabolic" = "Endocr.",
    "Mental and behavioural" = "Mental", "Nervous system" = "Nervous", "Eye and adnexa" = "Eye", "Ear and mastoid process" = "Ear", "Circulatory system" = "Circ.",
    "Respiratory system" = "Resp.", "Digestive system" = "Digest.", "Skin and subcutaneous tissue" = "Skin", "Musculoskeletal system and connective tissue" = "Muscul.",
    "Genitourinary system" = "Genito.", "Pregnancy, childbirth and the puerperium" = "Preg.", "Certain conditions originating in the perinatal period" = "Perinat.",
    "Congenital malformations and chromosomal abnormalities" = "Congenit.", "Symptoms, signs and abnormal laboratory findings" = "Sympt.", "Injury, poisoning" = "Injury"
  )
  ## Long format with abbreviations
  class_long_AB<-melt(class_by_catAB_wide, id.vars = c("catA", "catB"), variable.name = "bio_type", value.name = "pct")
  class_long_AB[, catA_short := abbrev[catA]]
  class_long_AB[, catB_short := abbrev[catB]]
  ## HHeatmap for one bio_type ##
  plot_catAB_heatmap<-function(dt, type_name, title_label, high_color = "#1B4521") {
    x<-dt[bio_type == type_name]
    if(nrow(x) == 0){message("No data for type: ", type_name); return(NULL)}
    ## Order axes by descending mean pct ##
    catA_order<-x[, .(m = mean(pct, na.rm = TRUE)), by = catA_short][order(-m), catA_short]
    catB_order<-x[, .(m = mean(pct, na.rm = TRUE)), by = catB_short][order(-m), catB_short]
    x[, catA_short := factor(catA_short, levels = catA_order)]
    x[, catB_short := factor(catB_short, levels = catB_order)]
    ## Complete grid (show NA cells explicitly) ##
    grid<-CJ(catA_short = levels(x$catA_short), catB_short = levels(x$catB_short))
    x_plot<-merge(grid, x, by = c("catA_short", "catB_short"), all.x = TRUE)
    ggplot(x_plot, aes(x = catB_short, y = catA_short, fill = pct)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = ifelse(!is.na(pct) & pct >= 5, sprintf("%.0f", pct), "")), size = 2.5, color = "grey20") +
      scale_fill_gradient(low = "white", high = high_color, name = "% of pairs", na.value = "grey93") +
      labs(x = "Secondary disease category (B)", y = "Index disease category (A)", title = title_label) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold",hjust=0.5),
        legend.position = "right"
      )
  }
  
  ## Omnipresent ##
  p_omni<-plot_catAB_heatmap(class_long_AB, "Omnipresent (all 9 windows)", "Omnipresent comorbidities by category pair A -> B\n(detectable in all 9 windows)", high_color = "#1B4521")
  if(!is.null(p_omni)){ggsave("ManuscriptFiles/Plots/Biological_clock_taxonomy_catAB_omnipresent.pdf", p_omni, width = 12, height = 10, useDingbats = FALSE)}
  
  ## Long-window only ##
  p_long<-plot_catAB_heatmap(class_long_AB, "Long-window only (late cumulative, no early signal)", "Long-window-only comorbidities by category pair A -> B\n(requires >=3 years follow-up, no signal in year 0-1)", high_color = "#A04000")
  if(!is.null(p_long)){ggsave("ManuscriptFiles/Plots/Biological_clock_taxonomy_catAB_longwindow.pdf",p_long, width = 12, height = 10, useDingbats = FALSE)}
  
  ## Persistent cumulative + partial conditional ##
  p_pers<-plot_catAB_heatmap(class_long_AB, "Persistent cumulative + partial conditional", "Persistent comorbidities with partial conditional risk by category pair A -> B", high_color = "#1B6CA8")
  if(!is.null(p_pers)){ggsave("ManuscriptFiles/Plots/Biological_clock_taxonomy_catAB_persistent_conditional.pdf", p_pers, width = 12, height = 10, useDingbats = FALSE)}
  
  ## Late cumulative + conditional risk ##
  p_late<-plot_catAB_heatmap( class_long_AB, "Late cumulative + conditional risk", "Late-emerging comorbidities with conditional risk by category pair A -> B", high_color = "#E8861A")
  if(!is.null(p_late)){ggsave("ManuscriptFiles/Plots/Biological_clock_taxonomy_catAB_late_conditional.pdf",p_late, width = 12, height = 10, useDingbats = FALSE)}
  
  #### Summary: top catA x catB pairs per bio_class (for paper/supplement) ####
  top_pairs_by_class<-class_by_catAB[!is.na(catA) & !is.na(catB), .SD[order(-pct)][1:5], by = bio_class]
  fwrite(top_pairs_by_class, "ManuscriptFiles/Results/Biological_clock_taxonomy_top_catAB_pairs.txt", sep = "\t", quote = FALSE)
}

#### Compare the results from our network to the one by Westergaard et al. ####
if (args[1] == "compare_westergaard") {
  if(!dir.exists("ManuscriptFiles/Plots/Westergaard")){dir.create("ManuscriptFiles/Plots/Westergaard", recursive = TRUE)}
  if(!dir.exists("ManuscriptFiles/Results/Westergaard")){dir.create("ManuscriptFiles/Results/Westergaard", recursive = TRUE)}
  code<-c(paste0("A0",0:9), paste0("A",10:99),
            paste0("B0",0:9), paste0("B",10:99),
            paste0("C0",0:9), paste0("C",10:99),
            paste0("D0",0:9), paste0("D",10:48), paste0("D",50:89),
            paste0("E0",0:9), paste0("E",10:99),
            paste0("F0",0:9), paste0("F",10:99),
            paste0("G0",0:9), paste0("G",10:99),
            paste0("H0",0:9), paste0("H",10:59), paste0("H",60:95),
            paste0("I0",0:9), paste0("I",10:99),
            paste0("J0",0:9), paste0("J",10:99),
            paste0("K0",0:9), paste0("K",10:93),
            paste0("L0",0:9), paste0("L",10:99),
            paste0("M0",0:9), paste0("M",10:99),
            paste0("N0",0:9), paste0("N",10:99),
            paste0("O0",0:9), paste0("O",10:99),
            paste0("P0",0:9), paste0("P",10:96),
            paste0("Q0",0:9), paste0("Q",10:99),
            paste0("R0",0:9), paste0("R",10:99),
            paste0("S0",0:9), paste0("S",10:99),
            paste0("T0",0:9), paste0("T",10:98))
  cate<-c(rep("I",200), rep("II",149), rep("III",40),
            rep("IV",100), rep("V",100), rep("VI",100),
            rep("VII",60), rep("VIII",36), rep("IX",100),
            rep("X",100), rep("XI",94), rep("XII",100),
            rep("XIII",100), rep("XIV",100), rep("XV",100),
            rep("XVI",97), rep("XVII",100), rep("XVIII",100),
            rep("XIX",199))
  
  catname<-c("Infectious","Neoplasms","Blood/Immune","Endocrine/Metabolic","Mental","Nervous","Eye","Ear","Circulatory","Respiratory",
               "Digestive","Skin","Musculoskeletal","Genitourinary","Pregnancy","Perinatal","Congenital","Symptoms","Injury")
  names(catname)<-unique(cate)
  names(cate)<-code
  ## Overlap metrics helper ##
  overlap_metrics<-function(A, B) {
    inter<-length(intersect(A, B))
    data.table(n_A = length(A), n_B = length(B), n_shared = inter, jaccard = inter/length(union(A, B)), A_in_B = inter/length(A), B_in_A = inter/length(B))
  }
  
  #### Load Westergaard data ####
  soren_raw<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt", stringsAsFactors = FALSE, sep = "\t")
  ## Parse into three populations following original script:
  ## both:  cols 1,2,3,6,7  -> A, B, number_both, rrmin_both, rr_both
  ## women: cols 1,2,4,8,9  -> A, B, number_women, rrmin_women, rr_women
  ## men:   cols 1,2,5,10,11-> A, B, number_men, rrmin_men, rr_men
  parse_soren<-function(raw, pop) {
    if (pop == "both") {
      dt<-raw[, .(disease_a = A, disease_b = B, number = nWomen + nMen, rrmin = adjustedRR2.5, rr = adjustedRR50, dir_low = adjustedDirection2.5, dir_mid = adjustedDirection50, dir_high = adjustedDirection97.5)]
    } else if (pop == "women") {
      dt<-raw[, .(disease_a = A, disease_b = B, number = nWomen, rrmin = womenRR2.5, rr = womenRR50, dir_low = womenDirection2.5, dir_mid = womenDirection50, dir_high = womenDirection97.5)]
    } else {
      dt<-raw[, .(disease_a = A, disease_b = B, number = nMen, rrmin = menRR2.5, rr = menRR50, dir_low = menDirection2.5, dir_mid = menDirection50, dir_high = menDirection97.5)]
    }
    dt<-dt[!is.na(rrmin) & !is.na(number)]
    ## Westergaard filter: rrmin >= 1.01 & number >= 100
    dt<-dt[rrmin >= 1.01 & number >= 100]
    dt[, preferred_dir_den := !is.na(dir_low) & dir_low > 0.5]
    dt[, population := pop]
    dt[, pair     := paste(disease_a, disease_b, sep = "_")]
    dt[, pair_rev := paste(disease_b, disease_a, sep = "_")]
    dt[, catA := catname[cate[disease_a]]]
    dt[, catB := catname[cate[disease_b]]]
    dt
  }
  soren<-rbindlist(lapply(c("both","women","men"), parse_soren, raw = soren_raw))
  
  #### Load Catalonia networks (window 0-5, from prepare_networks) ####
  load_cat_net<-function(pop) {
    f<-sprintf("Results/networks/RR_net_%s_0_5.txt", pop)
    if (!file.exists(f)) {
      cat(sprintf("WARNING: %s not found, trying shrinkage file\n", f))
      f<-sprintf("Results/shrinkage_rr_events/RR_contingency_%s_0_5_shrunk.txt", pop)
    }
    dt<-fread(f)
    dt<-dt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    dt[, population := pop]
    dt[, pair := paste(disease_a, disease_b, sep = "_")]
    dt[, catA := catname[cate[disease_a]]]
    dt[, catB := catname[cate[disease_b]]]
    dt
  }
  
  #### Restrict to common diseases ####
  dis_cat_universe<-readRDS("CaseControlStudy/valid_diseases_any.rds")
  dis_den_universe<-unique(c(soren_raw$A, soren_raw$B))
  dis_common_universe<-intersect(dis_cat_universe, dis_den_universe)
  cat_net_c<-cat_net[disease_a %in% dis_common_universe & disease_b %in% dis_common_universe]
  soren_c<-soren[disease_a %in% dis_common_universe & disease_b %in% dis_common_universe]
  cat(sprintf("After restriction — Cat both: %d, Den both: %d\n", cat_net_c[population == "both", .N], soren_c[  population == "both", .N]))
  
  #### Overlap metrics by population ####
  pops<-c("both","women","men")
  overlap_results<-rbindlist(lapply(pops, function(pop) {
    cat_pairs<-cat_net_c[population == pop, pair]
    den_pairs<-soren_c[  population == pop, pair]
    om<-overlap_metrics(cat_pairs, den_pairs)
    om[, population := pop]
    om
  }))
  fwrite(overlap_results, "ManuscriptFiles/Results/Westergaard/Overlap_metrics_Cat_vs_Den.txt", sep = "\t", quote = FALSE)
  
  #### Test if overlap is greater than expected by chance ####
  ## Universe of possible directed pairs (excluding self-pairs)
  n_universe<-length(dis_common_universe) * (length(dis_common_universe) - 1)
  cat(sprintf("Directed pair universe: %d\n", n_universe))
  overlap_hyper<-rbindlist(lapply(pops, function(pop) {
    cat_pairs<-cat_net_c[population == pop, pair]
    den_pairs<-soren_c[  population == pop, pair]
    n_cat<-length(cat_pairs)
    n_den<-length(den_pairs)
    n_shared<-length(intersect(cat_pairs, den_pairs))
    p_hyper<-phyper(n_shared - 1, n_cat, n_universe - n_cat, n_den, lower.tail = FALSE)
    ## Expected overlap by chance
    expected<-(n_cat / n_universe) * n_den
    fold_enrichment<-n_shared / expected
    cat(sprintf(" %s — observed=%d, expected=%.1f, fold=%.1fx, p=%.2e\n", pop, n_shared, expected, fold_enrichment, p_hyper))
    data.table(population = pop, n_universe = n_universe, n_cat = n_cat, n_den = n_den, n_shared_obs = n_shared, n_shared_exp = round(expected, 1), fold_enrichment = round(fold_enrichment, 2), p_hypergeom = p_hyper)
  }))
  print(overlap_hyper)
  fwrite(overlap_hyper,"ManuscriptFiles/Results/Westergaard/Overlap_hypergeometric_test.txt", sep = "\t", quote = FALSE)
  
  #### UpSet plot - both population  ####
  upset_list<-list(
    "Catalonia (both)" = cat_net_c[population == "both",  pair],
    "Denmark (both)" = soren_c[  population == "both",  pair],
    "Catalonia (women)" = cat_net_c[population == "women", pair],
    "Denmark (women)" = soren_c[  population == "women", pair],
    "Catalonia (men)" = cat_net_c[population == "men",   pair],
    "Denmark (men)" = soren_c[  population == "men",   pair]
  )
  upset_data<-fromList(upset_list)
  pdf("ManuscriptFiles/Plots/Westergaard/UpSet_Catalonia_vs_Denmark.pdf", width = 14, height = 7)
    upset(upset_data,sets = colnames(upset_data),keep.order = TRUE,order.by = "freq",decreasing = TRUE,mainbar.y.label = "N comorbidity pairs",sets.x.label = "Total pairs per set",text.scale = 1.3)
  dev.off()
  
  #### RR concordance for shared pairs ####
  concordance_results<-rbindlist(lapply(pops, function(pop) {
    cat_p<-cat_net_c[population == pop, .(pair, RR_cat = RR_shrunk, catA, catB)]
    den_p<-soren_c[  population == pop, .(pair, RR_den = rr)]
    shared<-merge(cat_p, den_p, by = "pair")
    if(nrow(shared) == 0){return(NULL)}
    cor_all<-cor.test(shared$RR_cat, shared$RR_den, method = "spearman", exact = FALSE)
    cor_log<-cor.test(log(shared$RR_cat), log(shared$RR_den), method = "spearman", exact = FALSE)
    cat(sprintf("%s — n_shared=%d, rho=%.3f (p=%.2e), rho_log=%.3f (p=%.2e)\n", pop, nrow(shared), cor_all$estimate, cor_all$p.value, cor_log$estimate, cor_log$p.value))
    ## Save shared pairs with both RR estimates ##
    fwrite(shared, sprintf("ManuscriptFiles/Results/Westergaard/Shared_pairs_%s.txt", pop), sep = "\t", quote = FALSE)
    data.table(population = pop, n_shared = nrow(shared), rho_rr = cor_all$estimate, p_rho_rr = cor_all$p.value, rho_log_rr = cor_log$estimate, p_rho_log_rr = cor_log$p.value)
  }))
  
  fwrite(concordance_results, "ManuscriptFiles/Results/Westergaard/RR_concordance_summary.txt", sep = "\t", quote = FALSE)
  
  ## Plot concordance: scatter Cat vs Den RR for shared pairs (both) ##
  shared_both<-fread("ManuscriptFiles/Results/Westergaard/Shared_pairs_both.txt")
  shared_both[, catA := catname[cate[gsub("_.+","",pair)]]]
  p_conc<-ggplot(shared_both, aes(x = log(RR_cat), y = log(RR_den), color = catA)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    scale_color_manual(values = setNames(c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C","#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54","#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89"),catname), name = "Category A") +
    labs(x = "log(RR) Catalonia",y = "log(RR) Denmark",
      title = "RR concordance for shared comorbidities (window 0-5)",
      subtitle = sprintf("Spearman rho (log scale) = %.3f, n = %d shared pairs", concordance_results[population == "both", rho_log_rr], concordance_results[population == "both", n_shared])
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(), legend.position  = "right")
  ggsave("ManuscriptFiles/Plots/Westergaard/RR_concordance_both.pdf", p_conc, width = 12, height = 9, useDingbats = FALSE)
  
  #### Enrichment catA x catB — shared vs Catalonia-only ####
  enrich_catAB_overlap_v2<-function(cat_pairs_dt, den_pairs_dt) {
    cat_pairs_dt<-cat_pairs_dt[!is.na(catA) & !is.na(catB)]
    den_pairs_dt<-den_pairs_dt[!is.na(catA) & !is.na(catB)]
    ## All pairs in either dataset ##
    all_pairs<-unique(rbindlist(list(
      cat_pairs_dt[, .(pair, catA, catB, in_cat = TRUE,  in_den = pair %in% den_pairs_dt$pair)],
      den_pairs_dt[!pair %in% cat_pairs_dt$pair, .(pair, catA, catB, in_cat = FALSE, in_den = TRUE)]), fill = TRUE))
    all_pairs[is.na(in_cat), in_cat := FALSE]
    n_shared <-all_pairs[in_cat == TRUE & in_den == TRUE, .N]
    n_cat_only<-all_pairs[in_cat == TRUE & in_den == FALSE, .N]
    n_den_only<-all_pairs[in_cat == FALSE & in_den == TRUE, .N]
    n_neither <-0L
    cat(sprintf("Shared: %d, Cat-only: %d, Den-only: %d\n", n_shared, n_cat_only, n_den_only))
    pairs_AB<-unique(all_pairs[, .(catA, catB)])
    
    res<-rbindlist(lapply(seq_len(nrow(pairs_AB)), function(i) {
      ca<-pairs_AB$catA[i]
      cb<-pairs_AB$catB[i]
      a<-all_pairs[catA == ca & catB == cb & in_cat == TRUE & in_den == TRUE,  .N]
      b<-all_pairs[catA == ca & catB == cb & in_cat == TRUE & in_den == FALSE, .N]
      c<-all_pairs[catA == ca & catB == cb & in_cat == FALSE & in_den == TRUE,  .N]
      b1<-n_shared - a
      b2<-n_cat_only - b
      b3<-n_den_only - c
      ## Fisher: overrepresentation in shared vs rest ##
      in_catAB <-a + b + c
      not_in_catAB<-(n_shared + n_cat_only + n_den_only) - in_catAB
      mat<-matrix(c(a, b + c, n_shared - a, not_in_catAB - (b + c)), nrow = 2, byrow = TRUE)
      if(any(mat < 0) || any(!is.finite(mat))){return(NULL)}
      ft<-tryCatch(fisher.test(mat), error = function(e) NULL)
      if(is.null(ft)){return(NULL)}
      data.table(catA = ca, catB = cb, n_shared = a, n_cat_only = b, n_den_only = c, n_total = a + b + c, odds = unname(ft$estimate), p = ft$p.value)
    }), fill = TRUE)
    
    res<-res[!is.na(p)]
    res[, p_adj := p.adjust(p, method = "BH")]
    res[, log_odds := log(odds)]
    res[, sig := p_adj <= 0.05]
    res
  }
  ## Run for each population ##
  enrich_list_v2<-lapply(pops, function(pop) {
    cat_p<-copy(cat_net_c[population == pop])
    den_p<-copy(soren_c[population == pop])
    den_p[, catA := catname[cate[gsub("_.+","",pair)]]]
    den_p[, catB := catname[cate[gsub(".+_","",pair)]]]
    enrich_catAB_overlap_v2(cat_p, den_p)
  })
  names(enrich_list_v2)<-pops
  for(pop in pops){if(!is.null(enrich_list_v2[[pop]])){fwrite(enrich_list_v2[[pop]],sprintf("ManuscriptFiles/Results/Westergaard/CatAB_enrichment_shared_%s.txt", pop), sep = "\t", quote = FALSE)}}
  
  ## Plot heatmap for both population ##
  plot_catAB_enrich_heatmap<-function(dt, title_label) {
    all_catA<-sort(unique(dt$catA))
    all_catB<-sort(unique(dt$catB))
    grid<-CJ(catA = all_catA, catB = all_catB)
    grid<-merge(grid, dt, by = c("catA","catB"), all.x = TRUE)
    grid[, fill_value := fifelse(is.na(sig) | sig == FALSE, 0, log_odds)]
    ggplot(grid, aes(x = catB, y = catA, fill = fill_value)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(midpoint = 0, low = "#C53030", mid = "white", high = "#2B6CB0", na.value = "grey93", name = "log(OR)\n(sig only)") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold",hjust=0.5)) +
      labs(x = "Secondary disease category (B)", y = "Index disease category (A)", title = title_label)
  }
  p_enr_both_v2<-plot_catAB_enrich_heatmap(enrich_list_v2[["both"]], "Category pairs enriched among shared Catalonia-Denmark comorbidities (both)\n[vs union of both networks]")
  ggsave("ManuscriptFiles/Plots/Westergaard/CatAB_enrichment_shared_both.pdf", p_enr_both_v2, width = 13, height = 11, useDingbats = FALSE)
  
  #### Catalonia-only vs Denmark-only enrichment ####
  enrich_catAB_population<-function(cat_dt, den_dt, pop_label) {
    cat_dt[, source := "Catalonia"]
    den_dt[, source := "Denmark"]
    ## All pairs from both sources
    all_pairs<-rbindlist(list(cat_dt[, .(pair, catA, catB, source)], den_dt[, .(pair, catA, catB, source)]), fill = TRUE)
    all_pairs<-all_pairs[!is.na(catA) & !is.na(catB)]
    n_cat<-all_pairs[source == "Catalonia", .N]
    n_den<-all_pairs[source == "Denmark", .N]
    if(n_cat == 0 || n_den == 0){return(NULL)}
    pairs_AB<-unique(all_pairs[, .(catA, catB)])
    res<-rbindlist(lapply(seq_len(nrow(pairs_AB)), function(i) {
      ca<-pairs_AB$catA[i]
      cb<-pairs_AB$catB[i]
      a1<-all_pairs[source == "Catalonia" & catA == ca & catB == cb, .N]
      a0<-all_pairs[source == "Denmark"   & catA == ca & catB == cb, .N]
      b1<-n_cat - a1
      b0<-n_den - a0
      if (any(c(a1,a0,b1,b0) < 0)) return(NULL)
      ft<-tryCatch(fisher.test(matrix(c(a1,b1,a0,b0), nrow = 2, byrow = TRUE)), error = function(e) NULL)
      if (is.null(ft)) return(NULL)
      data.table(catA = ca, catB = cb, n_cat = a1, n_den = a0, odds = unname(ft$estimate), p = ft$p.value)
    }), fill = TRUE)
    res<-res[!is.na(p)]
    res[, p_adj := p.adjust(p, method = "BH")]
    res[, log_odds := log(odds)]
    res[, sig := p_adj <= 0.05]
    res[, population := pop_label]
    res
  }
  enrich_pop_list<-lapply(pops, function(pop) {enrich_catAB_population(copy(cat_net_c[population == pop]), copy(soren_c[  population == pop]), pop)})
  enrich_pop_all<-rbindlist(enrich_pop_list, fill = TRUE)
  fwrite(enrich_pop_all, "ManuscriptFiles/Results/Westergaard/CatAB_enrichment_Cat_vs_Den.txt", sep = "\t", quote = FALSE)
  
  ## Plot: heatmap for both, positive = more in Catalonia, negative = more in Denmark ##
  plot_catAB_pop_heatmap<-function(dt, pop, title_label) {
    x<-dt[population == pop & !is.na(sig)]
    all_catA<-sort(unique(x$catA))
    all_catB<-sort(unique(x$catB))
    grid<-CJ(catA = all_catA, catB = all_catB)
    grid<-merge(grid, x, by = c("catA","catB"), all.x = TRUE)
    grid[, fill_value := fifelse(is.na(sig) | sig == FALSE, 0, log_odds)]
    ggplot(grid, aes(x = catB, y = catA, fill = fill_value)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(midpoint = 0, low = "#8B1A1A", mid = "white", high = "#1A5276", na.value = "grey93", name = "log(OR)\n(sig only)") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x  = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold",hjust=0.5)) +
      labs(x = "Secondary disease category (B)", y = "Index disease category (A)", title = title_label)
  }
  p_pop_both<-plot_catAB_pop_heatmap(enrich_pop_all, "both", "Category pairs enriched in Catalonia vs Denmark (both, window 0-5)")
  ggsave("ManuscriptFiles/Plots/Westergaard/CatAB_Cat_vs_Den_both.pdf", p_pop_both, width = 13, height = 11, useDingbats = FALSE)
  
  ## Sex comparison panel (women vs men side by side) ##
  p_pop_w<-plot_catAB_pop_heatmap(enrich_pop_all, "women", "Catalonia vs Denmark — Women")
  p_pop_m<-plot_catAB_pop_heatmap(enrich_pop_all, "men", "Catalonia vs Denmark — Men")
  p_pop_sex<-p_pop_w | p_pop_m
  ggsave("ManuscriptFiles/Plots/Westergaard/CatAB_Cat_vs_Den_sex.pdf", p_pop_sex, width = 22, height = 11, useDingbats = FALSE)
  
  #### RR concordance by catA x catB for shared pairs ####
  shared_both_full<-merge(
    cat_net_c[population == "both", .(pair, RR_cat = RR_shrunk, catA, catB)],
    soren_c[  population == "both", .(pair, RR_den = rr)],
    by = "pair"
  )
  ## Spearman rho per catA x catB ##
  rho_by_catAB<-shared_both_full[
    !is.na(catA) & !is.na(catB),
    {
      if (.N < 10) .(rho = NA_real_, p = NA_real_, n = .N)
      else {
        ct<-cor.test(log(RR_cat), log(RR_den), method = "spearman", exact = FALSE)
        .(rho = ct$estimate, p = ct$p.value, n = .N)
      }
    },
    by = .(catA, catB)
  ]
  rho_by_catAB[!is.na(p), p_adj := p.adjust(p, method = "BH")]
  fwrite(rho_by_catAB, "ManuscriptFiles/Results/Westergaard/RR_concordance_by_catAB.txt", sep = "\t", quote = FALSE)
  ## Heatmap of rho
  all_catA<-sort(unique(rho_by_catAB$catA))
  all_catB<-sort(unique(rho_by_catAB$catB))
  rho_grid<-CJ(catA = all_catA, catB = all_catB)
  rho_grid<-merge(rho_grid, rho_by_catAB, by = c("catA","catB"), all.x = TRUE)
  p_rho_heatmap<-ggplot(rho_grid, aes(x = catB, y = catA, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(data = rho_grid[!is.na(rho) & !is.na(p_adj) & p_adj <= 0.05], aes(label = sprintf("%.2f", rho)), size = 2.5, color = "grey20") +
    scale_fill_gradient2(midpoint = 0, low = "#C53030", mid = "white", high = "#1B6CA8", na.value = "grey93", limits = c(-1, 1), name = "Spearman rho\n(log RR)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), plot.title = element_text(face = "bold")) +
    labs(x = "Secondary disease category (B)", y = "Index disease category (A)",
      title = "RR concordance between Catalonia and Denmark by category pair",
      subtitle = "Values shown only for FDR-significant correlations (>=10 shared pairs)"
    )
  ggsave("ManuscriptFiles/Plots/Westergaard/RR_concordance_catAB_heatmap.pdf", p_rho_heatmap, width = 13, height = 11, useDingbats = FALSE)
  
  #### Directional concordance for shared pairs ####
  cat_dir<-fread("Results/networks/RR_net_both_0_5.txt")
  cat_dir<-cat_dir[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
  cat_dir[, pair := paste(disease_a, disease_b, sep = "_")]
  cat_dir[, pair_rev := paste(disease_b, disease_a, sep = "_")]
  ## Shared pairs with direction info from both datasets ##
  shared_dir<-merge( cat_dir[, .(pair, preferred_dir_cat = preferred_direction, theta_cat = theta)], soren_c[population == "both", .(pair, preferred_dir_den, dir_mid_den = dir_mid)], by = "pair")
  ## Concordance ##
  shared_dir[, dir_concordant := preferred_dir_cat == preferred_dir_den]
  cat(sprintf("Shared directed pairs: %d\n  Direction concordant: %d (%.1f%%)\n", nrow(shared_dir), shared_dir[dir_concordant == TRUE, .N], mean(shared_dir$dir_concordant, na.rm = TRUE) * 100))
  ## Directional reversals: Cat prefers A->B, Denmark prefers B->A ##
  soren_both_dir<-soren_c[population == "both"]
  reversed_dir<-merge(cat_dir[preferred_direction == TRUE, .(pair, pair_rev, theta_cat = theta)], soren_both_dir[preferred_dir_den == TRUE, .(pair_rev = pair, theta_den_rev = dir_mid)], by = "pair_rev")
  reversed_dir[, catA := catname[cate[gsub("_.+","",pair)]]]
  reversed_dir[, catB := catname[cate[gsub(".+_","",pair)]]]
  cat(sprintf("Pairs with directional reversal Cat vs Den: %d\n", nrow(reversed_dir)))
  fwrite(shared_dir, "ManuscriptFiles/Results/Westergaard/Directional_concordance_both.txt", sep = "\t", quote = FALSE)
  fwrite(reversed_dir, "ManuscriptFiles/Results/Westergaard/Directional_reversals_Cat_vs_Den.txt", sep = "\t", quote = FALSE)
  
  #### Summary table for paper ####
  summary_westergaard<-rbindlist(lapply(pops, function(pop) {
    cat_p<-cat_net_c[population == pop, pair]
    den_p<-soren_c[  population == pop, pair]
    shared<-intersect(cat_p, den_p)
    om<-overlap_metrics(cat_p, den_p)
    conc<-concordance_results[population == pop]
    data.table(Population = pop, N_Catalonia = length(cat_p), N_Denmark = length(den_p), N_shared = length(shared),
      Jaccard = round(om$jaccard, 4), Pct_Cat_in_Den = round(om$A_in_B * 100, 2), Pct_Den_in_Cat = round(om$B_in_A * 100, 2),
      Spearman_rho_logRR = round(conc$rho_log_rr, 3), p_Spearman = signif(conc$p_rho_log_rr, 3))
  }))
  print(summary_westergaard)
  fwrite(summary_westergaard, "ManuscriptFiles/Results/Westergaard/Summary_table_Cat_vs_Den.txt", sep = "\t", quote = FALSE)
}

if(args[1]=="create_network_for_sharing"){
  if("ComorbidityNetworks"%in%list.files()==FALSE){dir.create("ComorbidityNetworks")}
  ficheros<-list.files("Results/networks/")
  for(a in 1:length(ficheros)){
    # a<-1
    tt<-fread(paste("Results/networks/",ficheros[a],sep=""),stringsAsFactors = F)
    tt<-tt[lfsr < 0.05 & CI_low_RR_shrunk >= 1.01 & cases_event >= 100]
    tt<-tt[,c(1,2,23,27,30,33)]
    write.table(tt,paste("ComorbidityNetworks/",ficheros[a],sep=""),quote=F,sep="\t",row.names=F)
    print(paste(round((a/length(ficheros))*100,2),"%",sep=""))
  }
}