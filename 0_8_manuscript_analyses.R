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
})

## Provide the arguments ##
args<-commandArgs(trailingOnly = TRUE)

## create directories ##
if("Bayes"%in%list.files("ManuscriptFiles/")==FALSE){dir.create("ManuscriptFiles/Bayes")}
if("BayesSummaries"%in%list.files("ManuscriptFiles/")==FALSE){dir.create("ManuscriptFiles/BayesSummaries")}
if("BayesDirectional"%in%list.files("ManuscriptFiles/")==FALSE){dir.create("ManuscriptFiles/BayesDirectional")}
if("BayesNetworks"%in%list.files("ManuscriptFiles/")==FALSE){dir.create("ManuscriptFiles/BayesNetworks")}
if("Plots"%in%list.files("ManuscriptFiles/")==FALSE){dir.create("ManuscriptFiles/Plots")}
if("Results"%in%list.files("ManuscriptFiles/")==FALSE){dir.create("ManuscriptFiles/Results")}
if("Disease_prevalence_comorbidity_number"%in%list.files("ManuscriptFiles/Plots/")==FALSE){dir.create("ManuscriptFiles/Plots/Disease_prevalence_comorbidity_number")}
if("DenmarkOverlap"%in%list.files("ManuscriptFiles/Plots/")==FALSE){dir.create("ManuscriptFiles/Plots/DenmarkOverlap")}
if("DenmarkCorrelation"%in%list.files("ManuscriptFiles/Plots/")==FALSE){dir.create("ManuscriptFiles/Plots/DenmarkCorrelation")}
if("UpsetTimeWindows"%in%list.files("ManuscriptFiles/Plots/")==FALSE){dir.create("ManuscriptFiles/Plots/UpsetTimeWindows")}

## load functions ##
overlapmetrics <- function(A, B) {
  inter <- intersect(A, B)
  union <- union(A, B)
  data.frame(
    jaccard = length(inter) / length(union),
    C_in_D = length(inter) / length(A),
    D_in_C = length(inter) / length(B)
  )
}

## Remove from the tables generated the diseases that do not have enough control samples ##
if(args[1]=="Remove_diseases_without_enough_controls"){
  ## Create the directory where all the files will be saved ##
  if("ManuscriptFiles"%in%list.files()==FALSE){dir.create("ManuscriptFiles")}
  if("IntermediateFiles"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/IntermediateFiles")}
  if("IntermediateFilesByAge"%in%list.files("ManuscriptFiles")==FALSE){dir.create("ManuscriptFiles/IntermediateFilesByAge")}
  ## We have followed two approaches originally, for each disease identify the set of potential controls of interest,
  ## and then select for each patient the first 5 individuals from the control group matching the selected features.
  ## Alternatively, we have done the same but selecting randomly 5 individuals from the control group matching the features.
  ## Since there are some diseases with not enough control individuals we are going to remove them from the final network
  ## To do so, we will remove those for which we don't have controls for more than 1% of de individuals
  ## Remember that when selecting controls randomly we are only including those patients for which we have 5 control samples
  ## Load the number of patients per disease - We select this table as NOW we are focused on comorbidity relationships, so only considering diagnoses ##
  ## done with enough time before the patient leaves the study ##
  enfs<-read.csv2("Data/List_of_diseases_with_diagnoses_5_years_before_leaving_the_study.txt",stringsAsFactors = F,sep="\t")
  enfs<-enfs[order(enfs[,2],decreasing = T),]
  diseases<-enfs$cod
  tablaresumen_random<-c()
  ## Start a loop to anotate, for each disease, the number of individuals with 5 controls assigned vs. the total number of patients
  for(a in 1:length(diseases)){
    ## List the tables containing, per row, the number of controls identified for each patient (it will be always 5)
    files<-list.files(paste("Data/Numbers_first_5years_random/",diseases[a],sep=""))
    if(length(files)==0){tablaresumen_random<-rbind(tablaresumen_random,c(NA,NA,NA,NA))}
    if(length(files)>0){
      vect<-c()
      ## Since we are working with several files per diseases for the different ages and genders, put all of them together ##
      for(b in 1:length(files)){
        tt<-fread(paste("Data/Numbers_first_5years_random/",diseases[a],"/",files[b],sep=""),stringsAsFactors = F,header=F)
        vect<-c(vect,tt$V1)
      }
      tablaresumen_random<-rbind(tablaresumen_random,c(length(vect),mean(vect),max(vect),min(vect)))
    }
    print(paste(round((a/length(diseases))*100,2),"%",sep=""))
  }
  tablaresumen_random<-as.data.table(tablaresumen_random)
  ## Add the name of the disease and the total number of patients with the disease
  tablaresumen_random<-cbind(diseases,enfs$npatients,tablaresumen_random)
  colnames(tablaresumen_random)<-c("diseases","patients","patientwithcontrols","mean","max","min")
  sin_missing_controls<-tablaresumen_random[-which(tablaresumen_random$patients!=tablaresumen_random$patientwithcontrols)]
  write.table(sin_missing_controls,"ManuscriptFiles/Number_patients_per_disease.txt",quote=F,sep="\t",row.names=F)
  ## Diseases for which we don't have enough control samples ##
  missingcontrols<-tablaresumen_random[which(tablaresumen_random$patients!=tablaresumen_random$patientwithcontrols)]
  ## Get percentage of missing patients due to the lack of controls ##
  tutu<-cbind(missingcontrols$diseases,((missingcontrols$patients-missingcontrols$patientwithcontrols)/missingcontrols$patients)*100)
  tutu<-as.data.table(tutu)
  ## Diseases to be removed (those for which we have missed more than 1% of the patients due to not finding enough controls) ##
  removediseases<-tutu[which(tutu$V2>1)]$V1
  if(args[2]=="Global"){
    ## Load the tables generated ##
    ficheros<-list.files("FinalNetworks/AllGlobal/")
    summarytable<-c()
    for(a in 1:length(ficheros)){
      tt<-fread(paste("FinalNetworks/AllGlobal/",ficheros[a],sep=""),stringsAsFactors = F,sep="\t")
      ## get the indexes of the diseases to be removed ##
      rone<-c() ; rtwo<-c()
      for(b in 1:length(removediseases)){
        rone<-c(rone,which(tt$Disease1==removediseases[b]))
        rtwo<-c(rtwo,which(tt$Disease2==removediseases[b]))
      }
      toberemoved<-unique(c(rone,rtwo))
      if(length(toberemoved)>0){tt<-tt[-toberemoved]}
      FDR<-tt$pvalue*length(tt$pvalue)
      if(length(which(FDR>1))>0){FDR[which(FDR>1)]<-1}
      tt$FDR<-FDR
      write.table(tt,paste("ManuscriptFiles/IntermediateFiles/",ficheros[a],sep=""),quote=F,sep="\t",row.names = F)
      ## get the number of nodes and edges in the network ##
      sigtt<-tt[which(tt$FDR<=0.05)]
      summarytable<-rbind(summarytable,c(length(unique(c(tt$Disease1,tt$Disease2))),length(unique(c(sigtt$Disease1,sigtt$Disease2))),length(tt$Disease1),length(sigtt$Disease1),max(table(c(sigtt$Disease1,sigtt$Disease2))),names(which(table(c(sigtt$Disease1,sigtt$Disease2))==max(table(c(sigtt$Disease1,sigtt$Disease2)))))[1]))
      print(paste(round((a/length(ficheros))*100,2),"%",sep=""))
    }
    colnames(summarytable)<-c("diseases","sdiseases","comorbidities","scomorbidities","maxdegree","maxdegreedisease")
    rownames(summarytable)<-gsub(".txt","",ficheros)
    write.table(summarytable,"ManuscriptFiles/Number_edges_and_nodes_per_network.txt",quote=F,sep="\t")
  }
  if(args[2]=="ByAge"){
    ## Load the tables generated ##
    ficheros<-list.files("FinalNetworks/AllByAge/")
    ficheros<-ficheros[-grep("Adjusted",ficheros)]
    tablas<-unique(gsub("Adjusted_.+","Adjusted",gsub("_Men.+","_Men",gsub("_Women.+","_Women",ficheros))))
    for(z in 1:length(tablas)){
      # z<-1
      subfich<-ficheros[grep(tablas[z],ficheros)]
      subfich<-subfich[order(as.numeric(gsub("-.+","",gsub(paste(tablas[z],"_",sep=""),"",subfich))),decreasing=F)]
      summarytable<-c()
      for(a in 1:length(subfich)){
        tt<-fread(paste("FinalNetworks/AllByAge/",subfich[a],sep=""),stringsAsFactors = F,sep="\t")
        ## get the indexes of the diseases to be removed ##
        rone<-c() ; rtwo<-c()
        for(b in 1:length(removediseases)){
          rone<-c(rone,which(tt$Disease1==removediseases[b]))
          rtwo<-c(rtwo,which(tt$Disease2==removediseases[b]))
        }
        toberemoved<-unique(c(rone,rtwo))
        if(length(toberemoved)>0){tt<-tt[-toberemoved]}
        ## Concatenate all the tables indicate the gender and age analysed ##
        if(length(grep("_Men_",subfich[a]))>0){summarytable<-rbind(summarytable,cbind(tt[,c(1,2)],"Men",gsub(".+en_","",gsub(".txt","",subfich[a])),tt[,c(8:11)]))}
        if(length(grep("_Women_",subfich[a]))>0){summarytable<-rbind(summarytable,cbind(tt[,c(1,2)],"Women",gsub(".+en_","",gsub(".txt","",subfich[a])),tt[,c(8:11)]))}
        # print(paste(round((a/length(subfich))*100,2),"%",sep=""))
      }
      colnames(summarytable)[3:4]<-c("Gender","Age")
      ## Create a numeric column for the lower bound of Age
      summarytable[, AgeLower := as.numeric(sub("-.*", "", Age))]
      ## Order by Disease pair, Gender (Women first), and AgeLower
      setorder(summarytable, 
               Disease1, 
               Disease2, 
               AgeLower
      )
      ## remove helper column
      summarytable[, AgeLower := NULL]
      write.table(summarytable,paste("ManuscriptFiles/IntermediateFilesByAge/",tablas[z],".txt",sep=""),quote=F,sep="\t",row.names = F)
      print(paste(round((z/length(tablas))*100,2),"%",sep=""))
    }
    ## Create the adjusted table also ##
    ficheros<-list.files("ManuscriptFiles/IntermediateFilesByAge/")
    if(max(table(gsub("_Adjusted","",gsub("_Men","",gsub("_Women","",gsub(".txt","",ficheros))))))==2){
      fich<-unique(gsub("_Men","",gsub("_Women","",gsub(".txt","",ficheros))))
      for(a in 1:length(fich)){
        # a<-1
        files<-ficheros[grep(fich[a],ficheros)]
        tabla<-c()
        for(b in 1:length(files)){
          tt<-fread(paste("ManuscriptFiles/IntermediateFilesByAge/",files[b],sep=""),stringsAsFactors = F)
          tabla<-rbind(tabla,tt)
        }
        ## Create a numeric column for the lower bound of Age
        tabla[, AgeLower := as.numeric(sub("-.*", "", Age))]
        ## Order by Disease pair, Gender (Women first), and AgeLower
        setorder(tabla, 
                 Disease1, 
                 Disease2, 
                 AgeLower
        )
        ## remove helper column
        tabla[, AgeLower := NULL]
        write.table(tabla,paste("ManuscriptFiles/IntermediateFilesByAge/",fich[a],"_Adjusted.txt",sep=""),quote=F,sep="\t",row.names = F)
        print(paste(round((a/length(fich))*100,2),"%",sep=""))
      }
    }
  }
}

#### Generate comorbidity networks using a bayesian shrinkage approach ####
if(args[1]=="build_networks"){
  #### Apply bayesian shrinkage approach - for disease co-ocurrences ####
  ficheros<-list.files("ManuscriptFiles/IntermediateFiles/")
  ficheros<-setdiff(ficheros,list.files("ManuscriptFiles/Bayes/"))
  for(z in 1:length(ficheros)){
    print(paste("Starting with:",gsub(".txt","",ficheros[z])))
    # z<-13
    tt<-fread(paste("ManuscriptFiles/IntermediateFiles/",ficheros[z],sep=""),stringsAsFactors = F,sep="\t",quote=F)
    pair_counts<-tt[,c(1,2,8:11)]
    colnames(pair_counts)<-c("disease_A", "disease_B", "a", "b", "c", "d")
    setDT(pair_counts)
    setnames(
      pair_counts,
      old = names(pair_counts),
      new = c("A", "B", "a", "b", "c", "d")
    )
    ## Correct for 0s (in the number of patients with or without the disease of interest) ##
    cc <- 0.5
    pair_counts[, `:=`(
      a_cc = a + cc,
      b_cc = b + cc,
      c_cc = c + cc,
      d_cc = d + cc
    )]
    ## Calculate log(OR) and standard error ##
    pair_counts[, `:=`(
      log_or_hat = log((a_cc * d_cc) / (b_cc * c_cc)),
      se_hat = sqrt(1/a_cc + 1/b_cc + 1/c_cc + 1/d_cc)
    )]
    ## Remove infinite values ##
    pair_counts <- pair_counts[
      is.finite(log_or_hat) & is.finite(se_hat)
    ]
    ## Estimate Empirical Bayesian shrinkage (ashr) ##
    ash_fit <- ash(
      pair_counts$log_or_hat,
      pair_counts$se_hat,
      mixcompdist = "normal"
    )
    ## Extract posterior summaries ##
    pair_counts[, `:=`(
      log_or_shrunk = get_pm(ash_fit),
      log_or_se = get_psd(ash_fit),
      lfsr = get_lfsr(ash_fit)
    )]
    ## Back to OR scale ##
    pair_counts[, `:=`(
      OR_shrunk = exp(log_or_shrunk),
      OR_low = exp(log_or_shrunk - 1.96 * log_or_se),
      OR_high = exp(log_or_shrunk + 1.96 * log_or_se)
    )]
    ## Add prevalence of the secondary diseas ##
    pair_counts[, prev_B := (a + c) / (a + b + c + d)]
    ## Select the most relevant comorbidity relationships ##
    strong_pairs <- pair_counts[
      lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100
    ]
    ## Save the resulting network for backbone analysis ##
    strong_pairs<-strong_pairs[,c(1:3,16)]
    colnames(strong_pairs)<-c("Disease1","Disease2","N","OR")
    write.table(strong_pairs,paste("ManuscriptFiles/BayesNetworks/",ficheros[z],sep=""),quote=F,sep="\t",row.names = F)
    ## Quality checks ##
    summaryhat<-summary(pair_counts$log_or_hat)
    summaryshrunk<-summary(pair_counts$log_or_shrunk)
    summaries<-rbind(summaryhat,summaryshrunk)
    
    ## Save the files ##
    write.table(pair_counts,paste("ManuscriptFiles/Bayes/",ficheros[z],sep=""),quote=F,sep="\t",row.names = F)
    write.table(summaries,paste("ManuscriptFiles/BayesSummaries/",ficheros[z],sep=""),quote=F,sep="\t")
    
    #### Now go for significant directionality ####
    ## Emparejar A→B con su inverso B→A ##
    pairs<-pair_counts ; colnames(pairs)[1:2]<-c("disease_A","disease_B")
    dir_pairs <- merge(
      pairs,
      pairs,
      by.x = c("disease_A", "disease_B"),
      by.y = c("disease_B", "disease_A"),
      suffixes = c("_AB", "_BA")
    )
    ## Diferencia de efectos direccionales ##
    dir_pairs[, `:=`(
      delta_log_or = log_or_shrunk_AB - log_or_shrunk_BA,
      delta_se = sqrt(log_or_se_AB^2 + log_or_se_BA^2)
    )]
    ## Shrinkage bayesiano sobre direccionalidad ##
    ash_dir <- ash(
      dir_pairs$delta_log_or,
      dir_pairs$delta_se,
      mixcompdist = "normal"
    )
    dir_pairs[, `:=`(
      delta_shrunk = get_pm(ash_dir),
      lfsr_dir = get_lfsr(ash_dir)
    )]
    ## Direccionalidades robustas ##
    directional_pairs <- dir_pairs[
      lfsr_dir < 0.05 & abs(delta_shrunk) > log(1.3)
    ]
    ## Write directional associations ##
    write.table(dir_pairs,paste("ManuscriptFiles/BayesDirectional/",ficheros[z],sep=""),quote=F,sep="\t",row.names = F)
    #### Now follow Soren Brunaks approach ####
    pairs2<-pairs[,1:3]
    pairs2$a<-as.numeric(pairs2$a)
    pairs2[
      pairs2,
      on = .(disease_A = disease_B, disease_B = disease_A),
      opposite_a := fcoalesce(i.a, 0)
    ]
    if(length(which(is.na(pairs2$opposite_a)))>0){pairs2$opposite_a[which(is.na(pairs2$opposite_a))]<-0}
    ## Calculate significance of directionality ##
    pairs2[, p_value :=
             binom.test(a, a + opposite_a, p = 0.5, alternative = "greater")$p.value,
           by = .(disease_A, disease_B)
    ]
    pairs<-cbind(pair_counts,pairs2$p_value)
    colnames(pairs)[length(colnames(pairs))]<-"directionality_pval"
    write.table(pairs,paste("ManuscriptFiles/Bayes/",ficheros[z],sep=""),quote=F,sep="\t",row.names = F)
    print(paste(round((z/length(ficheros))*100,2),"%",sep=""))
  }
}

#### Analyse gender-related differences on disease prevalence and age of diagnoses, including comparison with Denmark ####
#### Plot gender-associated disease prevalence differences in Catalonia and Denmark and look for differences between them ####
if(args[1]=="compare_prevalences"){
  #### Catalonia ####
  ## @ @ @@ @@ @ @ ##
  dt<-fread("../Entregables/Diagnoses_20080101_20181231_first_diagnoses_of_each_disease.txt",stringsAsFactors = F,sep="|")
  prevalences <- dcast(dt,cod ~ sexe,value.var = "idp",fun.aggregate = uniqueN)
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
  ## disease category ##
  catename<-as.character(catname[cate[code]])
  names(catename)<-code
  
  colcod<-c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C",
            "#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54",
            "#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
  names(colcod)<-unique(cate)
  
  catcol<-as.character(colcod) ; names(catcol)<-catname
  
  codcol<-names(colcod) ; names(codcol)<-as.character(colcod)
  distocol<-as.character(colcod[as.character(cate)]) ; names(distocol)<-names(cate)
  
  prevalences<-cbind(prevalences,as.character(catname[cate[prevalences$disease]]),as.character(distocol[prevalences$disease]))
  colnames(prevalences)[7:8]<-c("diseasecategory","color")
  if(length(which(is.na(prevalences$diseasecategory)))>0){prevalences<-prevalences[-which(is.na(prevalences$diseasecategory))]}
  prevalences<-prevalences[,c(1,7,3:6,8)]
  ## Look for differences in gender prevalence disease by disease ##
  fisher_sex_test <- function(n_w, n_m, N_w, N_m){
    mat <- matrix(c(n_w, N_w - n_w,n_m, N_m - n_m),nrow = 2,byrow = TRUE)
    test <- fisher.test(mat)
    return(c(p_value = test$p.value,OR = unname(test$estimate)))
  }
  res <- t(mapply(fisher_sex_test,prevalences$women,prevalences$men,prevalences$nwomen,prevalences$nmen))
  prevalences$p_value <- res[, "p_value"]
  prevalences$OR      <- res[, "OR"]
  # Correct for multiple testing
  prevalences$p_adj <- p.adjust(prevalences$p_value, method = "BH")
  
  length(intersect(which(prevalences$p_adj<=0.05),which(prevalences$OR>=1.3)))
  length(intersect(which(prevalences$p_adj<=0.05),which(prevalences$OR<=(1/1.3))))
  
  # Define thresholds
  or_threshold <- 1.3
  p_threshold  <- 0.05
  
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
  prevalences_plot <- prevalences2[!is.na(OR) & !is.na(total_cases)]
  
  ## Stablish the order of the disease categories ##
  disease_order <- unique(prevalences_plot$diseasecategory)
  prevalences_plot[, diseasecategory := factor(diseasecategory, levels = disease_order)]
  
  ## Plot the association
  pdf("ManuscriptFiles/Plots/Gender_bias_by_disease.pdf",width = 10,height = 10)
  ggplot(
    prevalences_plot,
    aes(x = total_cases,y = OR,color = color)) +
    geom_point(size = 3,alpha = 0.6) +
    scale_color_identity(
      name = "disease category",
      guide = "legend",
      breaks = prevalences_plot$color[
        match(levels(prevalences_plot$diseasecategory),
              prevalences_plot$diseasecategory)
      ],
      labels = levels(prevalences_plot$diseasecategory)
    ) +
    xlab("Log (total cases, women + men)") +
    ylab("Log (OR)") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  dev.off()
  
  #### Make the plot interactive ####
  dig3dis<-fread("Data/ICD10_three_digits_names.txt",stringsAsFactors = F)
  tdig3dis<-dig3dis$Name ; names(tdig3dis)<-dig3dis$Code
  prevalences2<-cbind(prevalences2,as.character(tdig3dis[prevalences2$disease]))
  colnames(prevalences2)[13]<-c("diseasename")
  
  prevalences_plot <- prevalences2[!is.na(OR) & !is.na(total_cases)]
  # Fijar el orden de las categorías en la leyenda
  disease_order <- unique(prevalences_plot$diseasecategory)
  prevalences_plot[, diseasecategory := factor(diseasecategory, levels = disease_order)]
  # Crear ggplot con tooltip
  p <- ggplot(
    prevalences_plot,
    aes(
      x = total_cases,
      y = OR,
      color = color,
      text = paste0(
        "Disease name: ", diseasename, "<br>",
        "Disease: ", disease, "<br>",
        "Category: ", diseasecategory, "<br>",
        "OR: ", round(exp(OR), 0), "<br>",
        "Women: ", women,"<br>",
        "Men: ", men,"<br>"
      )
    )
  ) +
    geom_point(size = 3, alpha = 0.6) +
    scale_color_identity(
      name = "disease category",
      guide = "legend",
      breaks = prevalences_plot$color[
        match(levels(prevalences_plot$diseasecategory),
              prevalences_plot$diseasecategory)
      ],
      labels = levels(prevalences_plot$diseasecategory)
    ) +
    xlab("Total cases (women + men)") +
    ylab("Odds Ratio (OR)") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  # Convertir a interactivo
  ggplotly(p, tooltip = "text")
  
  #### Enrichment analysis ####
  ## Number of analyzed diseases ##
  n_all <- length(prevalences[, unique(disease)])
  
  enrichment_by_category <- function(dt, category_col, bias_value) {
    # Diseases biased
    biased <- dt[sex_bias == bias_value, unique(disease)]
    # If non-biased diseases --> EXIT
    if (length(biased) == 0) return(NULL)
    
    results <- dt[, {
      diseases_in_cat <- unique(disease)
      a <- sum(diseases_in_cat %in% biased) # biased inside the category
      b <- length(biased) - a # biased outside the category
      c <- length(diseases_in_cat) - a # non-biased inside the category
      d <- n_all - a - b - c # non-biased outside the category
      mat <- matrix(c(a, b, c, d), nrow = 2,dimnames = list(Bias = c("InCategory", "OutCategory"),Status = c("Biased", "NotBiased")))
      ft <- fisher.test(mat)
      list(n_biased_in_category = a,n_diseases_in_category = length(diseases_in_cat),OR_enrichment = unname(ft$estimate),p_value = ft$p.value)
    }, by = category_col]
    
    results[, sex_bias := bias_value]
    return(results)
  }
  # Enrichment for each gender separately
  enrichment_results <- rbindlist(list(enrichment_by_category(prevalences,category_col = "diseasecategory",bias_value = "women"), enrichment_by_category(prevalences,category_col = "diseasecategory",bias_value = "men")), use.names = TRUE)
  # Adjust for multiple testing by gender
  enrichment_results[, p_adj := p.adjust(p_value, method = "BH"),by = sex_bias]
  
  ## Plot gender differences
  plot_data <- enrichment_results[
    , .(
      diseasecategory,
      sex_bias,
      OR_enrichment,
      log_OR = log(OR_enrichment),
      n_biased_in_category,
      p_adj
    )
  ]
  ## Add color palette ##
  plot_data<-cbind(plot_data,as.character(catcol[plot_data$diseasecategory]))
  colnames(plot_data)[7]<-"color"
  ## Find significant enrichments ##
  plot_data[, significant := p_adj <= 0.05]
  plot_data$significant[which(plot_data$significant==TRUE)]<-"Yes"
  plot_data$significant[which(plot_data$significant==FALSE)]<-"No"
  plot_data[, diseasecategory := factor(diseasecategory,levels = rev(unique(diseasecategory)))]
  
  categories <- plot_data %>%
    distinct(diseasecategory, .keep_all = TRUE) %>%
    mutate(disease_label = paste0("<span style='color:", color, "'>\u25CF</span> ", diseasecategory))
  
  categories <- categories %>%
    mutate(disease_label = factor(as.character(disease_label), levels = disease_label))
  
  
  pdf("ManuscriptFiles/Plots/Enrichment_gender_bias_period_prevalence_by_disease_categories.pdf",width = 10,height = 10)
    ggplot(plot_data,
           aes(x = log_OR,
               y = diseasecategory,
               color = sex_bias,
               size = n_biased_in_category,
               alpha = significant)) +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_alpha_manual(values = c("Yes" = 1, "No" = 0.3)) +
      scale_color_manual(
        values = c("women" = "#A94C40", "men" = "#5474A5"),
        labels = c("Men-biased", "Women-biased")
      ) +
      scale_size_continuous(name = "Number of diseases") +
      labs(
        x = "Log odds ratio of enrichment",
        y = "Disease category (ICD-10)", ## NULL en el otro
        color = "Gender bias",
        title = "Enrichment of disease categories among diseases\nwith significant gender-biased diagnosis rates"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", hjust=0.5)
      )
  dev.off()
  
  pdf("ManuscriptFiles/Plots/LabelColors.pdf",width = 10,height = 12)
    plot(rep(1,19),1:19,col=plot_data$color[19:1],pch=16,cex=2)
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
  ## disease category ##
  catename<-as.character(catname[cate[code]])
  names(catename)<-code
  
  colcod<-c("#DD3232","#FAC6CC","#FFB6AD","#FBAF5F","#FACF63","#FFEF6C","#CED75C",
            "#D2EFDB","#A1CE5E","#1A86A8","#00619C","#0065A9","#002B54",
            "#985396","#805462","#E3DFD6","#B2B1A5","#58574B","#C5AB89")
  names(colcod)<-unique(cate)
  
  catcol<-as.character(colcod) ; names(catcol)<-catname
  
  codcol<-names(colcod) ; names(codcol)<-as.character(colcod)
  distocol<-as.character(colcod[as.character(cate)]) ; names(distocol)<-names(cate)
  
  denprevalences<-cbind(denprevalences,as.character(catname[cate[denprevalences$disease]]),as.character(distocol[denprevalences$disease]))
  colnames(denprevalences)[7:8]<-c("diseasecategory","color")
  if(length(which(is.na(denprevalences$diseasecategory)))>0){denprevalences<-denprevalences[-which(is.na(denprevalences$diseasecategory))]}
  denprevalences<-denprevalences[,c(1,7,3:6,8)]
  ## Look for differences in gender prevalence disease by disease ##
  fisher_sex_test <- function(n_w, n_m, N_w, N_m){
    mat <- matrix(c(n_w, N_w - n_w,n_m, N_m - n_m),nrow = 2,byrow = TRUE)
    test <- fisher.test(mat)
    return(c(p_value = test$p.value,OR = unname(test$estimate)))
  }
  res <- t(mapply(fisher_sex_test,denprevalences$women,denprevalences$men,denprevalences$nwomen,denprevalences$nmen))
  denprevalences$p_value <- res[, "p_value"]
  denprevalences$OR      <- res[, "OR"]
  # Correct for multiple testing
  denprevalences$p_adj <- p.adjust(denprevalences$p_value, method = "BH")
  
  length(intersect(which(denprevalences$p_adj<=0.05),which(denprevalences$OR>=1.3)))
  length(intersect(which(denprevalences$p_adj<=0.05),which(denprevalences$OR<=(1/1.3))))
  
  # Define thresholds
  or_threshold <- 1.3
  p_threshold  <- 0.05
  
  denprevalences[, sex_bias := fifelse(
    p_adj <= p_threshold & OR >= or_threshold, "women",
    fifelse(p_adj <= p_threshold & OR <= 1/or_threshold, "men", "none")
  )]
  
  #### Gender bias by disease ####
  denprevalences2<-denprevalences[-which(denprevalences$sex_bias=="none")]
  denprevalences2[, total_cases := women + men]
  denprevalences2$OR<-log(denprevalences2$OR)
  denprevalences2$total_cases<-log(denprevalences2$total_cases)
  # denprevalences2$OR[which(denprevalences2$OR==Inf)]<-11
  # denprevalences2$OR[which(denprevalences2$OR==-Inf)]<--11
  
  ## Check if OR and total cases are valid ## 
  denprevalences_plot <- denprevalences2[!is.na(OR) & !is.na(total_cases)]
  
  ## Stablish the order of the disease categories ##
  disease_order <- unique(denprevalences_plot$diseasecategory)
  denprevalences_plot[, diseasecategory := factor(diseasecategory, levels = disease_order)]
  
  ## Plot the association
  pdf("ManuscriptFiles/Plots/Gender_bias_by_disease_Denmark.pdf",width = 10,height = 10)
  ggplot(
    denprevalences_plot,
    aes(x = total_cases,y = OR,color = color)) +
    geom_point(size = 3,alpha = 0.6) +
    scale_color_identity(
      name = "disease category",
      guide = "legend",
      breaks = denprevalences_plot$color[
        match(levels(denprevalences_plot$diseasecategory),
              denprevalences_plot$diseasecategory)
      ],
      labels = levels(denprevalences_plot$diseasecategory)
    ) +
    xlab("Log (total cases, women + men)") +
    ylab("Log (OR)") +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  dev.off()
  
  #### Enrichment analysis ####
  ## Number of analyzed diseases ##
  n_all <- length(denprevalences[, unique(disease)])
  
  enrichment_by_category <- function(dt, category_col, bias_value) {
    # Diseases biased
    biased <- dt[sex_bias == bias_value, unique(disease)]
    # If non-biased diseases --> EXIT
    if (length(biased) == 0) return(NULL)
    
    results <- dt[, {
      diseases_in_cat <- unique(disease)
      a <- sum(diseases_in_cat %in% biased) # biased inside the category
      b <- length(biased) - a # biased outside the category
      c <- length(diseases_in_cat) - a # non-biased inside the category
      d <- n_all - a - b - c # non-biased outside the category
      mat <- matrix(c(a, b, c, d), nrow = 2,dimnames = list(Bias = c("InCategory", "OutCategory"),Status = c("Biased", "NotBiased")))
      ft <- fisher.test(mat)
      list(n_biased_in_category = a,n_diseases_in_category = length(diseases_in_cat),OR_enrichment = unname(ft$estimate),p_value = ft$p.value)
    }, by = category_col]
    
    results[, sex_bias := bias_value]
    return(results)
  }
  # Enrichment for each gender separately
  denenrichment_results <- rbindlist(list(enrichment_by_category(denprevalences,category_col = "diseasecategory",bias_value = "women"), enrichment_by_category(denprevalences,category_col = "diseasecategory",bias_value = "men")), use.names = TRUE)
  # Adjust for multiple testing by gender
  denenrichment_results[, p_adj := p.adjust(p_value, method = "BH"),by = sex_bias]
  
  ## Plot gender differences
  plot_data_den <- denenrichment_results[
    , .(
      diseasecategory,
      sex_bias,
      OR_enrichment,
      log_OR = log(OR_enrichment),
      n_biased_in_category,
      p_adj
    )
  ]
  ## Add color palette ##
  plot_data_den<-cbind(plot_data_den,as.character(catcol[plot_data_den$diseasecategory]))
  colnames(plot_data_den)[7]<-"color"
  ## Find significant enrichments ##
  plot_data_den[, significant := p_adj <= 0.05]
  plot_data_den$significant[which(plot_data_den$significant==TRUE)]<-"Yes"
  plot_data_den$significant[which(plot_data_den$significant==FALSE)]<-"No"
  plot_data_den[, diseasecategory := factor(diseasecategory,levels = rev(unique(diseasecategory)))]
  
  categories <- plot_data_den %>%
    distinct(diseasecategory, .keep_all = TRUE) %>%
    mutate(disease_label = paste0("<span style='color:", color, "'>\u25CF</span> ", diseasecategory))
  
  categories <- categories %>%
    mutate(disease_label = factor(as.character(disease_label), levels = disease_label))
  
  
  pdf("ManuscriptFiles/Plots/Enrichment_gender_bias_period_prevalence_by_disease_categories_Denmark.pdf",width = 10,height = 10)
  ggplot(plot_data_den,
         aes(x = log_OR,
             y = diseasecategory,
             color = sex_bias,
             size = n_biased_in_category,
             alpha = significant)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_alpha_manual(values = c("Yes" = 1, "No" = 0.3)) +
    scale_color_manual(
      values = c("women" = "#A94C40", "men" = "#5474A5"),
      labels = c("Men-biased", "Women-biased")
    ) +
    scale_size_continuous(name = "Number of diseases") +
    labs(
      x = "Log odds ratio of enrichment",
      y = "Disease category (ICD-10)", ## NULL en el otro
      color = "Gender bias",
      title = "Enrichment of disease categories among diseases\nwith significant gender-biased diagnosis rates"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust=0.5)
    )
  dev.off()
  
  pdf("ManuscriptFiles/Plots/LabelColors_Denmark.pdf",width = 10,height = 12)
    plot(rep(1,19),1:19,col=plot_data_den$color[19:1],pch=16,cex=2)
  dev.off()
  
  #### Compare with Denmark ####
  commondis<-intersect(prevalences$disease,denprevalences$disease)
  denwom<-intersect(commondis,denprevalences$disease[which(denprevalences$sex_bias=="women")])
  catwom<-intersect(commondis,prevalences$disease[which(prevalences$sex_bias=="women")])
  denmen<-intersect(commondis,denprevalences$disease[which(denprevalences$sex_bias=="men")])
  catmen<-intersect(commondis,prevalences$disease[which(prevalences$sex_bias=="men")])
  dennone<-intersect(commondis,denprevalences$disease[which(denprevalences$sex_bias=="none")])
  catnone<-intersect(commondis,prevalences$disease[which(prevalences$sex_bias=="none")])
  biasprevalences<-list("Denmark - Women preference"=denwom,"Denmark - Men preference"=denmen,
                        "Catalonia - Women preference"=catwom,"Catalonia - Men preference"=catmen)
  
  biasprevalences2<-list("Denmark - Women preference"=denwom,"Denmark - Men preference"=denmen,
                        "Catalonia - Women preference"=catwom,"Catalonia - Men preference"=catmen,
                        "Denmark - No preference"=dennone,"Catalonia - No preference"=catnone)
  
  ## Upset plot ##
  ## Gender associated ##
  upset_data <- fromList(biasprevalences)
  pdf("ManuscriptFiles/Plots/Disease_prevalence_biases_in_Catalonia_and_Denmark.pdf",width = 10,height = 6)
    upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = "")
  dev.off()
  venn.diagram(biasprevalences,filename = "ManuscriptFiles/Plots/Disease_prevalence_gender_biases_in_Catalonia_and_Denmark.png",category.names = c("Denmark\nWomen", "Denmark\nMen","Catalonia\nWomen","Catalonia\nMen"),cex = 0.8, main = "Diseases with a gender preference",disable.logging = TRUE)
  
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
  dt<-fread("Diagnoses_20080101_20181231_first_diagnoses_of_each_disease.txt",stringsAsFactors = F,sep="|")
  ## Load the function ##
  ttest_cod <- function(dat) {
    ## Check that both genders are present ##
    if(length(unique(dat$sexe)) != 2) {return(list(diff_mean = NA_real_, p_value = NA_real_))}
    dat_W <- dat[sexe == "D", edad]
    dat_M <- dat[sexe == "H", edad]
    ## At least two individuals by gender
    if(length(dat_W) < 2 | length(dat_M) < 2) {return(list(diff_mean = NA_real_, p_value = NA_real_))}
    ## Test t Welch
    t_test <- t.test(dat_M, dat_W, var.equal = FALSE)
    list(diff_mean = t_test$estimate[1] - t_test$estimate[2],p_value   = t_test$p.value)
  }
  results <- dt[, ttest_cod(.SD), by = cod]
  ## Add relevant columns ##
  summary_cod <- dt[, .(nMen    = sum(sexe == "H"),nWomen    = sum(sexe == "D"),
    MenMean = mean(edad[sexe == "H"], na.rm = TRUE),WomenMean = mean(edad[sexe == "D"], na.rm = TRUE),
    MenSD   = sd(edad[sexe == "H"], na.rm = TRUE),WomenSD   = sd(edad[sexe == "D"], na.rm = TRUE)
  ), by = cod][nWomen >= 2 & nMen >= 2]
  ## Put both tables together and correct for multiple testing ##
  catalonia_age <- merge(summary_cod, results, by = "cod", all.x = TRUE)
  catalonia_age[, fdr := p.adjust(p_value, method = "fdr")]
  colnames(catalonia_age)[c(1,8:10)]<-c("Code","DiffMean","pval","FDR")
  catalonia_age<-cbind(catalonia_age[,1:8],catalonia_age[,10])
  ## Save the table ##
  write.table(catalonia_age,"ManuscriptFiles/Results/SupplementaryData_gender_related_ageofdiagnosis_differences_Catalonia.txt",quote=F,sep="\t",row.names=F)
  
  ## Which diseases are diagnosed later in women? ##
  latterwomen<-catalonia_age[intersect(which(catalonia_age$FDR<=0.05),which(catalonia_age$DiffMean<0))] ; latterwomen<-latterwomen[order(abs(latterwomen$DiffMean),decreasing=T)]
  lattermen<-catalonia_age[intersect(which(catalonia_age$FDR<=0.05),which(catalonia_age$DiffMean>0))] ; lattermen<-lattermen[order(abs(lattermen$DiffMean),decreasing=T)]
  
  #### Boxplot age ####
  ## Add disease category ##
  dt <- merge(catalonia_age,prevalences[, .(disease, diseasecategory)],by.x = "Code",by.y = "disease",all.x = TRUE)
  ## Add disease category color ##
  dt[, color := distocol[Code]]
  dt[, diseasecategory := factor(diseasecategory)]
  if(length(which(is.na(dt$color)))>0){dt<-dt[-which(is.na(dt$color))]}
  
  ## All the points ##
  pdf(file="ManuscriptFiles/Plots/Gender_differences_agediagnosis_all.pdf",width = 12,height = 8)
  ggplot(dt, aes(x = DiffMean, y = diseasecategory, fill = diseasecategory)) +
    geom_violin(alpha = 0.4, color = NA) +                 # densidad
    geom_boxplot(width = 0.1, alpha = 0.6, outlier.shape = NA) +  # mediana y rango
    geom_jitter(aes(color = diseasecategory), 
                size = 2, alpha = 0.5, width = 0.2) +      # puntos individuales
    scale_fill_manual(values = unique(dt$color[order(dt$diseasecategory)])) +
    scale_color_manual(values = unique(dt$color[order(dt$diseasecategory)])) +
    labs(
      x = "DiffMean (Women − Men)",
      y = "Disease category",
      fill = "Category",
      color = "Category"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    )
  dev.off()
  
  ## Only FDR<=0.05 ##
  dt_sig <- dt[FDR <= 0.05]
  pdf(file="ManuscriptFiles/Plots/Gender_differences_agediagnosis_significant.pdf",width = 12,height = 8)
  ggplot(dt_sig, aes(x = DiffMean, y = diseasecategory, fill = diseasecategory)) +
    geom_violin(alpha = 0.4, color = NA) +                 # densidad
    geom_boxplot(width = 0.1, alpha = 0.6, outlier.shape = NA) +  # mediana y rango
    geom_jitter(aes(color = diseasecategory), 
                size = 2, alpha = 0.5, width = 0.2) +      # puntos individuales
    scale_fill_manual(values = unique(dt_sig$color[order(dt_sig$diseasecategory)])) +
    scale_color_manual(values = unique(dt_sig$color[order(dt_sig$diseasecategory)])) +
    labs(
      x = "DiffMean (Women − Men)",
      y = "Disease category",
      fill = "Category",
      color = "Category"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    )
  dev.off()
  
  #### Age vs. Prevalence (sum of the number of patients) ####
  dt <- merge(catalonia_age,prevalences[, .(disease, diseasecategory)],by.x = "Code",by.y = "disease",all.x = TRUE)
  ## Total prevalence ## 
  dt[, prevalence := nMen + nWomen]
  ## Add the color ##
  dt[, color := distocol[Code]]
  ## Identify significances ##
  dt[, signif := ifelse(FDR <= 0.05, "FDR ≤ 0.05", "Not significant")]
  dt[, signif := factor(signif, levels = c("Not significant", "FDR ≤ 0.05"))]
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
  dt <- merge(catalonia_age,prevalences,by.x = "Code",by.y = "disease",all = FALSE)
  ## Log(OR) ##
  dt[, logOR := log(OR)]
  ## Add disease category name ##
  dt[, catename := diseasecategory]
  ## Add the color of each disease ##
  dt[, color := distocol[Code]]
  ## Remove NAs ##
  xlims <- range(dt$DiffMean, na.rm = TRUE)
  ylims <- range(dt$logOR, na.rm = TRUE)
  ## Create a color vector ##
  cat_colors <- unique(dt[, .(catename = diseasecategory, color)])
  cat_colors_vec <- setNames(cat_colors$color, cat_colors$catename)
  ## Plot ##
  p <- ggplot(dt, aes(x = DiffMean, y = logOR)) +
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
  
  p <- ggplot(dt, aes(x = DiffMean, y = logOR)) +
    geom_point(
      aes(
        color = catename,
        text = paste0(
          "Code: ", Code, "<br>",
          "Disease: ", disease_name, "<br>",
          "Category: ", catename, "<br>",
          "DiffMean: ", round(DiffMean, 2), "<br>",
          "log(OR): ", round(logOR, 2)
        )
      ),
      size = 3.5,
      alpha = 0.65
    ) +
    scale_color_manual(values = cat_colors_vec) +
    labs(
      x = "DiffMean (Women − Men)",
      y = "log(OR) (Women vs. Men)",
      color = "Disease category"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90")
    )
  
  ggplotly(p, tooltip = "text")
  
  #### Compare with Denmark ####
  denmark_age<-fread("Epidemiology/41467_2019_8475_MOESM5_ESM.txt",stringsAsFactors = F) ; setkey(denmark_age,"Code")
  setkey(catalonia_age,"Code")
  commondis<-intersect(denmark_age$Code,catalonia_age$Code)
  cdenmark_age<-denmark_age[commondis] ; ccatalonia_age<-catalonia_age[commondis]
  ## Ensure both tables present the same ordering ##
  stopifnot(all(cdenmark_age$Code == ccatalonia_age$Code))
  ## Combine both tables ##
  dt <- data.table(Code = cdenmark_age$Code,Diff_Denmark   = cdenmark_age$DiffMean,Diff_Catalonia = ccatalonia_age$DiffMean,FDR_Denmark    = cdenmark_age$FDR,FDR_Catalonia  = ccatalonia_age$FDR)
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
  lista<-list("Later Women (Catalonia)"=ccatalonia_age$Code[intersect(which(ccatalonia_age$FDR<=0.05),which(ccatalonia_age$DiffMean<0))],
              "Later Men (Catalonia)"=ccatalonia_age$Code[intersect(which(ccatalonia_age$FDR<=0.05),which(ccatalonia_age$DiffMean>0))],
              "Later Women (Denmark)"=cdenmark_age$Code[intersect(which(cdenmark_age$FDR<=0.05),which(cdenmark_age$DiffMean<0))],
              "Later Men (Denmark)"=cdenmark_age$Code[intersect(which(cdenmark_age$FDR<=0.05),which(cdenmark_age$DiffMean>0))])
  upset_data <- fromList(lista)
  pdf("ManuscriptFiles/Plots/Disease_age_gender_biases_in_Catalonia_and_Denmark.pdf",width = 10,height = 6)
    upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = "")
  dev.off()
  venn.diagram(lista,filename = "ManuscriptFiles/Plots/Disease_age_gender_biases_in_Catalonia_and_Denmark.png",category.names = c("Catalonia\nLatter in Women","Catalonia\nLatter in Men","Denmark\nLatter in Women", "Denmark\nLatter in Men"),cex = 0.8, main = "Diseases with different age of diagnosis",disable.logging = TRUE)
  
}


#### Analyse networks ####
#### Calculate correlations between disease prevalence and number of comorbidities (separately for source and sink nodes) ####
if(args[1]=="calculate_correlations_prevalence_number_comorbidities"){
  #### Correlations between number of comorbidities and disease prevalence ####
  diseasesanalyzed<-fread("ManuscriptFiles/Number_patients_per_disease.txt",stringsAsFactors = F)
  prevalences<-fread("Data/DiseaseIncidencesGeneral.txt",stringsAsFactors = F,sep="\t")
  colnames(prevalences)[1]<-"diseases" ; setkey(prevalences,"diseases")
  prevalences<-prevalences[intersect(diseasesanalyzed$diseases,prevalences$diseases)]
  ficheros<-list.files("ManuscriptFiles/Bayes/")
  comparisons<-c("Adjusted","Women","Men")
  ## Incremental_time_windows ##
  ficheros1<-ficheros[grep("N0_",ficheros)]
  for(z in comparisons){
    # z<-"Adjusted"
    fichs<-ficheros1[grep(paste("_",z,sep=""),ficheros1)]
    lesult<-list()
    for(a in fichs){
      # a<-"N0_1_Adjusted.txt"
      ## Load the three tables ##
      bayesall<-fread(paste("ManuscriptFiles/Bayes/",a,sep=""),stringsAsFactors = F,sep="\t")
      ## Get the significant associations ##
      sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
      #### Correlations ####
      if(gsub(".txt","",gsub(".+_","",a))=="Adjusted"){patsbydisease<-prevalences[,1:2] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub(".txt","",gsub(".+_","",a))=="Women"){patsbydisease<-prevalences[,c(1,3)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub(".txt","",gsub(".+_","",a))=="Men"){patsbydisease<-prevalences[,c(1,4)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      ## significant bayes all ##
      ## @@ @@ @@ @ @ @@ @@ @@ ##
      ## Counts in Disease1 ##
      count_d1 <- sbayesall %>% dplyr::count(A, name = "n_Disease1") %>% rename(diseases = A)
      ## Counts in Disease1 ##
      count_d2 <- sbayesall %>% dplyr::count(B, name = "n_Disease2") %>% rename(diseases = B)
      ## Add to patsbydisease
      resulta <- patsbydisease %>% left_join(count_d1, by = "diseases") %>% left_join(count_d2, by = "diseases") %>%
        mutate(n_Disease1 = coalesce(n_Disease1, 0L),n_Disease2 = coalesce(n_Disease2, 0L))
      if(length(intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)))>0){resulta<-resulta[-intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)),]}
      lesult[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]]<-resulta
    }
    ## Plot ##
    pdf(file=paste("ManuscriptFiles/Plots/Disease_prevalence_comorbidity_number/",z,"_incremental_timewindows.pdf",sep=""),width = 18,height = 8)
      par(mfrow = c(2, 5))  # 2 filas, 3 columnas
      ## 1 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 2 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 3 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 4 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 5 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 2 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_2`$patientwithcontrols,lesult$`0_2`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 3 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_3`$patientwithcontrols,lesult$`0_3`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 4 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_4`$patientwithcontrols,lesult$`0_4`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 10 ##
      plot(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 5 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_5`$patientwithcontrols,lesult$`0_5`$n_Disease2)$p.value,6),sep=""),cex=0.7)
    dev.off()
    print(paste(z,"finished!"))
  }
  ## Continuous_time_windows ##
  ficheros2<-c(ficheros[grep("N0_1",ficheros)],ficheros[grep("N1_2",ficheros)],ficheros[grep("N2_3",ficheros)],ficheros[grep("N3_4",ficheros)],ficheros[grep("N4_5",ficheros)])
  for(z in comparisons){
    # z<-"Adjusted"
    fichs<-ficheros2[grep(paste("_",z,sep=""),ficheros2)]
    lesult<-list()
    for(a in fichs){
      # a<-"N0_1_Adjusted.txt"
      ## Load the three tables ##
      bayesall<-fread(paste("ManuscriptFiles/Bayes/",a,sep=""),stringsAsFactors = F,sep="\t")
      ## Get the significant associations ##
      sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
      #### Correlations ####
      if(gsub(".txt","",gsub(".+_","",a))=="Adjusted"){patsbydisease<-prevalences[,1:2] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub(".txt","",gsub(".+_","",a))=="Women"){patsbydisease<-prevalences[,c(1,3)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(gsub(".txt","",gsub(".+_","",a))=="Men"){patsbydisease<-prevalences[,c(1,4)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      ## significant bayes all ##
      ## @@ @@ @@ @ @ @@ @@ @@ ##
      ## Counts in Disease1 ##
      count_d1 <- sbayesall %>% dplyr::count(A, name = "n_Disease1") %>% rename(diseases = A)
      ## Counts in Disease1 ##
      count_d2 <- sbayesall %>% dplyr::count(B, name = "n_Disease2") %>% rename(diseases = B)
      ## Add to patsbydisease
      resulta <- patsbydisease %>% left_join(count_d1, by = "diseases") %>% left_join(count_d2, by = "diseases") %>%
        mutate(n_Disease1 = coalesce(n_Disease1, 0L),n_Disease2 = coalesce(n_Disease2, 0L))
      if(length(intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)))>0){resulta<-resulta[-intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)),]}
      lesult[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]]<-resulta
    }
    ## Plot ##
    pdf(file=paste("ManuscriptFiles/Plots/Disease_prevalence_comorbidity_number/",z,"_consecutive_timewindows.pdf",sep=""),width = 18,height = 8)
      par(mfrow = c(2, 5))  # 2 filas, 3 columnas
      ## 1 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 1 - 2 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 2 - 3 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 3 - 4 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 4 - 5 years\nsource disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 0 - 1 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`0_1`$patientwithcontrols,lesult$`0_1`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 1 - 2 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`1_2`$patientwithcontrols,lesult$`1_2`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 2 - 3 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`2_3`$patientwithcontrols,lesult$`2_3`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 3 - 4 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`3_4`$patientwithcontrols,lesult$`3_4`$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## 10 ##
      plot(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Time window: 4 - 5 years\nsink disease",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$`4_5`$patientwithcontrols,lesult$`4_5`$n_Disease2)$p.value,6),sep=""),cex=0.7)
    dev.off()
    print(paste(z,"finished!"))
  }
}
#### Analyze correlations and overlaps with Soren Brunak's work, study shared comorbidities between genders, create upset plots for overlaps by timewindow ####
if(args[1]=="all_time_windows_supplementary"){
  if(args[2]=="correlations_and_denmark_overlaps"){
    #### Overlaps with Soren Brunak's networks ####
    soren<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
    soren<-soren[,c(1,2,5,6,7,8,13,14,19,20)]
    nwom<-soren[,3] ; nwom[which(is.na(nwom))]<-0
    nmen<-soren[,4] ; nmen[which(is.na(nmen))]<-0
    soren<-cbind(soren[,1:2],nwom+nmen,soren[,3:10])
    colnames(soren)<-c("A","B","nadj","nwom","nmen","rradj_min","rradj","rrwom_min","rrwom","rrmen_min","rrmen")
    comparisons<-c("Adjusted","Women","Men")
    ## Incremental_time_windows ##
    ficheros<-list.files("ManuscriptFiles/Bayes/")
    ficheros1<-ficheros[grep("N0_",ficheros)]
    dir_overlap_significance<-c()
    direction_information<-c()
    for(z in comparisons){
      # z<-"Adjusted"
      fichs<-ficheros1[grep(paste("_",z,sep=""),ficheros1)]
      lesult<-list()
      vennfiles<-list()
      vennfiles_nodir<-list()
      for(a in fichs){
        # a<-"N0_1_Adjusted.txt"
        ## Load the three tables ##
        bayesall<-fread(paste("ManuscriptFiles/Bayes/",a,sep=""),stringsAsFactors = F,sep="\t")
        ## Select the epidemiology of interest ##
        if(gsub(".txt","",gsub(".+_","",a))=="Adjusted"){epidemiology<-soren[,c(1:3,6,7)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
        if(gsub(".txt","",gsub(".+_","",a))=="Women"){epidemiology<-soren[,c(1:2,4,8,9)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
        if(gsub(".txt","",gsub(".+_","",a))=="Men"){epidemiology<-soren[,c(1:2,5,10,11)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
        ## Remove NAs and select the number of nodes ##
        if(length(which(is.na(epidemiology$rrmin)))>0){epidemiology<-epidemiology[-which(is.na(epidemiology$rrmin))]}
        
        ## Get the nodes ##
        bayesallnodes<-unique(c(bayesall$A,bayesall$B))
        epidemiologynodes<-unique(c(epidemiology$A,epidemiology$B))
        
        ## Common diseases analyzed ##
        bayesallcommon<-intersect(bayesallnodes,epidemiologynodes)
        
        ## Select comorbidities between shared nodes in Catalonia ##
        bayesall <- bayesall[bayesall$A %in% bayesallcommon & bayesall$B %in% bayesallcommon,]
        
        ## Select comorbidities between shared nodes in Catalonia ##
        epidemiology <- epidemiology[epidemiology$A %in% bayesallcommon & epidemiology$B %in% bayesallcommon,]
        
        ## Get the significant associations ##
        sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
        sbayesall<-sbayesall[,c(1,2,3,16)] ; colnames(sbayesall)<-c("A","B","number","rr")
        sepidemiology<-epidemiology[rrmin >=1.01 & number >=100,c(1:3,5)]
        
        ## Plot Venn Diagram of shared interactions ##
        sba<-paste(sbayesall$A,sbayesall$B,sep="_")
        sep<-paste(sepidemiology$A,sepidemiology$B,sep="_")
        vennfiles[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]] <- venn.diagram(x = list(sba, sep),category.names = c("C", "D"),filename = NULL,fill = c("#19787F", "#440F53"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",gsub(paste("_",z,".+",sep=""),"",a))),paste("JI =",round(overlapmetrics(sba, sep)[1],4)),paste("C in D =",round(overlapmetrics(sba, sep)[2]*100,2)),paste("D in C =",round(overlapmetrics(sba, sep)[3]*100,2)),sep="\n"),disable.logging = TRUE)
        
        ## Calculate significance of the overlap ##
        all_pairs <- as.vector(outer(bayesallcommon, bayesallcommon, paste, sep = "_"))
        all_pairs <- all_pairs[!grepl("^(.+)_\\1$", all_pairs)]
        N <- length(all_pairs) ; K <- length(sba) ; n <- length(sep) ; k <- length(intersect(sba, sep))
        pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
        expected <- (K * n) / N ; enrichment <- k / expected
        dir_overlap_significance<-rbind(dir_overlap_significance,c(gsub(".txt","",a),k,expected,N,enrichment,pval))
        
        ## Get the number of patients with diseases and the number of comorbidities ##
        setkey(sbayesall, A, B)
        setkey(sepidemiology, A, B)
        lesult[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]] <- sbayesall[sepidemiology, nomatch = 0]
        
        ## Plot Venn Diagram removing directionality ##
        sbas<-t(apply(sbayesall[,1:2], 1, sort))[-which(duplicated(t(apply(sbayesall[,1:2], 1, sort)))),]
        ssba<-paste(sbas[,1],sbas[,2],sep="_")
        seps<-t(apply(sepidemiology[,1:2], 1, sort))[-which(duplicated(t(apply(sepidemiology[,1:2], 1, sort)))),]
        ssep<-paste(seps[,1],seps[,2],sep="_")
        vennfiles_nodir[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]] <- venn.diagram(x = list(ssba, ssep),category.names = c("C", "D"),filename = NULL,fill = c("#19787F", "#440F53"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",gsub(paste("_",z,".+",sep=""),"",a))),paste("JI =",round(overlapmetrics(ssba, ssep)[1],4)),paste("C in D =",round(overlapmetrics(ssba, ssep)[2]*100,2)),paste("D in C =",round(overlapmetrics(ssba, ssep)[3]*100,2)),sep="\n"),disable.logging = TRUE)
        
        #### Detailed information on directionality ####
        Cat_pairs     <- sba
        Den_pairs     <- sep
        
        ## Comunes ##
        comunes<-intersect(Cat_pairs,Den_pairs)
        
        ## En ambas direcciones ##
        comundt<-cbind(gsub("_.+","",comunes),gsub(".+_","",comunes))
        scomundt<-t(apply(comundt, 1, sort))
        bidirectional<-length(which(duplicated(scomundt)))*2
        
        Cat_only<-setdiff(Cat_pairs,comunes)
        Den_only<-setdiff(Den_pairs,comunes)
        
        Cat_only_rev<-paste(gsub(".+_","",Cat_only),gsub("_.+","",Cat_only),sep="_")
        Den_only_rev<-paste(gsub(".+_","",Den_only),gsub("_.+","",Den_only),sep="_")
        
        CatonlyAB_Comrev<-intersect(Cat_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
        DenonlyAB_Comrev<-intersect(Den_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
        
        CatonlyAB_DenBArev<-intersect(Cat_only,Den_only_rev)
        DenonlyAB_CatBArev<-intersect(Den_only,Cat_only_rev)
        ## Save numbers ##
        direction_information<-rbind(direction_information,c(a,length(Cat_only),length(Den_only),length(comunes),bidirectional,length(CatonlyAB_DenBArev),length(DenonlyAB_CatBArev),length(CatonlyAB_Comrev),length(DenonlyAB_Comrev),
                                                             length(setdiff(ssba,intersect(ssba,ssep))),length(setdiff(ssep,intersect(ssba,ssep))),length(intersect(ssba,ssep))))
      }
      ## Plot the Venn Diagram ##
      pdf(paste("ManuscriptFiles/Plots/DenmarkOverlap/",z,"_incremental_timewindows.pdf",sep=""), width = 20, height = 4)
        grid.arrange(grobs = list(grobTree(vennfiles$`0_1`),grobTree(vennfiles$`0_2`),grobTree(vennfiles$`0_3`),
                                  grobTree(vennfiles$`0_4`),grobTree(vennfiles$`0_5`)),ncol = 5)
      dev.off()
      ## Plot the Venn Diagram no-direction ##
      pdf(paste("ManuscriptFiles/Plots/DenmarkOverlap/NoDirection_",z,"_incremental_timewindows.pdf",sep=""), width = 20, height = 4)
        grid.arrange(grobs = list(grobTree(vennfiles_nodir$`0_1`),grobTree(vennfiles_nodir$`0_2`),grobTree(vennfiles_nodir$`0_3`),
                                  grobTree(vennfiles_nodir$`0_4`),grobTree(vennfiles_nodir$`0_5`)),ncol = 5)
      dev.off()
      
      ## Plot the correlations##
      pdf(file=paste("ManuscriptFiles/Plots/DenmarkCorrelation/",z,"_incremental_timewindows.pdf",sep=""),width = 18,height = 8)
      par(mfrow = c(2, 5))
      ## 1 ##
      plot(lesult$`0_1`$number,lesult$`0_1`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 1 years\nComorbidity prevalence")
      text(max(lesult$`0_1`$number)/3*2,max(lesult$`0_1`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_1`$number,lesult$`0_1`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_1`$number,lesult$`0_1`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$`0_2`$number,lesult$`0_2`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 2 years\nComorbidity prevalence")
      text(max(lesult$`0_2`$number)/3*2,max(lesult$`0_2`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_2`$number,lesult$`0_2`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_2`$number,lesult$`0_2`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$`0_3`$number,lesult$`0_3`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 3 years\nComorbidity prevalence")
      text(max(lesult$`0_3`$number)/3*2,max(lesult$`0_3`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_3`$number,lesult$`0_3`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_3`$number,lesult$`0_3`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$`0_4`$number,lesult$`0_4`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 4 years\nComorbidity prevalence")
      text(max(lesult$`0_4`$number)/3*2,max(lesult$`0_4`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_4`$number,lesult$`0_4`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_4`$number,lesult$`0_4`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$`0_5`$number,lesult$`0_5`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 5 years\nComorbidity prevalence")
      text(max(lesult$`0_5`$number)/3*2,max(lesult$`0_5`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_5`$number,lesult$`0_5`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_5`$number,lesult$`0_5`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$`0_1`$rr,lesult$`0_1`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 1 years\nComorbidity risk")
      text(max(lesult$`0_1`$rr)/3*2,max(lesult$`0_1`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_1`$rr,lesult$`0_1`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_1`$rr,lesult$`0_1`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(lesult$`0_2`$rr,lesult$`0_2`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 2 years\nComorbidity risk")
      text(max(lesult$`0_2`$rr)/3*2,max(lesult$`0_2`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_2`$rr,lesult$`0_2`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_2`$rr,lesult$`0_2`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(lesult$`0_3`$rr,lesult$`0_3`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 3 years\nComorbidity risk")
      text(max(lesult$`0_3`$rr)/3*2,max(lesult$`0_3`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_3`$rr,lesult$`0_3`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_3`$rr,lesult$`0_3`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(lesult$`0_4`$rr,lesult$`0_4`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 4 years\nComorbidity risk")
      text(max(lesult$`0_4`$rr)/3*2,max(lesult$`0_4`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_4`$rr,lesult$`0_4`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_4`$rr,lesult$`0_4`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 10 ##
      plot(lesult$`0_5`$rr,lesult$`0_5`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 5 years\nComorbidity risk")
      text(max(lesult$`0_5`$rr)/3*2,max(lesult$`0_5`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_5`$rr,lesult$`0_5`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_5`$rr,lesult$`0_5`$i.rr)$p.value,6),sep=""),cex=0.7)
      dev.off()
      print(paste(z,"finished!"))
    }
    colnames(direction_information)<-c("network","only_catalonia","only_denmark","commons","bidirectional","only_catalonia_denmark_reverse","only_denmark_catalonia_reverse","only_catalonia_common_reverse","only_denmark_common_reverse",
                                       "only_catalonia_nodirection","only_denmark_nodirection","common_nodirection")
    colnames(dir_overlap_significance)<-c("network","common","expected","all_posible","enrichment","pval_overlap")
    write.table(direction_information,"ManuscriptFiles/Results/Incremental_timewindow_directionality_overlap.txt",quote=F,sep="\t",row.names=F)
    write.table(dir_overlap_significance,"ManuscriptFiles/Results/Incremental_timewindow_catalonia_denmark_overlap_significance.txt",quote=F,sep="\t",row.names=F)
    ## Continuous_time_windows ##
    ficheros2<-c(ficheros[grep("N0_1",ficheros)],ficheros[grep("N1_2",ficheros)],ficheros[grep("N2_3",ficheros)],ficheros[grep("N3_4",ficheros)],ficheros[grep("N4_5",ficheros)])
    dir_overlap_significance<-c()
    direction_information<-c()
    for(z in comparisons){
      # z<-"Adjusted"
      fichs<-ficheros2[grep(paste("_",z,sep=""),ficheros2)]
      lesult<-list()
      vennfiles<-list()
      vennfiles_nodir<-list()
      for(a in fichs){
        # a<-"N0_1_Adjusted.txt"
        ## Load the three tables ##
        bayesall<-fread(paste("ManuscriptFiles/Bayes/",a,sep=""),stringsAsFactors = F,sep="\t")
        ## Select the epidemiology of interest ##
        if(gsub(".txt","",gsub(".+_","",a))=="Adjusted"){epidemiology<-soren[,c(1:3,6,7)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
        if(gsub(".txt","",gsub(".+_","",a))=="Women"){epidemiology<-soren[,c(1:2,4,8,9)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
        if(gsub(".txt","",gsub(".+_","",a))=="Men"){epidemiology<-soren[,c(1:2,5,10,11)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
        ## Remove NAs and select the number of nodes ##
        if(length(which(is.na(epidemiology$rrmin)))>0){epidemiology<-epidemiology[-which(is.na(epidemiology$rrmin))]}
        
        ## Get the nodes ##
        bayesallnodes<-unique(c(bayesall$A,bayesall$B))
        epidemiologynodes<-unique(c(epidemiology$A,epidemiology$B))
        
        ## Common diseases analyzed ##
        bayesallcommon<-intersect(bayesallnodes,epidemiologynodes)
        
        ## Select comorbidities between shared nodes in Catalonia ##
        bayesall <- bayesall[bayesall$A %in% bayesallcommon & bayesall$B %in% bayesallcommon,]
        
        ## Select comorbidities between shared nodes in Catalonia ##
        epidemiology <- epidemiology[epidemiology$A %in% bayesallcommon & epidemiology$B %in% bayesallcommon,]
        
        ## Get the significant associations ##
        sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
        sbayesall<-sbayesall[,c(1,2,3,16)] ; colnames(sbayesall)<-c("A","B","number","rr")
        sepidemiology<-epidemiology[rrmin >=1.01 & number >=100,c(1:3,5)]
        
        ## Plot Venn Diagram of shared interactions ##
        sba<-paste(sbayesall$A,sbayesall$B,sep="_")
        sep<-paste(sepidemiology$A,sepidemiology$B,sep="_")
        vennfiles[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]] <- venn.diagram(x = list(sba, sep),category.names = c("C", "D"),filename = NULL,fill = c("#19787F", "#440F53"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",gsub(paste("_",z,".+",sep=""),"",a))),paste("JI =",round(overlapmetrics(sba, sep)[1],4)),paste("C in D =",round(overlapmetrics(sba, sep)[2]*100,2)),paste("D in C =",round(overlapmetrics(sba, sep)[3]*100,2)),sep="\n"),disable.logging = TRUE)
        
        ## Calculate significance of the overlap ##
        all_pairs <- as.vector(outer(bayesallcommon, bayesallcommon, paste, sep = "_"))
        all_pairs <- all_pairs[!grepl("^(.+)_\\1$", all_pairs)]
        N <- length(all_pairs) ; K <- length(sba) ; n <- length(sep) ; k <- length(intersect(sba, sep))
        pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
        expected <- (K * n) / N ; enrichment <- k / expected
        dir_overlap_significance<-rbind(dir_overlap_significance,c(gsub(".txt","",a),k,expected,N,enrichment,pval))
        
        ## Get the number of patients with diseases and the number of comorbidities ##
        setkey(sbayesall, A, B)
        setkey(sepidemiology, A, B)
        lesult[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]] <- sbayesall[sepidemiology, nomatch = 0]
        
        ## Plot Venn Diagram removing directionality ##
        sbas<-t(apply(sbayesall[,1:2], 1, sort))[-which(duplicated(t(apply(sbayesall[,1:2], 1, sort)))),]
        ssba<-paste(sbas[,1],sbas[,2],sep="_")
        seps<-t(apply(sepidemiology[,1:2], 1, sort))[-which(duplicated(t(apply(sepidemiology[,1:2], 1, sort)))),]
        ssep<-paste(seps[,1],seps[,2],sep="_")
        vennfiles_nodir[[gsub("N","",gsub("_.txt","",gsub(z,"",a)))]] <- venn.diagram(x = list(ssba, ssep),category.names = c("C", "D"),filename = NULL,fill = c("#19787F", "#440F53"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",gsub(paste("_",z,".+",sep=""),"",a))),paste("JI =",round(overlapmetrics(ssba, ssep)[1],4)),paste("C in D =",round(overlapmetrics(ssba, ssep)[2]*100,2)),paste("D in C =",round(overlapmetrics(ssba, ssep)[3]*100,2)),sep="\n"),disable.logging = TRUE)
        
        #### Detailed information on directionality ####
        Cat_pairs     <- sba
        Den_pairs     <- sep
        
        ## Comunes ##
        comunes<-intersect(Cat_pairs,Den_pairs)
        
        ## En ambas direcciones ##
        comundt<-cbind(gsub("_.+","",comunes),gsub(".+_","",comunes))
        scomundt<-t(apply(comundt, 1, sort))
        bidirectional<-length(which(duplicated(scomundt)))*2
        
        Cat_only<-setdiff(Cat_pairs,comunes)
        Den_only<-setdiff(Den_pairs,comunes)
        
        Cat_only_rev<-paste(gsub(".+_","",Cat_only),gsub("_.+","",Cat_only),sep="_")
        Den_only_rev<-paste(gsub(".+_","",Den_only),gsub("_.+","",Den_only),sep="_")
        
        CatonlyAB_Comrev<-intersect(Cat_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
        DenonlyAB_Comrev<-intersect(Den_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
        
        CatonlyAB_DenBArev<-intersect(Cat_only,Den_only_rev)
        DenonlyAB_CatBArev<-intersect(Den_only,Cat_only_rev)
        ## Save numbers ##
        direction_information<-rbind(direction_information,c(a,length(Cat_only),length(Den_only),length(comunes),bidirectional,length(CatonlyAB_DenBArev),length(DenonlyAB_CatBArev),length(CatonlyAB_Comrev),length(DenonlyAB_Comrev),
                                                             length(setdiff(ssba,intersect(ssba,ssep))),length(setdiff(ssep,intersect(ssba,ssep))),length(intersect(ssba,ssep))))
      }
      ## Plot the Venn Diagram ##
      pdf(paste("ManuscriptFiles/Plots/DenmarkOverlap/",z,"_continuous_timewindows.pdf",sep=""), width = 20, height = 4)
        grid.arrange(grobs = list(grobTree(vennfiles$`0_1`),grobTree(vennfiles$`1_2`),grobTree(vennfiles$`2_3`),
                                  grobTree(vennfiles$`3_4`),grobTree(vennfiles$`4_5`)),ncol = 5)
      dev.off()
      
      ## Plot the Venn Diagram no-direction ##
      pdf(paste("ManuscriptFiles/Plots/DenmarkOverlap/NoDirection_",z,"_continuous_timewindows.pdf",sep=""), width = 20, height = 4)
        grid.arrange(grobs = list(grobTree(vennfiles_nodir$`0_1`),grobTree(vennfiles_nodir$`1_2`),grobTree(vennfiles_nodir$`2_3`),
                                  grobTree(vennfiles_nodir$`3_4`),grobTree(vennfiles_nodir$`4_5`)),ncol = 5)
      dev.off()
      
      ## Plot the correlations##
      pdf(file=paste("ManuscriptFiles/Plots/DenmarkCorrelation/",z,"_continuous_timewindows.pdf",sep=""),width = 18,height = 8)
      par(mfrow = c(2, 5))
      ## 1 ##
      plot(lesult$`0_1`$number,lesult$`0_1`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 1 years\nComorbidity prevalence")
      text(max(lesult$`0_1`$number)/3*2,max(lesult$`0_1`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_1`$number,lesult$`0_1`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_1`$number,lesult$`0_1`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$`1_2`$number,lesult$`1_2`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 1 - 2 years\nComorbidity prevalence")
      text(max(lesult$`1_2`$number)/3*2,max(lesult$`1_2`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`1_2`$number,lesult$`1_2`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`1_2`$number,lesult$`1_2`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$`2_3`$number,lesult$`2_3`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 2 - 3 years\nComorbidity prevalence")
      text(max(lesult$`2_3`$number)/3*2,max(lesult$`2_3`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`2_3`$number,lesult$`2_3`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`2_3`$number,lesult$`2_3`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$`3_4`$number,lesult$`3_4`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 3 - 4 years\nComorbidity prevalence")
      text(max(lesult$`3_4`$number)/3*2,max(lesult$`3_4`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`3_4`$number,lesult$`3_4`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`3_4`$number,lesult$`3_4`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$`4_5`$number,lesult$`4_5`$i.number,xlab="Catalonia",ylab="Denmark",main="Time window: 4 - 5 years\nComorbidity prevalence")
      text(max(lesult$`4_5`$number)/3*2,max(lesult$`4_5`$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$`4_5`$number,lesult$`4_5`$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`4_5`$number,lesult$`4_5`$i.number)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$`0_1`$rr,lesult$`0_1`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 0 - 1 years\nComorbidity risk")
      text(max(lesult$`0_1`$rr)/3*2,max(lesult$`0_1`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`0_1`$rr,lesult$`0_1`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`0_1`$rr,lesult$`0_1`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(lesult$`1_2`$rr,lesult$`1_2`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 1 - 2 years\nComorbidity risk")
      text(max(lesult$`1_2`$rr)/3*2,max(lesult$`1_2`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`1_2`$rr,lesult$`1_2`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`1_2`$rr,lesult$`1_2`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(lesult$`2_3`$rr,lesult$`2_3`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 2 - 3 years\nComorbidity risk")
      text(max(lesult$`2_3`$rr)/3*2,max(lesult$`2_3`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`2_3`$rr,lesult$`2_3`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`2_3`$rr,lesult$`2_3`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(lesult$`3_4`$rr,lesult$`3_4`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 3 - 4 years\nComorbidity risk")
      text(max(lesult$`3_4`$rr)/3*2,max(lesult$`3_4`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`3_4`$rr,lesult$`3_4`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`3_4`$rr,lesult$`3_4`$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 10 ##
      plot(lesult$`4_5`$rr,lesult$`4_5`$i.rr,xlab="Catalonia",ylab="Denmark",main="Time window: 4 - 5 years\nComorbidity risk")
      text(max(lesult$`4_5`$rr)/3*2,max(lesult$`4_5`$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$`4_5`$rr,lesult$`4_5`$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$`4_5`$rr,lesult$`4_5`$i.rr)$p.value,6),sep=""),cex=0.7)
      dev.off()
      print(paste(z,"finished!"))
    }
    colnames(direction_information)<-c("network","only_catalonia","only_denmark","commons","bidirectional","only_catalonia_denmark_reverse","only_denmark_catalonia_reverse","only_catalonia_common_reverse","only_denmark_common_reverse",
                                       "only_catalonia_nodirection","only_denmark_nodirection","common_nodirection")
    colnames(dir_overlap_significance)<-c("network","common","expected","all_posible","enrichment","pval_overlap")
    write.table(direction_information,"ManuscriptFiles/Results/Continuos_timewindow_directionality_overlap.txt",quote=F,sep="\t",row.names=F)
    write.table(dir_overlap_significance,"ManuscriptFiles/Results/Continuos_timewindow_catalonia_denmark_overlap_significance.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="shared_comorbidities_by_sex_catalonia_varying_timewindows_venn_diagrams"){
    ficheros<-list.files("ManuscriptFiles/Bayes/")
    #### Incremental_time_windows ####
    ficheros1<-ficheros[grep("N0_",ficheros)]
    rangos<-unique(gsub("_Men.+","",gsub("_Wom.+","",gsub("_Adj.+","",ficheros1))))
    comparisons<-c("Women","Adjusted","Men")
    direction_information<-c() ; ldirection_information<-list() ; lsdirection_information<-list()
    vennfiles<-list() ; vennfiles_directional<-list()
    lupsetincremental<-list() ; lupsetincrementalsigdir<-list()
    for(z in rangos){
      # z<-"N0_1"
      fichs<-sort(ficheros1[grep(paste(z,"_",sep=""),ficheros1)])
      lnetp<-list() ; lnet<-list() ; lnetdp<-list() ; lnetd<-list()
      for(a in fichs){
        # a<-"N0_5_Women.txt"
        ## Load the three tables ##
        bayesall<-fread(paste("ManuscriptFiles/Bayes/",a,sep=""),stringsAsFactors = F,sep="\t")
        
        ## Get the significant associations ##
        sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
        dsbayesall<-sbayesall[which(p.adjust(sbayesall$directionality_pval, method = "BH")<=0.05)]
        sbayesall<-sbayesall[,c(1,2,3,16)] ; colnames(sbayesall)<-c("A","B","number","rr")
        dsbayesall<-dsbayesall[,c(1,2,3,16)] ; colnames(dsbayesall)<-c("A","B","number","rr")
        lnet[[gsub(".txt","",gsub(".+_","",a))]]<-sbayesall
        lnetd[[gsub(".txt","",gsub(".+_","",a))]]<-dsbayesall
        lnetp[[gsub(".txt","",gsub(".+_","",a))]]<-paste(sbayesall$A,sbayesall$B,sep="_")
        lnetdp[[gsub(".txt","",gsub(".+_","",a))]]<-paste(dsbayesall$A,dsbayesall$B,sep="_")
        lupsetincremental[[gsub(".txt","",gsub(".+_","",a))]][[gsub("_","-",gsub("N","",z))]]<-paste(sbayesall$A,sbayesall$B,sep="_")
        lupsetincrementalsigdir[[gsub(".txt","",gsub(".+_","",a))]][[gsub("_","-",gsub("N","",z))]]<-paste(dsbayesall$A,dsbayesall$B,sep="_")
      }
      vennfiles[[gsub("_","-",gsub("N","",z))]] <- venn.diagram(x = lnetp,category.names = substr(names(lnetp),1,1),filename = NULL,fill = c("#19787F", "#4A75A9","#B6453A"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
      ## Plot Venn Diagram with directionality ##
      vennfiles_directional[[gsub("_","-",gsub("N","",z))]] <- venn.diagram(x = lnetdp,category.names = substr(names(lnetdp),1,1),filename = NULL,fill = c("#19787F", "#4A75A9","#B6453A"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
      
      #### Detailed information on directionality ####
      Wom_pairs     <- lnetp$Women
      Men_pairs     <- lnetp$Men
      ## Comunes ##
      comunes<-intersect(Wom_pairs,Men_pairs)
      ## En ambas direcciones ##
      comundt<-cbind(gsub("_.+","",comunes),gsub(".+_","",comunes))
      scomundt<-t(apply(comundt, 1, sort))
      bidirectional<-length(which(duplicated(scomundt)))*2
      Wom_only<-setdiff(Wom_pairs,comunes)
      Men_only<-setdiff(Men_pairs,comunes)
      Wom_only_rev<-paste(gsub(".+_","",Wom_only),gsub("_.+","",Wom_only),sep="_")
      Men_only_rev<-paste(gsub(".+_","",Men_only),gsub("_.+","",Men_only),sep="_")
      WomonlyAB_Comrev<-intersect(Wom_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
      MenonlyAB_Comrev<-intersect(Men_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
      WomonlyAB_MenBArev<-intersect(Wom_only,Men_only_rev)
      MenonlyAB_WomBArev<-intersect(Men_only,Wom_only_rev)
      ## Save numbers ##
      direction_information<-rbind(direction_information,c(z,length(Wom_only),length(Men_only),length(comunes),bidirectional,length(WomonlyAB_MenBArev),length(MenonlyAB_WomBArev),length(WomonlyAB_Comrev),length(MenonlyAB_Comrev),
                                                           length(setdiff(ssba,intersect(ssba,ssep))),length(setdiff(ssep,intersect(ssba,ssep))),length(intersect(ssba,ssep))))
      ldirection_information[[z]]$WomonlyAB_MenBArev<-WomonlyAB_MenBArev
      ldirection_information[[z]]$MenonlyAB_WomBArev<-MenonlyAB_WomBArev
      
      #### Detailed information on significant directionality ####
      Wom_pairs     <- lnetdp$Women
      Men_pairs     <- lnetdp$Men
      ## Comunes ##
      comunes<-intersect(Wom_pairs,Men_pairs)
      Wom_only<-setdiff(Wom_pairs,comunes)
      Men_only<-setdiff(Men_pairs,comunes)
      Wom_only_rev<-paste(gsub(".+_","",Wom_only),gsub("_.+","",Wom_only),sep="_")
      Men_only_rev<-paste(gsub(".+_","",Men_only),gsub("_.+","",Men_only),sep="_")
      WomonlyAB_MenBArev<-intersect(Wom_only,Men_only_rev)
      MenonlyAB_WomBArev<-intersect(Men_only,Wom_only_rev)
      ## Save numbers ##
      lsdirection_information[[z]]$WomonlyAB_MenBArev<-WomonlyAB_MenBArev
      lsdirection_information[[z]]$MenonlyAB_WomBArev<-MenonlyAB_WomBArev
      print(paste(z,"finished!"))
    }
    colnames(direction_information)<-c("network","only_women","only_men","commons","bidirectional","only_women_men_reverse","only_men_women_reverse","only_women_common_reverse","only_men_common_reverse",
                                       "only_women_nodirection","only_men_nodirection","common_nodirection")
    ## Save agreements in table and list ##
    write.table(direction_information,"ManuscriptFiles/Results/Women_Men_overlap_comorbidities_incremental_timewindow.txt",quote=F,sep="\t",row.names=F)
    save(ldirection_information,lsdirection_information,file="ManuscriptFiles/IntermediateFiles/List_comorbidities_gender_agreement_incremental.RData")
    ## Plot the Venn Diagram ##
    pdf("ManuscriptFiles/Plots/Sex_agreement_comorbidities_incremental_timewindows.pdf", width = 20, height = 4)
    grid.arrange(grobs = list(grobTree(vennfiles$`0-1`),grobTree(vennfiles$`0-2`),grobTree(vennfiles$`0-3`),
                              grobTree(vennfiles$`0-4`),grobTree(vennfiles$`0-5`)),ncol = 5)
    dev.off()
    ## Plot the Venn Diagram no-direction ##
    pdf(paste("ManuscriptFiles/Plots/Sex_agreement_significant_directional_comorbidities_incremental_timewindows.pdf",sep=""), width = 20, height = 4)
      grid.arrange(grobs = list(grobTree(vennfiles_directional$`0-1`),grobTree(vennfiles_directional$`0-2`),grobTree(vennfiles_directional$`0-3`),
                                grobTree(vennfiles_directional$`0-4`),grobTree(vennfiles_directional$`0-5`)),ncol = 5)
    dev.off()
    ## Save properly gender agreement diagram for Catalonia
    lista<-lnetp
    lnetp<-list("Women"=lnetp$Women,"Men"=lnetp$Men,"Adjusted"=lnetp$Adjusted)
    venn.diagram(x = lnetp,category.names = substr(names(lnetp),1,1),filename = "ManuscriptFiles/Plots/Gender_agreement_Catalonia.png",fill = c("#B6453A", "#4A75A9","#19787F"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
    lista<-lnetdp
    lnetdp<-list("Women"=lista$Women,"Men"=lista$Men,"Adjusted"=lista$Adjusted)
    venn.diagram(x = lnetdp,category.names = substr(names(lnetdp),1,1),filename = "ManuscriptFiles/Plots/Gender_agreement_disgnificant_directionality_Catalonia.png",fill = c("#B6453A", "#4A75A9","#19787F"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
    
    ## Plot gender agreement Denmark ##
    soren<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
    soren<-soren[,c(1,2,5,6,7,8,13,14,19,20)]
    nwom<-soren[,3] ; nwom[which(is.na(nwom))]<-0
    nmen<-soren[,4] ; nmen[which(is.na(nmen))]<-0
    soren<-cbind(soren[,1:2],nwom+nmen,soren[,3:10])
    colnames(soren)<-c("A","B","nadj","nwom","nmen","rradj_min","rradj","rrwom_min","rrwom","rrmen_min","rrmen")
    
    lista<-list("Women"=paste(soren[rrwom_min >=1.01 & nwom >=100][[1]],soren[rrwom_min >=1.01 & nwom >=100][[2]],sep="_"),
                "Men"=paste(soren[rrmen_min >=1.01 & nmen >=100][[1]],soren[rrmen_min >=1.01 & nmen >=100][[2]],sep="_"),
                "Adjusted"=paste(soren[rradj_min >=1.01 & nadj >=100][[1]],soren[rradj_min >=1.01 & nadj >=100][[2]],sep="_"))
    venn.diagram(x = lista,category.names = substr(names(lista),1,1),filename = "ManuscriptFiles/Plots/Gender_agreement_Denmark.png",fill = c("#B6453A", "#4A75A9","#19787F"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
    
    #### Continuous_time_windows ####
    ficheros2<-c(ficheros[grep("N0_1",ficheros)],ficheros[grep("N1_2",ficheros)],ficheros[grep("N2_3",ficheros)],ficheros[grep("N3_4",ficheros)],ficheros[grep("N4_5",ficheros)])
    rangos<-unique(gsub("_Men.+","",gsub("_Wom.+","",gsub("_Adj.+","",ficheros2))))
    comparisons<-c("Women","Adjusted","Men")
    direction_information<-c() ; ldirection_information<-list() ; lsdirection_information<-list()
    vennfiles<-list() ; vennfiles_directional<-list()
    lupsetcontinuous<-list() ; lupsetcontinuoussigdir<-list()
    for(z in rangos){
      # z<-"N0_1"
      fichs<-sort(ficheros2[grep(paste(z,"_",sep=""),ficheros2)])
      lnetp<-list() ; lnet<-list() ; lnetdp<-list() ; lnetd<-list()
      for(a in fichs){
        # a<-"N0_5_Women.txt"
        ## Load the three tables ##
        bayesall<-fread(paste("ManuscriptFiles/Bayes/",a,sep=""),stringsAsFactors = F,sep="\t")
        
        ## Get the significant associations ##
        sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
        dsbayesall<-sbayesall[which(p.adjust(sbayesall$directionality_pval, method = "BH")<=0.05)]
        sbayesall<-sbayesall[,c(1,2,3,16)] ; colnames(sbayesall)<-c("A","B","number","rr")
        dsbayesall<-dsbayesall[,c(1,2,3,16)] ; colnames(dsbayesall)<-c("A","B","number","rr")
        lnet[[gsub(".txt","",gsub(".+_","",a))]]<-sbayesall
        lnetd[[gsub(".txt","",gsub(".+_","",a))]]<-dsbayesall
        lnetp[[gsub(".txt","",gsub(".+_","",a))]]<-paste(sbayesall$A,sbayesall$B,sep="_")
        lnetdp[[gsub(".txt","",gsub(".+_","",a))]]<-paste(dsbayesall$A,dsbayesall$B,sep="_")
        lupsetcontinuous[[gsub(".txt","",gsub(".+_","",a))]][[gsub("_","-",gsub("N","",z))]]<-paste(sbayesall$A,sbayesall$B,sep="_")
        lupsetcontinuoussigdir[[gsub(".txt","",gsub(".+_","",a))]][[gsub("_","-",gsub("N","",z))]]<-paste(dsbayesall$A,dsbayesall$B,sep="_")
      }
      vennfiles[[gsub("_","-",gsub("N","",z))]] <- venn.diagram(x = lnetp,category.names = substr(names(lnetp),1,1),filename = NULL,fill = c("#19787F", "#4A75A9","#B6453A"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
      ## Plot Venn Diagram with directionality ##
      vennfiles_directional[[gsub("_","-",gsub("N","",z))]] <- venn.diagram(x = lnetdp,category.names = substr(names(lnetdp),1,1),filename = NULL,fill = c("#19787F", "#4A75A9","#B6453A"),alpha = 0.6,cex = 0.8, main = paste(gsub("_"," - ",gsub("N","",z)),"time window"),disable.logging = TRUE)
      
      #### Detailed information on directionality ####
      Wom_pairs     <- lnetp$Women
      Men_pairs     <- lnetp$Men
      ## Comunes ##
      comunes<-intersect(Wom_pairs,Men_pairs)
      ## En ambas direcciones ##
      comundt<-cbind(gsub("_.+","",comunes),gsub(".+_","",comunes))
      scomundt<-t(apply(comundt, 1, sort))
      bidirectional<-length(which(duplicated(scomundt)))*2
      Wom_only<-setdiff(Wom_pairs,comunes)
      Men_only<-setdiff(Men_pairs,comunes)
      Wom_only_rev<-paste(gsub(".+_","",Wom_only),gsub("_.+","",Wom_only),sep="_")
      Men_only_rev<-paste(gsub(".+_","",Men_only),gsub("_.+","",Men_only),sep="_")
      WomonlyAB_Comrev<-intersect(Wom_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
      MenonlyAB_Comrev<-intersect(Men_only,paste(gsub(".+_","",comunes),gsub("_.+","",comunes),sep="_"))
      WomonlyAB_MenBArev<-intersect(Wom_only,Men_only_rev)
      MenonlyAB_WomBArev<-intersect(Men_only,Wom_only_rev)
      ## Save numbers ##
      direction_information<-rbind(direction_information,c(z,length(Wom_only),length(Men_only),length(comunes),bidirectional,length(WomonlyAB_MenBArev),length(MenonlyAB_WomBArev),length(WomonlyAB_Comrev),length(MenonlyAB_Comrev),
                                                           length(setdiff(ssba,intersect(ssba,ssep))),length(setdiff(ssep,intersect(ssba,ssep))),length(intersect(ssba,ssep))))
      ldirection_information[[z]]$WomonlyAB_MenBArev<-WomonlyAB_MenBArev
      ldirection_information[[z]]$MenonlyAB_WomBArev<-MenonlyAB_WomBArev
      
      #### Detailed information on significant directionality ####
      Wom_pairs     <- lnetdp$Women
      Men_pairs     <- lnetdp$Men
      ## Comunes ##
      comunes<-intersect(Wom_pairs,Men_pairs)
      Wom_only<-setdiff(Wom_pairs,comunes)
      Men_only<-setdiff(Men_pairs,comunes)
      Wom_only_rev<-paste(gsub(".+_","",Wom_only),gsub("_.+","",Wom_only),sep="_")
      Men_only_rev<-paste(gsub(".+_","",Men_only),gsub("_.+","",Men_only),sep="_")
      WomonlyAB_MenBArev<-intersect(Wom_only,Men_only_rev)
      MenonlyAB_WomBArev<-intersect(Men_only,Wom_only_rev)
      ## Save numbers ##
      lsdirection_information[[z]]$WomonlyAB_MenBArev<-WomonlyAB_MenBArev
      lsdirection_information[[z]]$MenonlyAB_WomBArev<-MenonlyAB_WomBArev
      print(paste(z,"finished!"))
    }
    colnames(direction_information)<-c("network","only_women","only_men","commons","bidirectional","only_women_men_reverse","only_men_women_reverse","only_women_common_reverse","only_men_common_reverse",
                                       "only_women_nodirection","only_men_nodirection","common_nodirection")
    ## Save agreements in table and list ##
    write.table(direction_information,"ManuscriptFiles/Results/Women_Men_overlap_comorbidities_continuous_timewindow.txt",quote=F,sep="\t",row.names=F)
    save(ldirection_information,lsdirection_information,file="ManuscriptFiles/IntermediateFiles/List_comorbidities_gender_agreement_continuous.RData")
    ## Plot the Venn Diagram ##
    pdf("ManuscriptFiles/Plots/Sex_agreement_comorbidities_continuous_timewindows.pdf", width = 20, height = 4)
    grid.arrange(grobs = list(grobTree(vennfiles$`0-1`),grobTree(vennfiles$`1-2`),grobTree(vennfiles$`2-3`),
                              grobTree(vennfiles$`3-4`),grobTree(vennfiles$`4-5`)),ncol = 5)
    dev.off()
    ## Plot the Venn Diagram no-direction ##
    pdf(paste("ManuscriptFiles/Plots/Sex_agreement_significant_directional_comorbidities_continuous_timewindows.pdf",sep=""), width = 20, height = 4)
    grid.arrange(grobs = list(grobTree(vennfiles_directional$`0-1`),grobTree(vennfiles_directional$`1-2`),grobTree(vennfiles_directional$`2-3`),
                              grobTree(vennfiles_directional$`3-4`),grobTree(vennfiles_directional$`4-5`)),ncol = 5)
    dev.off()
    save(lupsetcontinuous,lupsetcontinuoussigdir,lupsetincremental,lupsetincrementalsigdir,file="ManuscriptFiles/IntermediateFiles/Lists_for_upsetplot_shared_interactions.RData")
  }
  if(args[2]=="upsetplot_comorbidities_by_timewindow"){
    load(file="ManuscriptFiles/IntermediateFiles/Lists_for_upsetplot_shared_interactions.RData") # lupsetcontinuous,lupsetcontinuoussigdir,lupsetincremental,lupsetincrementalsigdir
    for(a in 1:length(names(lupsetcontinuous))){
      upset_data <- fromList(lupsetcontinuous[[a]])
      pdf(paste("ManuscriptFiles/Plots/UpsetTimeWindows/Continuous_timewindows_",names(lupsetcontinuous)[a],".pdf",sep=""),width = 10,height = 6)
        print(upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = ""))
      dev.off()
    }
    for(a in 1:length(names(lupsetcontinuoussigdir))){
      upset_data <- fromList(lupsetcontinuoussigdir[[a]])
      pdf(paste("ManuscriptFiles/Plots/UpsetTimeWindows/Continuous_timewindows_significantdirectionality_",names(lupsetcontinuoussigdir)[a],".pdf",sep=""),width = 10,height = 6)
        print(upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = ""))
      dev.off()
    }
    for(a in 1:length(names(lupsetincremental))){
      upset_data <- fromList(lupsetincremental[[a]])
      pdf(paste("ManuscriptFiles/Plots/UpsetTimeWindows/Incremental_timewindows_",names(lupsetincremental)[a],".pdf",sep=""),width = 10,height = 6)
      print(upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = ""))
      dev.off()
    }
    for(a in 1:length(names(lupsetincrementalsigdir))){
      upset_data <- fromList(lupsetincrementalsigdir[[a]])
      pdf(paste("ManuscriptFiles/Plots/UpsetTimeWindows/Incremental_timewindows_significantdirectionality_",names(lupsetincrementalsigdir)[a],".pdf",sep=""),width = 10,height = 6)
      print(upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = ""))
      dev.off()
    }
  }
}

#### Comparison Denmark-Catalonia focused on 0-5 time window: shared comorbidities, correlations in comorbidity OR, prevalence and disease prevalence ####
if(args[1]=="0_5_time_window_main_manuscript"){
  #### Correlations between number of comorbidities and disease prevalence ####
  if(args[2]=="comorbidity_disease_prevalence_correlations"){
    diseasesanalyzed<-fread("ManuscriptFiles/Number_patients_per_disease.txt",stringsAsFactors = F)
    prevalences<-fread("Data/DiseaseIncidencesGeneral.txt",stringsAsFactors = F,sep="\t")
    colnames(prevalences)[1]<-"diseases" ; setkey(prevalences,"diseases")
    prevalences<-prevalences[intersect(diseasesanalyzed$diseases,prevalences$diseases)]
    ficheros<-list.files("ManuscriptFiles/Bayes/")
    comparisons<-c("Adjusted","Women","Men")
    ## Catalonia ##
    ficheros<-ficheros[grep("N0_5",ficheros)]
    lesult<-list()
    for(z in comparisons){
      # z<-"Adjusted"
      fichs<-ficheros[grep(paste("_",z,sep=""),ficheros)]
      ## Load the table ##
      bayesall<-fread(paste("ManuscriptFiles/Bayes/",fichs,sep=""),stringsAsFactors = F,sep="\t")
      ## Get the significant associations ##
      sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
      #### Correlations ####
      if(z=="Adjusted"){patsbydisease<-prevalences[,1:2] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(z=="Women"){patsbydisease<-prevalences[,c(1,3)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      if(z=="Men"){patsbydisease<-prevalences[,c(1,4)] ; colnames(patsbydisease)[2]<-"patientwithcontrols"}
      ## significant bayes all ##
      ## @@ @@ @@ @ @ @@ @@ @@ ##
      ## Counts in Disease1 ##
      count_d1 <- sbayesall %>% count(A, name = "n_Disease1") %>% rename(diseases = A)
      ## Counts in Disease1 ##
      count_d2 <- sbayesall %>% count(B, name = "n_Disease2") %>% rename(diseases = B)
      ## Add to patsbydisease
      resulta <- patsbydisease %>% left_join(count_d1, by = "diseases") %>% left_join(count_d2, by = "diseases") %>%
        mutate(n_Disease1 = coalesce(n_Disease1, 0L),n_Disease2 = coalesce(n_Disease2, 0L))
      if(length(intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)))>0){resulta<-resulta[-intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)),]}
      lesult[[z]]<-resulta
      print(paste(z,"finished!"))
    }
    
    ## Denmark ##
    soren<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
    soren<-soren[,c(1,2,5,6,7,8,13,14,19,20)]
    nwom<-soren[,3] ; nwom[which(is.na(nwom))]<-0
    nmen<-soren[,4] ; nmen[which(is.na(nmen))]<-0
    soren<-cbind(soren[,1:2],nwom+nmen,soren[,3:10])
    colnames(soren)<-c("A","B","nadj","nwom","nmen","rradj_min","rradj","rrwom_min","rrwom","rrmen_min","rrmen")
    sorendis<-fread("Epidemiology/41467_2019_8475_MOESM4_ESM.txt",stringsAsFactors = F,sep="\t")
    sorendis<-cbind(sorendis[,c(1,3,4)],sorendis[[3]]+sorendis[[4]])
    colnames(sorendis)<-c("diseases","men","women","adjusted")
    comparisons<-c("Adjusted","Women","Men")
    sesult<-list()
    for(z in comparisons){
      # z<-"Adjusted"
      if(z=="Adjusted"){network<-soren[,c(1,2,3,6,7)] ; snetwork<-network[rradj_min >=1.01 & nadj >=100] ; colnames(snetwork)<-c("A","B","number","rrmin","rr")}
      if(z=="Women"){network<-soren[,c(1,2,4,8,9)] ; snetwork<-network[rrwom_min >=1.01 & nwom >=100] ; colnames(snetwork)<-c("A","B","number","rrmin","rr")}
      if(z=="Men"){network<-soren[,c(1,2,5,10,11)] ; snetwork<-network[rrmen_min >=1.01 & nmen >=100] ; colnames(snetwork)<-c("A","B","number","rrmin","rr")}
      
      #### Correlations ####
      if(z=="Adjusted"){patsbydisease<-sorendis[,c(1,4)] ; colnames(patsbydisease)<-c("diseases","patientwithcontrols")}
      if(z=="Women"){patsbydisease<-sorendis[,c(1,3)] ; colnames(patsbydisease)<-c("diseases","patientwithcontrols")}
      if(z=="Men"){patsbydisease<-sorendis[,c(1,2)] ; colnames(patsbydisease)<-c("diseases","patientwithcontrols")}
      ## significant bayes all ##
      ## @@ @@ @@ @ @ @@ @@ @@ ##
      ## Counts in Disease1 ##
      count_d1 <- snetwork %>% count(A, name = "n_Disease1") %>% rename(diseases = A)
      ## Counts in Disease1 ##
      count_d2 <- snetwork %>% count(B, name = "n_Disease2") %>% rename(diseases = B)
      ## Add to patsbydisease
      resulta <- patsbydisease %>% left_join(count_d1, by = "diseases") %>% left_join(count_d2, by = "diseases") %>%
        mutate(n_Disease1 = coalesce(n_Disease1, 0L),n_Disease2 = coalesce(n_Disease2, 0L))
      if(length(intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)))>0){resulta<-resulta[-intersect(which(resulta$n_Disease1==0),which(resulta$n_Disease2==0)),]}
      sesult[[z]]<-resulta
      print(paste(z,"finished!"))
    }
    
    ## Plot correlations between prevalence and comorbidities ##
    pdf(file="ManuscriptFiles/Plots/Correlation_number_patients_comorbidity_relationships_CatDen.pdf",width = 8,height = 22)
      par(mfrow = c(6, 2))  # 6 filas, 2 columnas
      ## Adjusted - Source - Catalonia ##
      plot(lesult$Adjusted$patientwithcontrols,lesult$Adjusted$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Both genders (source disease)\nCatalonia",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$Adjusted$patientwithcontrols,lesult$Adjusted$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$Adjusted$patientwithcontrols,lesult$Adjusted$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## Adjusted - Source - Denmark ##
      plot(sesult$Adjusted$patientwithcontrols,sesult$Adjusted$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Both genders (source disease)\nDenmark",
           xlim=c(0,1e6))
      text(6e5,150,paste("cor = ",round(cor.test(sesult$Adjusted$patientwithcontrols,sesult$Adjusted$n_Disease1)$estimate,3),
                         "\npval = ",round(cor.test(sesult$Adjusted$patientwithcontrols,sesult$Adjusted$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## Adjusted - Sink - Catalonia ##
      plot(lesult$Adjusted$patientwithcontrols,lesult$Adjusted$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Both genders (sink disease)\nCatalonia",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$Adjusted$patientwithcontrols,lesult$Adjusted$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$Adjusted$patientwithcontrols,lesult$Adjusted$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## Adjusted - Sink - Denmark ##
      plot(sesult$Adjusted$patientwithcontrols,sesult$Adjusted$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Both genders (sink disease)\nDenmark",
           xlim=c(0,1e6))
      text(6e5,150,paste("cor = ",round(cor.test(sesult$Adjusted$patientwithcontrols,sesult$Adjusted$n_Disease2)$estimate,3),
                         "\npval = ",round(cor.test(sesult$Adjusted$patientwithcontrols,sesult$Adjusted$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## Women - Source - Catalonia ##
      plot(lesult$Women$patientwithcontrols,lesult$Women$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Women (source disease)\nCatalonia",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$Women$patientwithcontrols,lesult$Women$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$Women$patientwithcontrols,lesult$Women$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## Women - Source - Denmark ##
      plot(sesult$Women$patientwithcontrols,sesult$Women$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Women (source disease)\nDenmark",
           xlim=c(0,1e6))
      text(6e5,100,paste("cor = ",round(cor.test(sesult$Women$patientwithcontrols,sesult$Women$n_Disease1)$estimate,3),
                         "\npval = ",round(cor.test(sesult$Women$patientwithcontrols,sesult$Women$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## Women - Sink - Catalonia ##
      plot(lesult$Women$patientwithcontrols,lesult$Women$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Women (sink disease)\nCatalonia",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$Women$patientwithcontrols,lesult$Women$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$Women$patientwithcontrols,lesult$Women$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## Women - Sink - Denmark ##
      plot(sesult$Women$patientwithcontrols,sesult$Women$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Women (sink disease)\nDenmark",
           xlim=c(0,1e6))
      text(6e5,100,paste("cor = ",round(cor.test(sesult$Women$patientwithcontrols,sesult$Women$n_Disease2)$estimate,3),
                         "\npval = ",round(cor.test(sesult$Women$patientwithcontrols,sesult$Women$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## Men - Source - Catalonia ##
      plot(lesult$Men$patientwithcontrols,lesult$Men$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Men (source disease)\nCatalonia",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$Men$patientwithcontrols,lesult$Men$n_Disease1)$estimate,3),
                        "\npval = ",round(cor.test(lesult$Men$patientwithcontrols,lesult$Men$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## Men - Source - Denmark ##
      plot(sesult$Men$patientwithcontrols,sesult$Men$n_Disease1,xlab="Number of patients",ylab="Number of comorbidities",main="Men (source disease)\nDenmark",
           xlim=c(0,1e6))
      text(6e5,100,paste("cor = ",round(cor.test(sesult$Men$patientwithcontrols,sesult$Men$n_Disease1)$estimate,3),
                         "\npval = ",round(cor.test(sesult$Men$patientwithcontrols,sesult$Men$n_Disease1)$p.value,6),sep=""),cex=0.7)
      ## Men - Sink - Catalonia ##
      plot(lesult$Men$patientwithcontrols,lesult$Men$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Men (sink disease)\nCatalonia",
           xlim=c(0,5e5))
      text(3e5,20,paste("cor = ",round(cor.test(lesult$Men$patientwithcontrols,lesult$Men$n_Disease2)$estimate,3),
                        "\npval = ",round(cor.test(lesult$Men$patientwithcontrols,lesult$Men$n_Disease2)$p.value,6),sep=""),cex=0.7)
      ## Men - Sink - Denmark ##
      plot(sesult$Men$patientwithcontrols,sesult$Men$n_Disease2,xlab="Number of patients",ylab="Number of comorbidities",main="Men (sink disease)\nDenmark",
           xlim=c(0,1e6))
      text(6e5,100,paste("cor = ",round(cor.test(sesult$Men$patientwithcontrols,sesult$Men$n_Disease2)$estimate,3),
                         "\npval = ",round(cor.test(sesult$Men$patientwithcontrols,sesult$Men$n_Disease2)$p.value,6),sep=""),cex=0.7)
    dev.off()
  }
  #### Overlaps with Soren Brunak's networks and correlations between shared OR and comorbidity prevalences ####
  if(args[2]=="comorbidity_overlap_and_correlations_Catalonia_vs_Denmark"){
    soren<-fread("Epidemiology/41467_2019_8475_MOESM6_ESM.txt",stringsAsFactors = F,sep="\t")
    soren<-soren[,c(1,2,5,6,7,8,13,14,19,20)]
    nwom<-soren[,3] ; nwom[which(is.na(nwom))]<-0
    nmen<-soren[,4] ; nmen[which(is.na(nmen))]<-0
    soren<-cbind(soren[,1:2],nwom+nmen,soren[,3:10])
    colnames(soren)<-c("A","B","nadj","nwom","nmen","rradj_min","rradj","rrwom_min","rrwom","rrmen_min","rrmen")
    comparisons<-c("Adjusted","Women","Men")
    ficheros<-ficheros[grep("N0_5",ficheros)]
    lesult<-list()
    vennfiles<-list()
    for(z in comparisons){
      # z<-"Adjusted"
      fichs<-ficheros[grep(paste("_",z,sep=""),ficheros)]
      ## Load the three tables ##
      bayesall<-fread(paste("ManuscriptFiles/Bayes/",fichs,sep=""),stringsAsFactors = F,sep="\t")
      ## Select the epidemiology of interest ##
      if(z=="Adjusted"){epidemiology<-soren[,c(1:3,6,7)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
      if(z=="Women"){epidemiology<-soren[,c(1:2,4,8,9)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
      if(z=="Men"){epidemiology<-soren[,c(1:2,5,10,11)] ; colnames(epidemiology)<-c("A","B","number","rrmin","rr")}
      ## Remove NAs and select the number of nodes ##
      if(length(which(is.na(epidemiology$rrmin)))>0){epidemiology<-epidemiology[-which(is.na(epidemiology$rrmin))]}
      
      ## Get the nodes ##
      bayesallnodes<-unique(c(bayesall$A,bayesall$B))
      epidemiologynodes<-unique(c(epidemiology$A,epidemiology$B))
      
      ## Common diseases analyzed ##
      bayesallcommon<-intersect(bayesallnodes,epidemiologynodes)
      
      ## Select comorbidities between shared nodes in Catalonia ##
      bayesall <- bayesall[bayesall$A %in% bayesallcommon & bayesall$B %in% bayesallcommon,]
      
      ## Select comorbidities between shared nodes in Denmark ##
      epidemiology <- epidemiology[epidemiology$A %in% bayesallcommon & epidemiology$B %in% bayesallcommon,]
      
      ## Get the significant associations ##
      sbayesall<-bayesall[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
      sbayesall<-sbayesall[,c(1,2,3,16)] ; colnames(sbayesall)<-c("A","B","number","rr")
      sepidemiology<-epidemiology[rrmin >=1.01 & number >=100,c(1:3,5)]
      
      ## Plot Venn Diagram of shared interactions ##
      sba<-paste(sbayesall$A,sbayesall$B,sep="_")
      sep<-paste(sepidemiology$A,sepidemiology$B,sep="_")
      vennfiles[[z]] <- venn.diagram(x = list(sba, sep),category.names = c("C", "D"),filename = NULL,fill = c("#19787F", "#440F53"),alpha = 0.6,cex = 0.8, main = paste(z,paste("JI =",round(overlapmetrics(sba, sep)[1],4)),paste("C in D =",round(overlapmetrics(sba, sep)[2]*100,2)),paste("D in C =",round(overlapmetrics(sba, sep)[3]*100,2)),sep="\n"),disable.logging = TRUE)
      
      ## Get the number of patients with diseases and the number of comorbidities ##
      setkey(sbayesall, A, B)
      setkey(sepidemiology, A, B)
      lesult[[z]] <- sbayesall[sepidemiology, nomatch = 0]
      print(paste(z,"finished!"))
    }
    ## Get the number of individuals with each of the diseases ##
    ## Catalonia ##
    diseasesanalyzed<-fread("ManuscriptFiles/Number_patients_per_disease.txt",stringsAsFactors = F)
    prevalences<-fread("Data/DiseaseIncidencesGeneral.txt",stringsAsFactors = F,sep="\t")
    colnames(prevalences)[1]<-"diseases" ; setkey(prevalences,"diseases")
    prevalences<-prevalences[intersect(diseasesanalyzed$diseases,prevalences$diseases)]
    setkey(prevalences,"diseases")
    ## Denmark ##
    sorendis<-fread("Epidemiology/41467_2019_8475_MOESM4_ESM.txt",stringsAsFactors = F,sep="\t")
    sorendis<-cbind(sorendis[,c(1,3,4)],sorendis[[3]]+sorendis[[4]])
    colnames(sorendis)<-c("diseases","Men","Women","General")
    setkey(sorendis,"diseases")
    ## Common diseases ##
    commondis<-intersect(prevalences$diseases,sorendis$diseases)
    ## Select comorbidities between shared nodes in Catalonia ##
    catalonia <- prevalences[commondis]
    ## Select comorbidities between shared nodes in Catalonia ##
    denmark <- sorendis[commondis]
    
    
    ## Plot the Venn Diagram ##
    pdf(paste("ManuscriptFiles/Plots/VennDiagram_shared_comorbidities.pdf"), width = 12, height = 4)
    grid.arrange(grobs = list(grobTree(vennfiles$Adjusted),grobTree(vennfiles$Women),grobTree(vennfiles$Men)),ncol = 3)
    dev.off()
    
    ## Plot the correlations ##
    pdf(file="ManuscriptFiles/Plots/Correlations_comorbidities_Catalonia_Denmark.pdf",width = 12,height = 12)
      par(mfrow = c(3, 3))
      ## 1 ##
      plot(lesult$Adjusted$number,lesult$Adjusted$i.number,xlab="Catalonia",ylab="Denmark",main="Both genders\nComorbidity prevalence")
      text(max(lesult$Adjusted$number)/3*2,max(lesult$Adjusted$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$Adjusted$number,lesult$Adjusted$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$Adjusted$number,lesult$Adjusted$i.number)$p.value,6),sep=""),cex=0.7)
      ## 2 ##
      plot(lesult$Women$number,lesult$Women$i.number,xlab="Catalonia",ylab="Denmark",main="Women\nComorbidity prevalence")
      text(max(lesult$Women$number)/3*2,max(lesult$Women$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$Women$number,lesult$Women$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$Women$number,lesult$Women$i.number)$p.value,6),sep=""),cex=0.7)
      ## 3 ##
      plot(lesult$Men$number,lesult$Men$i.number,xlab="Catalonia",ylab="Denmark",main="Men\nComorbidity prevalence")
      text(max(lesult$Men$number)/3*2,max(lesult$Men$i.number)/3*2,
           paste("cor = ",round(cor.test(lesult$Men$number,lesult$Men$i.number)$estimate,3),
                 "\npval = ",round(cor.test(lesult$Men$number,lesult$Men$i.number)$p.value,6),sep=""),cex=0.7)
      ## 4 ##
      plot(lesult$Adjusted$rr,lesult$Adjusted$i.rr,xlab="Catalonia",ylab="Denmark",main="Both genders\nComorbidity risk")
      text(max(lesult$Adjusted$rr)/3*2,max(lesult$Adjusted$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$Adjusted$rr,lesult$Adjusted$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$Adjusted$rr,lesult$Adjusted$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 5 ##
      plot(lesult$Women$rr,lesult$Women$i.rr,xlab="Catalonia",ylab="Denmark",main="Women\nComorbidity risk")
      text(max(lesult$Women$rr)/3*2,max(lesult$Women$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$Women$rr,lesult$Women$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$Women$rr,lesult$Women$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 6 ##
      plot(lesult$Men$rr,lesult$Men$i.rr,xlab="Catalonia",ylab="Denmark",main="Men\nComorbidity risk")
      text(max(lesult$Men$rr)/3*2,max(lesult$Men$i.rr)/3*2,
           paste("cor = ",round(cor.test(lesult$Men$rr,lesult$Men$i.rr)$estimate,3),
                 "\npval = ",round(cor.test(lesult$Men$rr,lesult$Men$i.rr)$p.value,6),sep=""),cex=0.7)
      ## 7 ##
      plot(catalonia$General,denmark$General,xlab="Catalonia",ylab="Denmark",main="Both genders\nDisease prevalence")
      text(max(catalonia$General)/3*2,max(denmark$General)/3*2,
           paste("cor = ",round(cor.test(catalonia$General,denmark$General)$estimate,3),
                 "\npval = ",round(cor.test(catalonia$General,denmark$General)$p.value,6),sep=""),cex=0.7)
      ## 8 ##
      plot(catalonia$Female,denmark$Women,xlab="Catalonia",ylab="Denmark",main="Women\nDisease prevalence")
      text(max(catalonia$Female)/3*2,max(denmark$Women)/3*2,
           paste("cor = ",round(cor.test(catalonia$Female,denmark$Women)$estimate,3),
                 "\npval = ",round(cor.test(catalonia$Female,denmark$Women)$p.value,6),sep=""),cex=0.7)
      ## 9 ##
      plot(catalonia$Male,denmark$Men,xlab="Catalonia",ylab="Denmark",main="Men\nDisease prevalence")
      text(max(catalonia$Male)/3*2,max(denmark$Men)/3*2,
           paste("cor = ",round(cor.test(catalonia$Male,denmark$Men)$estimate,3),
                 "\npval = ",round(cor.test(catalonia$Male,denmark$Men)$p.value,6),sep=""),cex=0.7)
    dev.off()
  }
  #### Look for ICD10 categories enriched in overlapping comorbidities ####
  if(args[2]=="enrichment_overlapping_comorbidities"){
    abayes<-fread("ManuscriptFiles/Bayes/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
    wbayes<-fread("ManuscriptFiles/Bayes/N0_5_Women.txt",stringsAsFactors = F,sep="\t")
    mbayes<-fread("ManuscriptFiles/Bayes/N0_5_Men.txt",stringsAsFactors = F,sep="\t")
  }
}

#### Analyze tendencies in comorbidities by OR changes ####
if(args[3]=="compare_time_windows"){
  ficheros<-list.files("ManuscriptFiles/Bayes/")
  ficheros<-ficheros[grep("N0_",ficheros)]
  comparisons<-c("Adjusted","Women","Men")
  lnet<-list()
  lnet2<-list()
  lall<-list()
  for(a in comparisons){
    subfich<-ficheros[grep(paste("_",a,sep=""),ficheros)]
    for(b in subfich){
      tt<-fread(paste("ManuscriptFiles/Bayes/",b,sep=""),stringsAsFactors = F,sep="\t")
      tt<-tt[,c(1:3,13:18,20)]
      lall[[a]][[gsub(paste("_",a,".txt",sep=""),"",b)]]<-tt
      if(gsub("_.+","",gsub("N0_","",b))==1){athres<-20}
      if(gsub("_.+","",gsub("N0_","",b))==2){athres<-40}
      if(gsub("_.+","",gsub("N0_","",b))==3){athres<-60}
      if(gsub("_.+","",gsub("N0_","",b))==4){athres<-80}
      if(gsub("_.+","",gsub("N0_","",b))==5){athres<-100}
      stt1<-tt[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=athres]
      # print(paste(gsub(paste("_",a,".txt",sep=""),"",b),": ",length(stt1$A),sep=""))
      stt2<-tt[lfsr < 0.05 & abs(log_or_shrunk) > log(1.3) & a>=100]
      print(paste(gsub(paste("_",a,".txt",sep=""),"",b),": ",length(stt2$A),sep=""))
      lnet[[a]][[gsub(paste("_",a,".txt",sep=""),"",b)]]<-stt2
      lnet2[[a]][[gsub(paste("_",a,".txt",sep=""),"",b)]]<-stt1
    }
    print(a)
  }
  
  ## Plot the overlaps between the different time windows ##
  ## Venn Diagram ##
  lista<-list("0-1"=paste(lnet$Adjusted$N0_1$A,lnet$Adjusted$N0_1$B,sep="_"),
              "0-2"=paste(lnet$Adjusted$N0_2$A,lnet$Adjusted$N0_2$B,sep="_"),
              "0-3"=paste(lnet$Adjusted$N0_3$A,lnet$Adjusted$N0_3$B,sep="_"),
              "0-4"=paste(lnet$Adjusted$N0_4$A,lnet$Adjusted$N0_4$B,sep="_"),
              "0-5"=paste(lnet$Adjusted$N0_5$A,lnet$Adjusted$N0_5$B,sep="_"))
  venn.diagram(lista,filename = "ManuscriptFiles/Plots/Overlap_comorbidities_vaerying_time_window_VennDiagram.png",alpha = 0.6,cex = 0.8,disable.logging = TRUE)
  ## Upset plot ##
  upset_data <- fromList(lista)
  pdf("ManuscriptFiles/Plots/Overlap_comorbidities_vaerying_time_window_UpSetPlot.pdf",width = 10,height = 6)
  upset(upset_data,order.by = "freq",decreasing = TRUE,mainbar.y.label = "",sets.x.label = "")
  dev.off()
  
  #### Temporal stability of comorbidity relationships ####
  windows <- names(lnet$Adjusted)  # N0_1 ... N0_5
  ## Put all the time windows together ##
  dt_long <- rbindlist(
    lapply(seq_along(windows), function(i) {
      dt <- copy(lnet$Adjusted[[windows[i]]])
      dt[, window := i]           # ventana temporal (1–5)
      dt
    }),
    use.names = TRUE,
    fill = TRUE
  )
  
  ## Indicate if they have already been filtered ##
  dt_long[, present := 1L]
  
  ## Agregate the presence by pairs, so that each pair appears just once and the windows in which is present is shown ##
  stability <- dt_long[,.(n_windows = .N,windows_present = list(sort(window)),min_window = min(window),max_window = max(window),mean_log_or = mean(log_or_shrunk),sd_log_or   = sd(log_or_shrunk)),by = .(A, B)]
  
  patternvec<-sapply(stability$windows_present, function(x) paste(x, collapse = ","))
  temporal_pattern<-rep("other",length(stability$A))
  temporal_pattern[which(patternvec=="1")]<-"early-onset (one year)"
  temporal_pattern[which(patternvec=="1,2")]<-"early-onset (two years)"
  temporal_pattern[which(patternvec=="1,2,3")]<-"early-onset (three years)"
  temporal_pattern[which(patternvec=="5")]<-"late-onset (one year)"
  temporal_pattern[which(patternvec=="4,5")]<-"late-onset (two years)"
  temporal_pattern[which(patternvec=="3,4,5")]<-"late-onset (three years)"
  temporal_pattern[which(patternvec=="2,3,4")]<-"transient"
  temporal_pattern[which(patternvec=="1,2,3,4,5")]<-"persistent"
  
  stability<-cbind(stability,temporal_pattern)
  
  ## Summarize temporal patterns ##
  pattern_summary <- stability[, .N, by = temporal_pattern]
  
  print(pattern_summary)
  
  ## more representative examples by pattern ##
  top_early_only <- stability[
    temporal_pattern == "early-onset (one year)"
  ][order(-mean_log_or)][1:10]
  
  top_persistent <- stability[
    temporal_pattern == "persistent"
  ][order(-mean_log_or)][1:10]
  
  top_late_onset <- stability[
    temporal_pattern == "late-onset (one year)"
  ][order(-mean_log_or)][1:10]
  
  print(top_early_only)
  print(top_persistent)
  print(top_late_onset)
  
  
  #### Analyze the obtained results ####
  ## Are there diseases with significant changes in the number of significant comorbidities by time window? ##
  fit_nb_safe <- function(dt) {
    ## If there is no variation it cannot be computed ##
    if (length(unique(dt$count)) == 1) {
      return(list(beta = NA_real_,se   = NA_real_,z    = NA_real_,pval = NA_real_))
    }
    fit <- tryCatch(
      glm.nb(count ~ year, data = dt, control = glm.control(maxit = 50)),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(beta = NA_real_,se   = NA_real_,z    = NA_real_,pval = NA_real_))
    }
    coefs <- summary(fit)$coefficients
    list(
      beta = coefs["year", "Estimate"],
      se   = coefs["year", "Std. Error"],
      z    = coefs["year", "z value"],
      pval = coefs["year", "Pr(>|z|)"]
    )
  }
  ## Disease A ##
  ## @@ @ @ @@ ##
  comorb_by_A_window <- dt_long[,.N,by = .(A, window)]
  ## Widefy ##
  comorb_wide_A <- dcast(comorb_by_A_window,A ~ window,value.var = "N",fill = 0)
  ## Order ##
  comorb_wide_A<-comorb_wide_A[order(comorb_wide_A$`1`-comorb_wide_A$`5`,decreasing=T),]
  ## Put it in long mode and make the year an integer ##
  comorb_long_A <- melt(comorb_wide_A,id.vars = "A",variable.name = "year",value.name = "count")
  comorb_long_A[, year := as.integer(year)]
  res_nb_A <- comorb_long_A[, fit_nb_safe(.SD), by = A]
  res_nb_A[, p_adj := p.adjust(pval, method = "BH")]
  res_nb_A[, direction := fifelse(beta > 0, "Increase", "Decrease")]
  res_nb_A[, pct_change_per_year := (exp(beta) - 1) * 100]
  sig_nb_A <- res_nb_A[!is.na(p_adj) & p_adj < 0.05]
  sig_nb_A<-sig_nb_A[order(abs(sig_nb_A$beta),decreasing = T)]
  ## Add the number of comorbidities ##
  sig_nb_A_full <- merge(sig_nb_A,comorb_wide_A,by = "A",all.x = TRUE)
  sig_nb_A_full<-cbind(sig_nb_A_full,sig_nb_A_full$beta*abs(sig_nb_A_full$`5`-sig_nb_A_full$`1`))
  colnames(sig_nb_A_full)[14]<-c("beta*diff")
  sig_nb_A_full<-sig_nb_A_full[order(sig_nb_A_full$`beta*diff`,decreasing = T)]
  sig_nb_A_full<-cbind(sig_nb_A_full,as.character(distocol[sig_nb_A_full$A]))
  colnames(sig_nb_A_full)[15]<-"color"
  ## Plot tendencies ##
  ## Increasing ##
  sig_nb_A_full_i<-sig_nb_A_full[which(sig_nb_A_full$direction=="Increase")]
  plot_data_A_increase <- melt(
    sig_nb_A_full_i,
    id.vars = c("A", "color"),
    measure.vars = c("1","2","3","4","5"),
    variable.name = "year",
    value.name = "count"
  )
  plot_data_A_increase[, year := as.integer(as.character(year))]
  ## Decreasing ##
  sig_nb_A_full_d<-sig_nb_A_full[which(sig_nb_A_full$direction=="Decrease")]
  plot_data_A_decrease <- melt(
    sig_nb_A_full_d,
    id.vars = c("A", "color"),
    measure.vars = c("1","2","3","4","5"),
    variable.name = "year",
    value.name = "count"
  )
  plot_data_A_decrease[, year := as.integer(as.character(year))]
  p_decrease_A <- ggplot(plot_data_A_decrease, 
                         aes(x = year, y = count, group = A, color = color)) +
    geom_line(alpha = 0.6, size = 0.8) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_identity() +   # usa los códigos hex tal cual
    labs(x = "Time window duration (years)",
         y = "Number of comorbid conditions",
         title = "Decreasing trends") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  p_increase_A <- ggplot(plot_data_A_increase, 
                         aes(x = year, y = count, group = A, color = color)) +
    geom_line(alpha = 0.6, size = 0.8) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_identity() +
    labs(x = "Time window duration (years)",
         y = "Number of comorbid conditions",
         title = "Increasing trends") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  
  ## Disease B ##
  ## @@ @ @ @@ ##
  comorb_by_B_window <- dt_long[,.N,by = .(B, window)]
  ## Widefy ##
  comorb_wide_B <- dcast(comorb_by_B_window,B ~ window,value.var = "N",fill = 0)
  ## Order ##
  comorb_wide_B<-comorb_wide_B[order(comorb_wide_B$`1`-comorb_wide_B$`5`,decreasing=T),]
  ## Put it in long mode and make the year an integer ##
  comorb_long_B <- melt(comorb_wide_B,id.vars = "B",variable.name = "year",value.name = "count")
  comorb_long_B[, year := as.integer(year)]
  res_nb_B <- comorb_long_B[, fit_nb_safe(.SD), by = B]
  res_nb_B[, p_Adj := p.adjust(pval, method = "BH")]
  res_nb_B[, direction := fifelse(beta > 0, "Increase", "Decrease")]
  res_nb_B[, pct_change_per_year := (exp(beta) - 1) * 100]
  sig_nb_B <- res_nb_B[!is.na(p_Adj) & p_Adj < 0.05]
  sig_nb_B<-sig_nb_B[order(abs(sig_nb_B$beta),decreasing = T)]
  ## Add the number of comorbidities ##
  sig_nb_B_full <- merge(sig_nb_B,comorb_wide_B,by = "B",all.x = TRUE)
  sig_nb_B_full<-cbind(sig_nb_B_full,sig_nb_B_full$beta*abs(sig_nb_B_full$`5`-sig_nb_B_full$`1`))
  colnames(sig_nb_B_full)[14]<-c("beta*diff")
  sig_nb_B_full<-sig_nb_B_full[order(sig_nb_B_full$`beta*diff`,decreasing = T)]
  sig_nb_B_full<-cbind(sig_nb_B_full,as.character(distocol[sig_nb_B_full$B]))
  colnames(sig_nb_B_full)[15]<-"color"
  
  ## Plot tendencies ##
  ## Increasing ##
  sig_nb_B_full_i<-sig_nb_B_full[which(sig_nb_B_full$direction=="Increase")]
  plot_data_B_increase <- melt(
    sig_nb_B_full_i,
    id.vars = c("B", "color"),
    measure.vars = c("1","2","3","4","5"),
    variable.name = "year",
    value.name = "count"
  )
  plot_data_B_increase[, year := as.integer(as.character(year))]
  ## Decreasing ##
  sig_nb_B_full_d<-sig_nb_B_full[which(sig_nb_B_full$direction=="Decrease")]
  plot_data_B_decrease <- melt(
    sig_nb_B_full_d,
    id.vars = c("B", "color"),
    measure.vars = c("1","2","3","4","5"),
    variable.name = "year",
    value.name = "count"
  )
  plot_data_B_decrease[, year := as.integer(as.character(year))]
  p_decrease_B <- ggplot(plot_data_B_decrease, 
                         aes(x = year, y = count, group = B, color = color)) +
    geom_line(alpha = 0.6, size = 0.8) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_identity() +   # usa los códigos hex tal cual
    labs(x = "Time window duration (years)",
         y = "Number of comorbid conditions",
         title = "Decreasing trends") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  p_increase_B <- ggplot(plot_data_B_increase, 
                         aes(x = year, y = count, group = B, color = color)) +
    geom_line(alpha = 0.6, size = 0.8) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_identity() +
    labs(x = "Time window duration (years)",
         y = "Number of comorbid conditions",
         title = "Increasing trends") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
  
  ## Plot them together ##
  p_combined<-p_decrease_A + p_increase_A + p_decrease_B + p_increase_B + plot_layout(nrow=2,ncol=2)
  
  pdf(file="ManuscriptFiles/Plots/Number_of_comorbidities_change_by_timewindow.pdf",width = 12,height = 9)
  p_combined
  dev.off()
  
  sig_nb_B_full_d<-sig_nb_B_full_d[order(sig_nb_B_full_d$`1`,decreasing = T)]
  sig_nb_A_full_d<-sig_nb_A_full_d[order(sig_nb_A_full_d$`1`,decreasing = T)]
  sig_nb_B_full_i<-sig_nb_B_full_i[order(sig_nb_B_full_i$`5`,decreasing = T)]
  sig_nb_A_full_i<-sig_nb_A_full_i[order(sig_nb_A_full_i$`5`,decreasing = T)]
  
  head(sig_nb_A_full_d)
  head(as.character(tdig3dis[sig_nb_B_full_d$B]))
  head(sig_nb_B_full_d)
  head(sig_nb_A_full_i)
  head(sig_nb_B_full_i)
  
  ## Add the name ##
  dig3dis<-fread("Data/ICD10_three_digits_names.txt",stringsAsFactors = F)
  tdig3dis<-dig3dis$Name ; names(tdig3dis)<-dig3dis$Code
  
  ## Which is the jaccard index of increasing comorbidities as source and sink? and on decreasing?
  length(intersect(sig_nb_A_full_i$A,sig_nb_B_full_i$B))/length(c(sig_nb_A_full_i$A,sig_nb_B_full_i$B))
  length(intersect(sig_nb_A_full_d$A,sig_nb_B_full_d$B))/length(c(sig_nb_A_full_d$A,sig_nb_B_full_d$B))
  
  
  
  
  
  ## Enrichment analysis ##
  changelist<-list("A_significant_increase"=sig_nb_A_full_i$A,"A_significant_decrease"=sig_nb_A_full_d$A,
                   "B_significant_increase"=sig_nb_B_full_i$B,"B_significant_decrease"=sig_nb_B_full_d$B)
  
  universe<-unique(c(dt_long$A,dt_long$B))
  
  catename
  
  
  
  library(data.table)
  
  # Función que calcula enriquecimiento por lista de enfermedades
  enrichment_test <- function(sig_vector, universe, catename) {
    
    # solo usar enfermedades en el universo
    sig_vector <- sig_vector[sig_vector %in% universe]
    universe <- universe[universe %in% names(catename)]
    
    # categorias presentes
    categories <- unique(catename[universe])
    
    res <- lapply(categories, function(cat) {
      
      # A: número en lista y en categoría
      a <- sum(sig_vector %in% names(catename)[catename == cat])
      # B: número fuera de lista y en categoría
      b <- sum((universe %in% names(catename)[catename == cat]) & !(universe %in% sig_vector))
      # C: número en lista y fuera de categoría
      c <- sum(sig_vector %in% universe & catename[sig_vector] != cat)
      # D: número fuera de lista y fuera de categoría
      d <- length(universe) - a - b - c
      
      mat <- matrix(c(a, b, c, d), nrow=2)
      
      ft <- fisher.test(mat, alternative = "greater")  # test de enriquecimiento
      
      data.table(
        category = cat,
        count_in_list = a,
        count_in_universe = sum(universe %in% names(catename)[catename == cat]),
        pvalue = ft$p.value
      )
    })
    
    res <- rbindlist(res)
    # ajuste BH
    res[, p_adj := p.adjust(pvalue, method = "BH")]
    res <- res[order(p_adj)]
    return(res)
  }
  
  # Lista de aumentos en A
  enrichment_A_increase <- enrichment_test(
    sig_vector = changelist$A_significant_increase,
    universe = universe,
    catename = catename
  )
  
  # Lista de disminuciones en A
  enrichment_A_decrease <- enrichment_test(
    sig_vector = changelist$A_significant_decrease,
    universe = universe,
    catename = catename
  )
  
  # Puedes repetir para B
  enrichment_B_increase <- enrichment_test(
    sig_vector = changelist$B_significant_increase,
    universe = universe,
    catename = catename
  )
  
  enrichment_B_decrease <- enrichment_test(
    sig_vector = changelist$B_significant_decrease,
    universe = universe,
    catename = catename
  )
  
  
}


#### Caracterizar las redes de las distintas ventanas temporales: grado, betweeness, conexiones, nodos, ####














#### Concatenate comorbidities to get trajectories ####
## Nota: AllGlobal tiene todas las parejas de enfermedades, mientras que Global tiene solo las que tienen una co-ocurrencia significativa (FDR<=0.05)

if(args[1]=="trajectories"){
  if(args[2]=="get_patients_in_comorbidity_relationships"){
    if("Pats_by_comorbidity"%in%list.files("Data/")==FALSE){dir.create("Data/Pats_by_comorbidity")}
    ## Load comorbidity relationships
    comorbidities<-fread("FinalNetworks/Global/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
    sigcomdir<-comorbidities[which(comorbidities$binomial_FDR<=0.05)]
    sigcomdir<-sigcomdir[order(sigcomdir$CasesWith,decreasing = T)]
    ## Identify the number of patients with disease 1 and later on disease 2, difference not larger than 5 years ##
    npats<-rep(NA,length(sigcomdir$Disease1))
    # lpats<-list()
    for(a in 1:length(sigcomdir$Disease1)){
      # a<-1
      dis1<-sigcomdir$Disease1[a] ; dis2<-sigcomdir$Disease2[a]
      d1<-fread(paste("Data/DiseaseDiagnoses5years/",dis1,".txt",sep="")) ; setkey(d1,"idp")
      d2<-fread(paste("Data/DiseaseDiagnoses/",dis2,".txt",sep="")) ; setkey(d2,"idp")
      d1<-d1[intersect(d1$idp,d2$idp)] ; d2<-d2[intersect(d1$idp,d2$idp)]
      pats<-d1$idp[intersect(which(d1$dat<d2$dat),which(d2$dat-d1$dat<=50000))]
      npats[a]<-length(pats)
      # lpats[[paste(dis1,dis2,sep="_")]]<-pats
      write.table(pats,paste("Data/Pats_by_comorbidity/",dis1,"_",dis2,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names = F)
      print(paste(round((a/length(sigcomdir$Disease1))*100,2),"%",sep=""))
    }
    tabla<-cbind(sigcomdir$Disease1,sigcomdir$Disease2,npats)
    colnames(tabla)<-c("Disease1","Disease2","CasesWith")
    write.table(npats,"Data/Number_patients_with_comorbidities.txt",quote=F,sep="\t",row.names=F)
  }
  if(args[2]=="get_trajectories"){
    if(args[3]=="build_list"){
      ## Load comorbidity relationships
      comorbidities<-fread("FinalNetworks/Global/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
      sigcomdir<-comorbidities[which(comorbidities$binomial_FDR<=0.05)]
      sigcomdir<-sigcomdir[order(sigcomdir$CasesWith,decreasing = T)]
      ## load number of patients in each comorbidity ##
      numbers<-fread("Data/Number_patients_with_comorbidities.txt",stringsAsFactors = F,sep="\t")
      casenumbers<-cbind(sigcomdir[,1:2],as.numeric(numbers$x),paste(sigcomdir$Disease1,"_",sigcomdir$Disease2,".txt",sep="")) ; colnames(casenumbers)[3:4]<-c("CasesWith","filename")
      write.table(casenumbers,"Data/Number_of_patients_in_directional_comorbidities.txt",quote=F,sep="\t",row.names=F)
      ## load patients in comorbid relationships ##
      files<-list.files("Data/Pats_by_comorbidity/")
      lcomorbidities<-list()
      for(a in 1:length(casenumbers$filename)){
        # a<-1
        tt<-fread(paste("Data/Pats_by_comorbidity/",casenumbers$filename[a],sep=""),stringsAsFactors = F,header=F)$V1
        lcomorbidities[[gsub(".txt","",casenumbers$filename[a])]]<-tt
        print(paste(round((a/length(casenumbers$filename))*100,2),"%",sep=""))
      }
      save(lcomorbidities,file="Data/Patients_by_significant_directional_comorbidities.RData")
    }
    if(args[3]=="get_length_3"){
      ## Get 3-steps long trajectories ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## Load the list indicating the patients suffering specific comorbidity relationships ##
      load("Data/Patients_by_significant_directional_comorbidities.RData") # lcomorbidities
      ## Load the number of patients suffering from comorbidities ##
      casenumbers<-fread("Data/Number_of_patients_in_directional_comorbidities.txt",stringsAsFactors = F)
      ## load patients in comorbid relationships ##
      files<-list.files("Data/Pats_by_comorbidity/")
      lengththree<-c()
      lengththreel<-list()
      for(a in 1:length(casenumbers$filename)){
        # a<-1
        second<-casenumbers$Disease2[a]
        cuales<-which(casenumbers$Disease1==second)
        if(length(cuales)>0){
          for(b in cuales){
            # b<-19
            cuantos<-length(intersect(lcomorbidities[[a]],lcomorbidities[[b]]))
            lengththree<-rbind(lengththree,c(paste(gsub(".txt","-",casenumbers$filename[a]),gsub(".+_","",names(lcomorbidities)[b]),sep=""),casenumbers$CasesWith[a],cuantos))
            lengththreel[[paste(gsub(".txt","-",casenumbers$filename[a]),gsub(".+_","",names(lcomorbidities)[b]),sep="")]]<-intersect(lcomorbidities[[a]],lcomorbidities[[b]])
          }
        }
        print(paste(round((a/length(casenumbers$filename))*100,2),"%",sep=""))
      }
      save(lengththreel,file="Data/Patients_by_3_steps_trajectories.RData")
      colnames(lengththree)<-c("trajectory","patients_in_two","patients_in_three")
      write.table(lengththree,"Data/Number_of_patients_in_3_steps_long_trajectories.txt",quote=F,sep="\t",row.names=F)
    }
    if(args[3]=="get_length_4"){
      ini<-Sys.time()
      ## Get 4-steps long trajectories ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## Load the number of patients following a trajectory ##
      lengththree<-fread("Data/Number_of_patients_in_3_steps_long_trajectories.txt",stringsAsFactors = F)
      lengththree<-lengththree[which(lengththree$patients_in_three>=20)]
      ## Load the list indicating the patients following three steps long trajectories ##
      load("Data/Patients_by_3_steps_trajectories.RData") # lengththreel
      ## Load the list indicating the patients suffering specific comorbidity relationships ##
      load("Data/Patients_by_significant_directional_comorbidities.RData") # lcomorbidities
      lengthfourl<-list()
      for(a in 1:length(lengththree$trajectory)){
        # a<-1
        second<-gsub(".+-","",lengththree$trajectory[a])
        cuales<-grep(paste(second,"_",sep=""),names(lcomorbidities))
        if(length(cuales)>0){
          for(b in cuales){
            # b<-252
            cuantos<-length(intersect(lengththreel[[lengththree$trajectory[a]]],lcomorbidities[[b]]))
            if(cuantos>=20){
              recuento<-c(paste(gsub("-","_",lengththree$trajectory[a]),gsub(".+_","",names(lcomorbidities)[b]),sep="-"),lengththree[a,2:3],cuantos)
              write.table(recuento,"Data/Number_of_patients_in_4_steps_long_trajectories.txt",quote=F,sep="\t",row.names=F,col.names = F,append = T)
              lengthfourl[[paste(gsub("-","_",lengththree$trajectory[a]),gsub(".+_","",names(lcomorbidities)[b]),sep="-")]]<-intersect(lengththreel[[lengththree$trajectory[a]]],lcomorbidities[[b]])
            }
          }
        }
        print(paste(round((a/length(lengththree$trajectory))*100,3),"%",sep=""))
      }
      save(lengthfourl,file="Data/Patients_by_4_steps_trajectories.RData")
      fin<-Sys.time()
      print(paste("Finished in:",fin-ini))
    }
    if(args[3]=="get_length_5"){
      ini<-Sys.time()
      ## Get 5-steps long trajectories ##
      ## @@ @@ @@ @@ @@ @@ @@ @@ @@ @@ ##
      ## Load the number of patients following a trajectory ##
      lengthfour<-fread("Data/Number_of_patients_in_4_steps_long_trajectories.txt",stringsAsFactors = F)
      ## Load the list indicating the patients following four steps long trajectories ##
      load("Data/Patients_by_4_steps_trajectories.RData") # lengthfourl
      ## Load the list indicating the patients suffering specific comorbidity relationships ##
      load("Data/Patients_by_significant_directional_comorbidities.RData") # lcomorbidities
      lengthfivel<-list()
      for(a in 1:length(lengthfour$V1)){
        # a<-1
        second<-gsub(".+-","",lengthfour$V1[a])
        cuales<-grep(paste(second,"_",sep=""),names(lcomorbidities))
        if(length(cuales)>0){
          for(b in cuales){
            # b<-4449
            cuantos<-length(intersect(lengthfourl[[lengthfour$V1[a]]],lcomorbidities[[b]]))
            if(cuantos>=20){
              recuento<-t(c(paste(gsub("-","_",lengthfour$V1[a]),gsub(".+_","",names(lcomorbidities)[b]),sep="-"),as.numeric(lengthfour[a,2:4]),cuantos))
              write.table(recuento,"Data/Number_of_patients_in_5_steps_long_trajectories.txt",quote=F,sep="\t",row.names=F,col.names = F,append = T)
              lengthfivel[[paste(gsub("-","_",lengthfour$V1[a]),gsub(".+_","",names(lcomorbidities)[b]),sep="-")]]<-intersect(lengthfourl[[lengthfour$V1[a]]],lcomorbidities[[b]])
            }
          }
        }
        print(paste(round((a/length(lengthfour$V1))*100,3),"%",sep=""))
      }
      save(lengthfivel,file="Data/Patients_by_5_steps_trajectories.RData")
      fin<-Sys.time()
      print(paste("Finished in:",fin-ini))
    }
    if(args[3]=="verify_trajectories"){
      ## Load the number of patients following a trajectory ##
      lengthfour<-fread("Data/Number_of_patients_in_4_steps_long_trajectories.txt",stringsAsFactors = F)
      ## Check diabetes trajectories ##
      lengthfour[grep("E11",lengthfour$V1)][which(lengthfour[grep("E11",lengthfour$V1)]$V4>=100)]
      copds<-lengthfour[grep("J44",lengthfour$V1)][which(lengthfour[grep("J44",lengthfour$V1)]$V4>=100)]
      ## Load the list indicating the patients following four steps long trajectories ##
      load("Data/Patients_by_4_steps_trajectories.RData") # lengthfourl
      thepats<-lengthfourl[["I10_E78_E11-E66"]]
      pats<-fread("Diagnoses_20080101_20181231_first_diagnoses_of_each_disease.txt",stringsAsFactors = F)
      setkey(pats,"idp")
      ordentab<-c()
      for(a in 1:length(thepats)){
        subpat<-pats[thepats[a]]
        ordentab<-rbind(ordentab,c(which(subpat$cod=="I10"),which(subpat$cod=="E78"),which(subpat$cod=="E11"),which(subpat$cod=="E66")))
      }
      which(ordentab[,3]>ordentab[,4])
    }
  }
}



if(args[1]=="NetSci_2025"){
  redes<-list.files("FinalNetworks/ByAge/")[grep("N0_5",list.files("FinalNetworks/ByAge/"))]
  
  tt<-fread("FinalNetworks/Global/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
  tt<-fread("FinalNetworks/AllGlobal/N0_5_Adjusted.txt",stringsAsFactors = F,sep="\t")
  library("igraph")
  age1<-fread("FinalNetworks/ByAge/N0_5_Adjusted_20-25.txt",stringsAsFactors = F,sep="\t")
  age2<-fread("FinalNetworks/ByAge/N0_5_Adjusted_40-45.txt",stringsAsFactors = F,sep="\t")
  age3<-fread("FinalNetworks/ByAge/N0_5_Adjusted_65-70.txt",stringsAsFactors = F,sep="\t")
  
  
  ## get 10 most frequent nodes
  ## get 10 most frequent categories
  
  sort(table(c(age1$Disease1,age1$Disease2)),decreasing = T)[1:20]
  sort(table(c(age2$Disease1,age2$Disease2)),decreasing = T)[1:20]
  sort(table(c(age3$Disease1,age3$Disease2)),decreasing = T)[1:20]
  
  gage1<-graph_from_data_frame(age1[,1:3])
  gage2<-graph_from_data_frame(age2[,1:3])
  gage3<-graph_from_data_frame(age3[,1:3])
  
  cage1<-fread("NetworksChronic/5years/ByAge/General_20-25.txt",stringsAsFactors = F,sep="\t")
  cage2<-fread("NetworksChronic/5years/ByAge/General_40-45.txt",stringsAsFactors = F,sep="\t")
  cage3<-fread("NetworksChronic/5years/ByAge/General_65-70.txt",stringsAsFactors = F,sep="\t")
  gcage1<-graph_from_data_frame(cage1[,1:2],directed = T)
  gcage2<-graph_from_data_frame(cage2[,1:2],directed = T)
  gcage3<-graph_from_data_frame(cage3[,1:2],directed = T)
  
  V(gcage1)$in_degree
  grados<-degree(gcage1)
  
  sort(table(c(cage1$Disease1,cage1$Disease2)),decreasing = T)[1:20]
  sort(table(c(cage2$Disease1,cage2$Disease2)),decreasing = T)[1:20]
  sort(table(c(cage3$Disease1,cage3$Disease2)),decreasing = T)[1:20]
  
  #### Chronic conditions ####
  edades<-list.files("NetworksChronic/5years/ByAge/")[grep("General",list.files("NetworksChronic/5years/ByAge/"))]
  edades<-edades[order(as.numeric(gsub("-.+","",gsub("General_","",edades))),decreasing = F)]
  lindegree<-list() ; loutdegree<-list()
  enfermedades<-c()
  comors<-c() ; diseases<-c()
  for(a in 1:length(edades)){
    # a<-1
    tt<-fread(paste("NetworksChronic/5years/ByAge/",edades[a],sep=""),stringsAsFactors = F,sep="\t")
    comors<-c(comors,length(tt$Disease1))
    diseases<-c(diseases,length(unique(c(tt$Disease1,tt$Disease2))))
    grafo<-graph_from_data_frame(tt[,1:2],directed=T)
    lindegree[[gsub(".txt","",gsub("General_","",edades))[a]]]<-degree(grafo,mode="in")
    loutdegree[[gsub(".txt","",gsub("General_","",edades))[a]]]<-degree(grafo,mode="out")
    enfermedades<-unique(c(enfermedades,tt$Disease1,tt$Disease2))
  }
  matin<-matrix(ncol = length(edades),nrow=length(enfermedades),0) ; colnames(matin)<-gsub(".txt","",gsub("General_","",edades)) ; rownames(matin)<-enfermedades
  matout<-matin
  matinall<-matin ; matoutall<-matin
  for(a in 1:length(names(lindegree))){
    # a<-1
    matin[names(lindegree[[a]]),a]<-as.numeric(lindegree[[a]])
    matout[names(loutdegree[[a]]),a]<-as.numeric(loutdegree[[a]])
  }
  difference<-t(matout-matin)
  difference<-difference[,order(colnames(difference),decreasing = F)]
  rownames(difference)[20]<-">95"
  pdf(file="HeatMap_comorbidities_age.pdf",height = 8,width = 16)
    heatmap.2(difference,key=F,
              # hclustfun = function(d) hclust(d, method = "ward.D2"),
              # distfun = function(d) dist(d,method = "binary"),
              col=colorpanel(10000,"#C81E17","#FFFFFF","#405191"),
              dendrogram = "none",Colv=FALSE,Rowv=FALSE,cexRow = 2,
              scale="none", margins=c(1,6), trace="none",
              lwid=c(0.3,4), lhei=c(0.1,4),colsep = c(3,26,33,55,93,111,123,124,163,175,183,188,215,231,232,244,245),sepcolor="black")
  dev.off()
  pdf(file="Barplots_comorbidities_age.pdf")
    barplot(comors,border=F)
    barplot(diseases,border=F)
  dev.off()
  maximos<-c() ; minimos<-c()
  for(a in 1:length(difference[1,])){
    # a<-1
    maximos<-c(maximos,length(which(difference[,a]==max(difference))))
    minimos<-c(minimos,length(which(difference[,a]==min(difference))))
  }
  
  colnames(difference)[which(maximos>=1)]
  colnames(difference)[which(minimos>=1)]
  difference[,which(maximos>=1)]
  difference[,which(minimos>=1)]
  
  ## Do it with all the diseases but plot only results for chronic ones ##
  #### All conditions ####
  edades<-list.files("FinalNetworks/ByAge/")[grep("N0_5_Adjusted",list.files("FinalNetworks/ByAge/"))]
  edades<-edades[order(as.numeric(gsub("-.+","",gsub("N0_5_Adjusted_","",edades))),decreasing = F)]
  lindegree<-list() ; loutdegree<-list()
  enfermedades<-c()
  for(a in 1:length(edades)){
    # a<-1
    tt<-fread(paste("FinalNetworks/ByAge/",edades[a],sep=""),stringsAsFactors = F,sep="\t")
    grafo<-graph_from_data_frame(tt[,1:2],directed=T)
    lindegree[[gsub(".txt","",gsub("N0_5_Adjusted_","",edades))[a]]]<-degree(grafo,mode="in")[intersect(names(degree(grafo,mode="in")),rownames(matinall))]
    loutdegree[[gsub(".txt","",gsub("N0_5_Adjusted_","",edades))[a]]]<-degree(grafo,mode="out")[intersect(names(degree(grafo,mode="out")),rownames(matinall))]
    enfermedades<-unique(c(enfermedades,tt$Disease1,tt$Disease2))
  }
  for(a in 1:length(names(lindegree))){
    # a<-1
    
    matinall[names(lindegree[[a]]),a]<-as.numeric(lindegree[[a]])
    matoutall[names(loutdegree[[a]]),a]<-as.numeric(loutdegree[[a]])
  }
  differenceall<-t(matoutall-matinall)
  differenceall<-differenceall[,order(colnames(differenceall),decreasing = F)]
  heatmap.2(differenceall,key=F,
            # hclustfun = function(d) hclust(d, method = "ward.D2"),
            # distfun = function(d) dist(d,method = "binary"),
            col=colorpanel(10000,"#C81E17","#FFFFFF","#405191"),
            dendrogram = "none",Colv=FALSE,Rowv=FALSE,cexRow = 0.6,cexCol = 0.6,
            scale="none", margins=c(1,5.2), trace="none",
            lwid=c(0.3,4), lhei=c(0.1,4))
}












