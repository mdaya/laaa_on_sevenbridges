args <- commandArgs(trailingOnly = TRUE)

allele.dose.file <-  args[1]
afr.dose.file <-  args[2]
allele.afr.dose.file <- args[3]
pheno.file <- args[4]
r.model.file <- args[5]
min.maf <- args[6]
out.file.name <- args[7]

#Load the model functions
source(r.model.file)

#Read the dose frames and the phenotype frame
allele.frame <- read.delim(allele.dose.file)
afr.frame <- read.delim(afr.dose.file)
allele.afr.frame <- read.delim(allele.afr.dose.file)
pheno.frame <- read.delim(pheno.file)

#Create the output file
cat("position\tref\talt\talt_frq\tn\tunadj_allele_dose_beta\tunadj_allele_dose_beta_se\tunadj_allele_dose_p\tallele_dose_beta\tallele_dose_beta_se\tallele_dose_p\tafr_dose_beta\tafr_dose_beta_se\tafr_dose_p\tallele_afr_dose_beta\tallele_afr_dose_beta_se\tallele_afr_dose_p\tanova_p\n", sep="", file=out.file.name, append=F)

for (position in allele.frame$position) {
  #Merge in the allele dose
  allele.col.frame <- 
    data.frame(sample_id=names(allele.frame)[-c(1:3)],
               allele_dose=t(allele.frame[allele.frame$position == position,-c(1:3)]))
  names(allele.col.frame)[2] <- "allele_dose"
  model.frame <- merge(pheno.frame, allele.col.frame)
  
  #Merge in the afr dose
  afr.col.frame <- 
    data.frame(sample_id=names(afr.frame)[-c(1:3)],
               afr_dose=t(afr.frame[afr.frame$position == position,-c(1:3)]))
  names(afr.col.frame)[2] <- "afr_dose"
  model.frame <- merge(model.frame, afr.col.frame)
  
  #Merge in the allele.afr dose
  allele.afr.col.frame <- 
    data.frame(sample_id=names(allele.afr.frame)[-c(1:3)],
               allele.afr_dose=t(allele.afr.frame[allele.afr.frame$position == position,-c(1:3)]))
  names(allele.afr.col.frame)[2] <- "allele_afr_dose"
  model.frame <- merge(model.frame, allele.afr.col.frame)
  
  #Fit the model
  model <- runLaaaModelSummary(model.frame)
  m.null <- runNullModel(model.frame)
  m.laaa <- runLaaaModel(model.frame)
  anova_p <- anova(m.null, m.laaa)[2,6]
  allele.model <- runAlleleModelSummary(model.frame)
  
  #Get the ref and alt alleles
  ref <- allele.frame$ref[allele.frame$position == position]
  alt <- allele.frame$alt[allele.frame$position == position]
  
  #Get P-value colname (this could be t or z depending on which model is fit, so determine this dynamically)
  p.val.col <- colnames(model$coefficients)[grep("^Pr", colnames(model$coefficients))]
  
  #Estimate N and the alternate allele frequency
  n <- nrow(model.frame)
  frq <- sum(model.frame$allele_dose)/(n*2)
  
  if (frq >= min.maf) {
    
    #Get the model output
    if ("allele_dose" %in% rownames(model$coefficients)) {
      allele_dose_beta <- model$coefficients["allele_dose", "Estimate"]
      allele_dose_beta_se <- model$coefficients["allele_dose", "Std. Error"]
      allele_dose_p <- model$coefficients["allele_dose", p.val.col]      
    } else {
      allele_dose_beta <- NA
      allele_dose_beta_se <- NA
      allele_dose_p <- NA  
    }
    if ("afr_dose" %in% rownames(model$coefficients)) {
      afr_dose_beta <- model$coefficients["afr_dose", "Estimate"]
      afr_dose_beta_se <- model$coefficients["afr_dose", "Std. Error"]
      afr_dose_p <- model$coefficients["afr_dose", p.val.col]      
    } else {
      afr_dose_beta <- NA
      afr_dose_beta_se <- NA
      afr_dose_p <- NA  
    }
    if ("allele_afr_dose" %in% rownames(model$coefficients)) {
      allele_afr_dose_beta <- model$coefficients["allele_afr_dose", "Estimate"]
      allele_afr_dose_beta_se <- model$coefficients["allele_afr_dose", "Std. Error"]
      allele_afr_dose_p <- model$coefficients["allele_afr_dose", p.val.col]      
    } else {
      allele_afr_dose_beta <- NA
      allele_afr_dose_beta_se <- NA
      allele_afr_dose_p <- NA  
    }
    if ("allele_dose" %in% rownames(allele.model$coefficients)) {
      unadj_allele_dose_beta <- allele.model$coefficients["allele_dose", "Estimate"]
      unadj_allele_dose_beta_se <- allele.model$coefficients["allele_dose", "Std. Error"]
      unadj_allele_dose_p <- allele.model$coefficients["allele_dose", p.val.col]      
    } else {
      unadj_allele_dose_beta <- NA
      unadj_allele_dose_beta_se <- NA
      unadj_allele_dose_p <- NA  
    }
    
    #Write the output
    cat(position,"\t", ref, "\t", alt, "\t", frq, "\t", n, "\t",
        unadj_allele_dose_beta,"\t",
        unadj_allele_dose_beta_se,"\t",
        unadj_allele_dose_p,"\t",
        allele_dose_beta,"\t",
        allele_dose_beta_se,"\t",
        allele_dose_p,"\t",
        afr_dose_beta,"\t",
        afr_dose_beta_se,"\t",
        afr_dose_p,"\t",
        allele_afr_dose_beta,"\t",
        allele_afr_dose_beta_se,"\t",
        allele_afr_dose_p, "\t",
        anova_p, "\n", file=out.file.name, append=T, sep="")
  }

}


