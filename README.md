A docker environment and calling scripts to run a R version of the LAAA method PMID: 29226381 on Biodata catalyst

The docker build can be pulled from https://quay.io/repository/mdaya/laaa

To create the allele dose, African dose, and allele-African dose frames required
to run LAAA, run the following command from the docker instance:

```
bash /home/analyst/create_dose_frames.sh \
   alleles_rephased_file \
   vit_file \
   snp_info_file \
   sample_id_file \
   chr \
   begin_hg19_pos \
   end_hg19_pos
```

To run LAAA, run the following command from the docker instance:

```
bash /home/analyst/run_laaa.sh \
   allele_dose_file \
   afr_dose_file \
   allele_afr_dose_file \
   phenotype_file \
   r_model_file \
   min_maf \
   out_file_name
```

To create locus zoom plots of the output, run the following command from the docker instance:

```
bash /home/analyst/locus_zoom.sh \
   assoc_file \
   begin_hg19_pos \
   end_hg19_pos \
   out_file_prefix
```

## create_dose_frames parameters

All parameters are mandatory and should be specified in order

### alleles\_rephased\_file

Full path file name of the rephased allele file output by RFMix

### vit\_file

Full path file name of the Viterbi file output by RFMix

### snp\_info\_file

Full path file name of a file containing the (ordered) positions of the variants in the
alleles\_rephased\_file as the first column of the file. 

### sample\_id\_file

Full path file name of a file containing the sample IDs in the alleles\_rephased\_file (in
the same order as they appear in the alleles\_rephased\_file ) 

### chr

The number of the chroomsome from which dosages should be extracted

### begin\_hg19\_pos

The hg19 base pair position that marks the beginning of the genomic region to be
extracted

### end\_hg19\_pos

The hg19 base pair position that marks the end of the genomic region to be
extracted

## run\_laaa parameters

All parameters are mandatory and should be specified in order

### allele\_dose\_file 

The allele\_dose\_file output by create\_dose\_frames

### afr\_dose\_file 

The afr\_dose\_file output by create\_dose\_frames

###  allele\_afr\_dose\_file 

The allele\_afr\_dose\_file output by create\_dose\_frames

### phenotype\_file

A tab delimitted phenotype file that contains the outcome variable and all
covariates that should be included in the LAAAA model. The first column should
be the sample ID. This sample ID will be mapped to the IDs in the sample\_id
file. 

### r\_model\_file

A .R file containing the following functions, which will be used to fit the LAAA
model. The Y and covariates list should be updated as required for the specific
analysis, and variables should be named exactly according to the column names in
the phenotype\_file

```
runLaaaModelSummary <- function(model.frame) {
  return (summary(lm(Y ~ covariates + allele_dose + afr_dose + allele_afr_dose, data=model.frame)))
}

runNullModel <- function(model.frame) {
  return (lm(Y ~ covariates, data=model.frame))
}

runLaaaModel <- function(model.frame) {
  return (lm(Y ~ covariates + allele_dose + afr_dose + allele_afr_dose, data=model.frame))
}

runAlleleModelSummary <- function(model.frame) {
  return (summary(lm(Y ~ covariates + allele_dose, data=model.frame)))
}
```

### min\_maf

The minimum MAF a variant should have to be included in the model

### out\_file\_name

The output file name.

### Content of output file

The output file will contain the follwing columns:

* position: Base pair position of variant
* ref: Reference allele
* alt: Alternate allele (the beta coefficients reflect the association between the number of copies of this allele and the outcome)
* alt\_frq: Frequency of the alternate allele
* unadj\_allele\_dose\_beta: Beta coefficient of allelic dose not considering local ancestry
* unadj\_allele\_dose\_beta\_se: Standard error of unadj\_allele\_dose\_beta
* unadj\_allele\_dose\_beta\_p: P-value of unadj\_allele\_dose
* allele\_dose\_beta: Beta coefficient of allelic dose considering local ancestry (LAAA model)
* allele\_dose\_beta\_se: Standard error of allele\_dose\_beta
* allele\_dose\_beta\_p: P-value of allele\_dose
* afr\_dose\_beta: Beta coefficient of local ancestry expressed in terms of number of copies of African ancestry ("African" dose; LAAA model)
* afr\_dose\_beta\_se: Standard error of afr\_dose\_beta
* afr\_dose\_beta\_p: P-value of afr\_dose
* allele\_afr\_dose\_beta: Beta coefficient of the combination of allelic dose and African local ancestry (LAAA model)
* allele\_afr\_dose\_beta\_se: Standard error of allele\_afr\_dose\_beta
* allele\_afr\_dose\_beta\_p: P-value of allele\_afr\_dose
* anova\_p: P-value of the model improvement when considering allele\_dose, afr\_dose, allele\_afr\_dose compared to the model with only covariates (LAAA model p-value)

## locus\_zoom parameters

### assoc\_file

The full path name of the association output file produced by run\_laaa

### begin\_hg19\_pos

The hg19 base pair position that marks the beginning of the genomic region to
plot

### end\_hg19\_pos

The hg19 base pair position that marks the end of the genomic region to plot

### out\_file\_prefix

Prefix to be given to the locus zoom PDF output files
