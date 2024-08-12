#R script for the whole process of INTERVAL factor analysis and fine mapping paper (public version)
#Version: date20240807

## Background
#The *Script_INTERVAL.R:* describes the whole process of the main analysis in more detailed steps and provides the R codes for running 
#factor analysis and fine mapping:

## Overview of the whole process
# - Step_1: collect and combine the 99 (raw) blood cell traits in the INTERVAL cohort of UK blood donors [^1][^2][^3];
# - Step_2: delete the rows/individuals of sample with missing values (so the sample size is reduced to 18k);
# - Step_3: the reduced-size 99 blood cell traits are normalised using inverse-normal rank transformation;
# - Step_4: run factor analysis (FA) based on reduced-size 99 normalised blood cell traits; 
# - Step_5: get the optimal 25 FA latent factors with a matrix contains all loadings of 99 blood cell traits;
# - Step_6:  use inverse normal rank transformation for each of the 25 latent factors;
# - Step_7: link all 99 blood cell traits with the 25 FA latent factors and create a network/connection visualization;
# - Step_8: re-use BOLT-LMM GWAS of full-size 99 blood cell traits [^2][^3];
# - Step_9: run BOLT-LMM to get GWAS of reduced-size 99 blood cell traits and 25 latent factors; 
# - Step_10: compare GWAS signals between different blood cell traits and latent factors;
# - Step_11: check the (combination of) strengths and regions of GWAS signals;
# - Step_12: check correlations of genotypes between GWAS significant SNPs;
# - Step_13: implement conditional analyses by adding lead SNPs ([^2][^3]) into BOLT-LMM regressions;
# - Step_14: detect novel SNPs based on FA in comparison with blood cell traits (both full-size 43k and the reduced-size 18k samples);
# - Step_15: use VEP (build 37) and connections of latent factors with blood cell traits to better understand novelty;
# - Step_16: select regions for fine-mapping based on both FA latent factors and blood cell traits, compare the results.

## Reference
# [^1]: https://www.donorhealth-btru.nihr.ac.uk/studies/interval-study/
# [^2]: Astle, W.J., et al. (2016). The Allelic Landscape of Human Blood Cell Trait Variation and Links to Common Complex Disease. Cell 167, 1415–1429.e19.
# [^3]: Akbari, P., et al. (2023). A genome-wide association study of blood cell morphology identifies cellular proteins implicated in disease aetiology. Nat. Commun. 14, 5023.


##Stage_1: Factor Analysis and Inverse Normal Transformation----
#install related R packages
library(RNOmni)  #RNOmni: Rank Normal Transformation Omnibus Test (Inverse normal transformation (INT) based genetic association testing.)
library(psych)   #Implement factor analysis

# Step_1: collect and combine the 99 blood cell traits in the INTERVAL cohort of UK blood donors;
data_raw_full <- read.csv("Important_INTERVAL_63_36_traits_combined.csv") #the raw data cannot be shared publicly
data_raw_full = data_raw_full[-1,]
data_raw_99trait_names <- read.csv("Important_trait_names_99.csv")
data_raw_99trait_names = paste0(data_raw_99trait_names$Trait_name,'_gwas_normalised')
data_raw_99trait <- data_raw_full[,c('FID',data_raw_99trait_names)]
data_raw_99trait = as.data.frame(data_raw_99trait)
data_raw_99trait_all = data_raw_99trait
for (i in 2:dim(data_raw_99trait)[2]){
  data_raw_99trait_check = as.numeric(data_raw_99trait[,i])
  data_raw_99trait_all[,i] = data_raw_99trait_check
}
colnames(data_raw_99trait_all) = gsub("_gwas_normalised", "", colnames(data_raw_99trait_all))

# Step_2: delete the rows/individuals of sample with missing values (so the sample size is reduced to 18k);
data_raw_99trait_all_18k = data_raw_99trait_all[complete.cases(data_raw_99trait_all), ]

# Step_3: use inverse normal rank transformation on each of the 99 blood cell traits;
data_raw_99trait_all_18k_normalised = data_raw_99trait_all_18k
for(i in 2:100){
  data_raw_99trait_all_18k_normalised[,i] <- RankNorm(data_raw_99trait_all_18k[,i])
}

# Step_4: run factor analysis (FA) with varimax rotation on reduced-size 99 normal-transformed blood cell traits; 
# Step_5: get the optimal 25 FA latent factors with a matrix contains all loadings of 99 blood cell traits;
# Step_6: (optional) the 25 FA latent factors are also normalised but the effect/impact is small;
#library(psych)
yk = data_raw_99trait_all_18k_normalised[,-1]
#prop.fn <- function(x) abs(x)/sum(abs(x)) #check
prop.fn <- function(x) x^2/sum(x^2)
yz <- apply(yk,2,function(x) (x-mean(x,na.rm=T))/sqrt(var(x,na.rm=T)))
fa.parallel(yz,fm="ml")       #fm="ml" will do a maximum likelihood factor analysis; scree plot suggests selecting 25 latent factors
faY <- fa(r=yz, nfactors=25, fm="ml",rotate="varimax")
fYstar <- faY$scores
fYstar_df <- as.data.frame(fYstar)
fYstar_df$FID = data_raw_99trait_all_18k_normalised$FID
#create raw-latent 99-by-25 factor loading table based on prop.fn <- function(x) x^2/sum(x^2)
Pyload <- t(apply(faY$loadings,1,prop.fn))
#heatmap(Pyload, Rowv = NA, Colv = NA) 
Pyload_FA25 = as.data.frame(Pyload)
Pyload_FA25$cell_type = 0
raw_trait_names <- read.csv("Important_trait_names_99.csv")
for (i in 1:length(Pyload_FA25$ML1)){Pyload_FA25$cell_type[i] = raw_trait_names$Cell_type[raw_trait_names$Trait_name==rownames(Pyload_FA25)[i]]}

# Step_7: link all 99 blood cell traits with the 25 FA latent factors and create a network/connection visualization;
# The following dataset combines 25 FA latent factors and 99 blood cell traits, as well as PCs and the Clinic column for BOLT-LMM
Data2018_add_fa_18310_all99traits_normalised <- read.delim("Data2018_add_fa_18310_all99traits_normalised.txt")
# Columns are: FID, IID, missing, PC.1-10, clinic, 25 FA latent factors and 99 blood cell traits (18,310 individual ids)


##Stage_2: BOLT-LMM and Conditional Analyses----
# The following steps were implemented in hpc with BOLT-LMM software
# Step_8: re-use BOLT-LMM GWAS of full-size 99 blood cell traits [^2][^3];
# Step_9: run BOLT-LMM to get GWAS of reduced-size 99 blood cell traits and 25 latent factors; 
# Step_10: compare GWAS signals between different blood cell traits and latent factors;
# Step_11: check the (combination of) strengths and regions of GWAS signals;
# Step_12: check correlations of genotypes between GWAS significant SNPs;
# Step_13: implement conditional analyses by adding lead SNPs ([^2][^3]) into BOLT-LMM regressions;
# Step_14: detect novel SNPs based on FA in comparison with blood cell traits (both full-size 43k and the reduced-size 18k samples);
# Step_15: use VEP (build 37) and connections of latent factors with blood cell traits to better understand novelty;


##Stage_3: Fine-mapping----
# Step_16: select regions for fine-mapping based on both FA latent factors and blood cell traits, compare the results.
# sub_steps as follows: GET genotypes from INTERVAL bgen files
# sub_step_1: We have big size impute_20_interval.bgen (i.e. 5.52GB) and impute_20_interval.bgen.bgi files from INTERVAL project
# sub_step_2: Use "qctool" to get the regional subset small file, e.g. reg239_date0127.bgen
#(reference: https://www.chg.ox.ac.uk/~gav/qctool_v2/documentation/examples/filtering_variants.html)
#Terminal code: ./qctool -g impute_20_interval.bgen -s interval.samples -og reg239_date0127.bgen -incl-range 20:3875857-4422772 -ofiletype bgen_v1.2 -bgen-bits 8
#Note: in order to use snp_readBGEN() in R package "bigsnpr", we need put: "-ofiletype bgen_v1.2 -bgen-bits 8" in the qctool code
#(reference: https://privefl.github.io/bigsnpr/reference/snp_readBGEN.html)
# sub_step_3: Use the "BGEN" package (i.e. the bgenix tool) to get the index of reg239_date0127.bgen (i.e. reg239_date0127.bgen.bgi)
#(reference: https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md)
# sub_step_4: Now we have a smaller bgen file and its index file: reg239_date0127.bgen and reg239_date0127.bgen.bgi
library(bigsnpr)
library(vroom)
library(propagate)
library(RSQLite)
#library(R2BGLiMS)  #can be ignored if using library(flashfmZero)
#library(MGflashfm) #can be ignored if using library(flashfmZero)
#source("makegroups1trait_v20240222.R") #can be ignored if using library(flashfmZero)
#source("No-Cor-flashfm-v20240222.R")   #can be ignored if using library(flashfmZero)
library(flashfmZero) 

#testing region239 in chr20 as an example###
#setwd("region239_testing")
#get GWAS
load('reg239_chr20_gwas_corx_raf_flipcheck.RData') #cleaned GWAS inputs

#Or load all regions
# args=(commandArgs(TRUE))
# reg=as.numeric(args[1])  # region number
# fm_item5b_updated_csv <- read.csv("fm_item5b_updated_csv_sorted.csv")
# sel_chr <- fm_item5b_updated_csv$CHR[reg]
# sel_reg <- fm_item5b_updated_csv$Region_id[reg]
# sel_start <- fm_item5b_updated_csv$start_updated[reg]
# sel_end <- fm_item5b_updated_csv$end_updated[reg]
# 
# new_dir = paste0("./tmpDIR/regid_",sel_reg)
# if(!dir.exists(file.path(new_dir))) dir.create(file.path(new_dir))
# 
# ##get GWAS  
# load(paste0('reg',sel_reg,'_chr',sel_chr,'_gwas_corx_raf_flipcheck_231225.RData'))


#Note: save both reg239.bgen and reg239.bgen.bgi in the current working directory
bgenfiles = 'reg239.bgen' 
backingfile = 'reg239_output' 
list_snp_id=list()
list_snp_id[[1]] = paste0(GWAS_list_fa25$ML13$CHR,'_',
                          GWAS_list_fa25$ML13$BP,'_',
                          GWAS_list_fa25$ML13$ALLELE1,'_',
                          GWAS_list_fa25$ML13$ALLELE0)

rds = snp_readBGEN(
  bgenfiles,
  backingfile,
  list_snp_id,
  ind_row = NULL,
  bgi_dir = dirname(bgenfiles),
  read_as = "dosage",
  ncores = 1
)

chr=20; start=lb=3875857; end=ub=4422772 # region 239 as an example

########

test <- snp_attach(rds)

head(GWAS_list_fa25[[1]])
head(test$map)

### NOTEs: 
## 1. in Ensembl.B37 rs7509501 has freq 0.838 for G in EUR (and freq 0.159 for A allele)
## 2. So, in GWAS_list_fa25$ML13 the rs7509501 allele counted for A1FREQ is ALLELE1 because A1FREQ=0.847 and G = ALLELE1
## 3. and in test$map rs7509501  has freq=0.153, so must be counting A, which is allele2 here
## 4. So, we need snpinfo$ALLELE1 to match snpinfo2$allele2  

snpinfo2 = test$map[,c(2,4,5,6)]
snpinfo2 = as.data.frame(snpinfo2)
names(snpinfo2)[1:4] = c("rs","ps","allele1","allele2")

genotypes = test$genotypes[1:dim(test$genotypes)[1],1:dim(test$genotypes)[2]]
colnames(genotypes) = paste0("chr",chr,"_",snpinfo2$ps)
genotypes[1:20,1:5]

interval_sample_id <- read.csv("interval.samples", sep="")
rownames(genotypes) = as.character(interval_sample_id[2:length(interval_sample_id$ID_1),1])
#our 18k a subset of 43k (here we need the individual IIDs, which are saved in Data2018_add_fa_18310_all99traits_normalised)
genotypes = genotypes[c(as.character(Data2018_add_fa_18310_all99traits_normalised$IID)),]
genotypes[1:20,1:5]

rownames(snpinfo2) = paste0("chr",chr,"_",snpinfo2$ps)
rsnp <- which(snpinfo2$ps>=lb & snpinfo2$ps<=ub)
genotypes <- genotypes[,rsnp]
snpinfo2 <- snpinfo2[rsnp,]

# align alleles in RP and GWAS at intersecting snps
# do intersection in first trait of GWAS_list_raw99, since all traits have same variants
snpinfo <- GWAS_list_raw99[[1]]
rownames(snpinfo) <- paste0("chr",chr,"_",snpinfo$BP)
#snpkeep <- which(snpinfo$INFO > 0.8) # keep only INFO > 0.8
#snpinfo <- snpinfo[snpkeep,]
isnp <- intersect(rownames(snpinfo), colnames(genotypes) )
X <- genotypes[,isnp]
snpinfo2 <- snpinfo2[isnp,]
snpinfo <- snpinfo[isnp,]

#flip in raw99 GWAS and fa25 GWAS
gwas.list.interval.raw99 = list()
gwas.list.interval.fa25 = list()
flip1 <- which(snpinfo$ALLELE1 != snpinfo2$allele2)
if(length(flip1)>0) check <- which(snpinfo$ALLELE1[flip1] == snpinfo2$allele1[flip1] & snpinfo$ALLELE0[flip1] == snpinfo2$allele2[flip1]) 
if(length(check)>0){
  for(i in 1:length(GWAS_list_raw99)){
    #beta
    rownames(GWAS_list_raw99[[i]]) <- paste0("chr",GWAS_list_raw99[[i]]$CHR,"_",GWAS_list_raw99[[i]]$BP)
    GWAS_list_raw99[[i]] <- GWAS_list_raw99[[i]][rownames(snpinfo),]
    beta_sel = GWAS_list_raw99[[i]]$BETA
    names(beta_sel) = GWAS_list_raw99[[i]]$snp_id
    beta_sel[flip1[check]] = -beta_sel[flip1[check]]
    #eaf
    raf_sel <- GWAS_list_raw99[[i]]$A1FREQ
    names(raf_sel) = GWAS_list_raw99[[i]]$snp_id
    raf_sel[flip1[check]] = 1-raf_sel[flip1[check]]
    #gwas.list
    gwas.list.interval.raw99[[names(GWAS_list_raw99)[i]]]$rsID = names(beta_sel)
    gwas.list.interval.raw99[[names(GWAS_list_raw99)[i]]]$beta = unname(beta_sel)
    gwas.list.interval.raw99[[names(GWAS_list_raw99)[i]]]$EAF = unname(raf_sel) 
    #as.data.frame
    gwas.list.interval.raw99[[names(GWAS_list_raw99)[i]]] = as.data.frame(gwas.list.interval.raw99[[names(GWAS_list_raw99)[i]]])
  }
  for(i in 1:length(GWAS_list_fa25)){
    #beta
    rownames(GWAS_list_fa25[[i]]) <- paste0("chr",GWAS_list_fa25[[i]]$CHR,"_",GWAS_list_fa25[[i]]$BP)
    GWAS_list_fa25[[i]] <- GWAS_list_fa25[[i]][rownames(snpinfo),]
    beta_sel = GWAS_list_fa25[[i]]$BETA
    names(beta_sel) = GWAS_list_fa25[[i]]$snp_id
    beta_sel[flip1[check]] = -beta_sel[flip1[check]]
    #eaf
    raf_sel <- GWAS_list_fa25[[i]]$A1FREQ
    names(raf_sel) = GWAS_list_fa25[[i]]$snp_id
    raf_sel[flip1[check]] = 1-raf_sel[flip1[check]]
    #gwas.list
    gwas.list.interval.fa25[[names(GWAS_list_fa25)[i]]]$rsID = names(beta_sel)
    gwas.list.interval.fa25[[names(GWAS_list_fa25)[i]]]$beta = unname(beta_sel)
    gwas.list.interval.fa25[[names(GWAS_list_fa25)[i]]]$EAF = unname(raf_sel) 
    #as.data.frame
    gwas.list.interval.fa25[[names(GWAS_list_fa25)[i]]] = as.data.frame(gwas.list.interval.fa25[[names(GWAS_list_fa25)[i]]])
  }
}

qc <- function(g,theta=0.1,BestGuess=TRUE) {
  ind0 <- which(g<=theta)
  ind1 <- which(g>=1-theta & g<=1+theta)
  ind2 <- which(g>=2-theta)
  bg <- rep(NA,length(g))
  if(BestGuess) {
    if(length(ind0)>0) bg[ind0] <- 0
    if(length(ind1)>0) bg[ind1] <- 1
    if(length(ind2)>0) bg[ind2] <- 2 
  } else{
    if(length(ind0)>0) bg[ind0] <- g[ind0]
    if(length(ind1)>0) bg[ind1] <- g[ind1]
    if(length(ind2)>0) bg[ind2] <- g[ind2]   
  }  
  return(bg)
}

Xqc <- apply(X,2,qc,theta=0.2,BestGuess=TRUE)
#> dim(Xqc) reg 239 theta=0.1
#[1] 18310  2629

Nprop <- apply(Xqc,2,function(x) sum(!is.na(x)))/nrow(Xqc)
keep <- which(Nprop >= 0.80)

Xqc <- Xqc[,keep]

for(i in 1:length(gwas.list.interval.raw99)) gwas.list.interval.raw99[[i]] <- gwas.list.interval.raw99[[i]][keep,]
for(i in 1:length(gwas.list.interval.fa25)) gwas.list.interval.fa25[[i]] <- gwas.list.interval.fa25 [[i]][keep,] 

ldout <- bigcor(Xqc, size=min(2000, ncol(Xqc)), verbose=FALSE, fun= "cor", use="p" )
corX <- as.matrix(ldout[1:dim(ldout)[1],1:dim(ldout)[1]])
rownames(corX) <- colnames(corX) <- gwas.list.interval.raw99[[1]]$rsID

rafqc <- apply(Xqc,2,mean,na.rm=T)/2
names(rafqc) <- gwas.list.interval.raw99[[1]]$rsID

summary(abs(rafqc-gwas.list.interval.raw99[[1]]$EAF))

# X, snpinfo, snpinfo2 all have same variant ordering, based in isnp, so can use names from snpinfo
#raf <- apply(X,2,mean)/2
#names(raf) <- snpinfo$snp_id
#corX <- cor(X,use="p")
## sanity check for alignment of alleles in gwas and rp:
#summary(abs(raf-gwas.list.interval.raw99[[1]]$EAF))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 1.174e-07 2.458e-07 2.742e-07 4.042e-07 1.398e-06 




#FLASHFM0withJAMd - use for multiple latent factors
save.path="./results/tmpDIR"

#delete traits with no signals
remaining_fa25 = names(gwas.list.interval.fa25)
remaining_raw99 = names(gwas.list.interval.raw99)

id_nosignal = c()
for(j in 1:length(gwas.list.interval.fa25)){
  if(min(gwas.list.interval.fa25[[j]]$P_value)>1e-6) id_nosignal = c(j, id_nosignal)
}
#print("check id_nosignal is ", id_nosignal)
if(is.null(id_nosignal)){
  gwas.list.interval.fa25 <- gwas.list.interval.fa25
}else{
  gwas.list.interval.fa25 <- gwas.list.interval.fa25[-id_nosignal]
  remaining_fa25 <- remaining_fa25[-id_nosignal]
}

id_nosignal = c()
for(j in 1:length(gwas.list.interval.raw99)){
  if(min(gwas.list.interval.raw99[[j]]$P_value)>1e-6) id_nosignal = c(j, id_nosignal)
}
#print(paste0("check id_nosignal is ", id_nosignal))
if(is.null(id_nosignal)){
  gwas.list.interval.raw99 <- gwas.list.interval.raw99
}else{
  gwas.list.interval.raw99 <- gwas.list.interval.raw99[-id_nosignal]
  remaining_raw99 <- remaining_raw99[-id_nosignal]
}



# use this function for one latent factor or one blood cell trait
if(length(gwas.list.interval.fa25)==1){
  # for one latent factor 
  JAMdwithGroups.out.fa25.single <- JAMdwithGroups(gwas.list.interval.fa25[[1]], N=18310, corX, cred = 0.99,
                                                   jam.nM.iter =5, maxcv = 1, maxcv_stop = 20, 
                                                   min.mppi = 0.01, r2.minmerge = 0.8)
}else{
  JAMdwithGroups.out.fa25.single <- 0
}

if(length(gwas.list.interval.raw99)==1){
  #  for one blood cell trait
  JAMdwithGroups.out.raw99.single <- JAMdwithGroups(gwas.list.interval.raw99[[1]], N=18310, corX, cred = 0.99,
                                                    jam.nM.iter =5, maxcv = 1, maxcv_stop = 20, 
                                                    min.mppi = 0.01, r2.minmerge = 0.8)
}else{
  JAMdwithGroups.out.raw99.single <- 0
}


# JAM on multiple blood cell traits 
if(length(gwas.list.interval.raw99)>1){
  multiJAMd.out.raw99 <- multiJAMd(gwas.list.interval.raw99,  corX, N=18310, save.path,
                                   maxcv = 1, maxcv_stop = 20, jam.nM.iter =5, r2.minmerge=0.8, minsnpmppi = 0.01,
                                   #NCORES=length(gwas.list.interval.raw99) #if running on the hpc
                                   NCORES=1) 
  multiJAMd.CSraw99  <- multiJAMdCS(multiJAMd.out.raw99, cred = 0.99)
}else{
  multiJAMd.out.raw99 <- 0
  multiJAMd.CSraw99 <- 0
}




# multiple latent factors flashfmZero
if(length(gwas.list.interval.fa25)>1){
  FLASHFM0withJAMd.mt0g1.fa25 <- FLASHFMZEROwithJAMd(gwas.list.interval.fa25, 
                                                     corX, 
                                                     N = 18310, 
                                                     save.path, 
                                                     TOdds = 1,
                                                     cpp = 0.99, 
                                                     #NCORES=length(gwas.list.interval.fa25), #if running on the hpc
                                                     NCORES = 1,
                                                     maxcv = 1, 
                                                     maxcv_stop = 20, 
                                                     jam.nM.iter = 5, 
                                                     r2.minmerge = 0.8, 
                                                     minsnpmppi = 0.01)
  FLASHFM0withJAMd.mtg1CS.fa25 <- allcredsetsPP(FLASHFM0withJAMd.mt0g1.fa25$mpp.pp,cred=.99)
}else{
  FLASHFM0withJAMd.mt0g1.fa25 <- 0
  FLASHFM0withJAMd.mtg1CS.fa25 <- 0
}

save(FLASHFM0withJAMd.mt0g1.fa25, 
     FLASHFM0withJAMd.mtg1CS.fa25, 
     multiJAMd.out.raw99, 
     multiJAMd.CSraw99, 
     JAMdwithGroups.out.fa25.single,
     JAMdwithGroups.out.raw99.single,
     gwas.list.interval.fa25,
     gwas.list.interval.raw99,
     corX,
     remaining_fa25,
     remaining_raw99,
     GWAS_list_fa25,
     GWAS_list_raw99,
     Data2018_IID,
     fm_item5b_updated_csv,
     interval_sample_id,
     sel_reg,
     sel_chr,
     sel_start,
     sel_end,
     file = paste0("reg",id_reg,"_chr",id_chr,"_fm_results_20240223.RData"))


##END----




