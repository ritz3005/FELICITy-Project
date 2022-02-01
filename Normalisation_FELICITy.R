#Normalisation of Methylation data
#Analyst: James Macdonald 
#         Ritika Sharma
#Updated: 31.01.2022

#####################################
options(stringsAsFactors = FALSE)
rm(list=ls())

##Loading required packages for analysis################
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
library(parallel)
library(readxl)
library(affycoretools)
library(BiocStyle)
library(xtable)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
data(SNPs.147CommonSingle)
library(Homo.sapiens)
library(AnnotationHub)

##paths###

idatdir <-  "N:/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/infos/ScanData/"
phenodir <-  "N:/users/Ritika/Projects/FELICITy_Project/FELICITy_data_anal_Jim/Updated/"
rdir <-  "N:/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"
list.files(phenodir)


###function for fixing one id###
fixIt <- function(tofix) {
  tmp <- strsplit(tofix, "_")
  tmp <- lapply(tmp, function(x) {x[2] <- sprintf("%03d", as.numeric(x[2])); return(x)})
  out <- do.call(c, lapply(tmp, paste, collapse = "_"))
  out
}

samps <- read.xlsx(paste0(phenodir, "EPIC_plate_FELICITy_21090919.xlsx"))
samps$Patient.Code <- fixIt(samps$Patient.Code)

#read in the sample sheet for the experiment
targets <- read.metharray.sheet(idatdir, "M00936_Pl1_3_n114_Felicity.csv") 
targets$Sample_Name <- fixIt(targets$Sample_Name)

samps <- samps[match(targets$Sample_Name, samps$Patient.Code),]

##Phenotype data
phenos <- read.xlsx(paste0(phenodir, "Felicity_phenotypes_114.xlsx"))
phenos <- phenos[match(targets$Sample_Name, phenos$Code),]
phenos$Group <- factor(phenos$Group)


samps <- cbind(samps[,-c(11,13)], phenos[,-1])
samps$Stress.group <- factor(samps$Stress.group)
samps$Illumina.Plate <- factor(samps$Illumina.Plate)
samps$AChE.BChE <- samps$fet_AChE/samps$fet_BChE
samps$Stress.group <- relevel(samps$Stress.group, "2")
samps$Gender <- factor(samps$Gender)
samps$Smoking <- factor(samps$Smoking)

##read in the raw data from the IDAT files
RGset <- read.metharray.exp(targets = targets)

#quality control and normalisation of the data

if(!file.exists("methdataall.Rdata")){
  dat <- RGset
  dp <- detectionP(dat)
  failed <- dp > 0.01
  cpgs <- row.names(dp)[rowSums(failed) >= 5]
  snps <- row.names(SNPs.147CommonSingle)[apply(SNPs.147CommonSingle[,c(4,6)], 1, function(x) any(!is.na(x)))]
  eset <- preprocessFunnorm(dat, ratioConvert = FALSE)
  eset <- eset[!row.names(eset) %in% c(snps, cpgs),]
  infind <- is.infinite(rowSums(getM(eset)))
  eset <- eset[!infind,]
  set.seed(0xabeef)
  eset.nosex <- subsetByOverlaps(eset, keepSeqlevels(rowRanges(eset), paste0("chr", 1:22), pruning.mode = "coarse"))
  save(dat, eset, eset.nosex,  file = paste0(idatdir,"methdataall.Rdata"))
} else {
  load("methdataall.Rdata")
}

dat 
#########################################################################################


