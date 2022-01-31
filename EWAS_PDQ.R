#Association between Methylation and PDQ: EWAS
#Analyst: Ritika Sharma
#Updated: 31.01.2022

#####################################
options(stringsAsFactors = FALSE)
rm(list=ls())

##Packages################
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
library(sva)
library(limma)
library(RColorBrewer)
library(DMRcate)
library(Homo.sapiens)
library(missMethyl)

##define paths###################

idatdir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/infos/ScanData/"
rdir <- "N:/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Results_models_FELICITy/" ##results 
basedir <- "N:/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"  #basedirectory

##loading the phenotype data######################
samps <-  read.csv(paste0(basedir,"Samps_FELICITy.csv"),fill=T , stringsAsFactors=F) 

##loading normalised methylation data
load(paste0(idatdir, "methdataall.Rdata"))
dat 

annoEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annoEPIC)

#matching cpg site ronames of mvals and annoepic
annoEPICSub <- annoEPIC[match(rownames(getM(eset.nosex)),annoEPIC$Name),                    
                        c(1:4,12:19,22:ncol(annoEPIC))]                       


################DMPs:-PDQ (adjusted model)##########################################################################################################
##checking variables##
sum(is.na(samps$PDQ))#0

##Creating a design model matrix
dim(model.matrix(~ PDQ + Gender  + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth + Illumina.Plate, samps))

design <- model.matrix(~ PDQ + Gender  + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth +Illumina.Plate, samps)

##Generating surrogate variables
sv <- sva(getM(eset.nosex)[,!is.na(samps$Cortisol)], design,vfilter = 5e4) #3

##combining the deisgn matrix and surrogate variables
design <- cbind(design, sv$sv)

##Finally fitting in the limma model
fit <- lmFit(getM(eset.nosex), design)
fit2 <- eBayes(fit)

dim(fit2)
#[1] 808554      11

colnames(fit2)

topTable(fit2, 2)

DMPs <- topTable(fit2, coef=2,confint = TRUE, genelist=annoEPICSub, number= 808554) 

bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #0
sum(DMPs$adj.P.Val < 0.05) #7

##saving the results
write.table(DMPs,paste0(rdir, "/PDQ/PDQ_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)

####Lambda correction######

##current lambda##
chisq = qchisq(DMPs$P.Value,1,lower.tail=FALSE)
lamb <- median(chisq, na.rm = T) / qchisq(0.5,1)
lamb
# 1.264856  #high!

##BACON correction##
library(bacon)

bc <- bacon(DMPs$t)
bc

estimates(bc) 

inflation(bc)
#sigma.0
#[1] 1.08681

bias(bc)
#mu.0
#[1] -0.1629976

head(pval(bc))

head(tstat(bc))

#qqplot##
png(paste0(basedir,"qqplot_bc2.png")) ##bc plot
plot(bc, type="qq")
dev.off()

chisq = qchisq(pval(bc),1,lower.tail=FALSE)
lamb <- median(chisq, na.rm = T) / qchisq(0.5,1)
lamb
#[1] 1.047653 ##new!

BC_pvals <- pval(bc)

DMPs1 <- (cbind(DMPs, BC_pvals)) 
dim(DMPs)

DMPs1 <- mutate(DMPs1,fdr = p.adjust(BC_pvals,method="fdr"), bonfp = pmin(BC_pvals*nrow(DMPs1),1))
sum(DMPs1$bonfp < 0.05,na.rm=T) #[1] 1
sum(DMPs1$fdr   < 0.05,na.rm=T) #[1] 1

DMPs1 <- DMPs1[order(DMPs1$ BC_pvals), ]

##saving the results
write.table(DMPs,paste0(rdir, "/PDQ/PDQ_corrected_FELICITy_DMPs.csv"), sep=",", row.names=FALSE)

save(DMPs1, samps,  file = paste0(rdir,"DMPS_PDQ_adjustedmod.Rdata"))

##DMRs:- DNA Methylated Regions for PDQ################################################################
library(DMRcate)

myAnnotation <- cpg.annotate(object = getM(eset.nosex), datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = FALSE, 
                             coef = 2, 
                             arraytype = "EPIC")

#Your contrast returned 1 individually significant probes; a small but real effect. 
#Consider manually setting the value of pcutoff to return more DMRs, but be warned that doing this increases 
#the risk of Type I errors.

str(myAnnotation)
#List of 7
#$ ID    : chr [1:808554] "cg14817997" "cg26928153" "cg16269199" "cg13869341" ...
#$ stat  : num [1:808554] 0.456 -0.75 0.437 -0.461 0.651 ...
#$ CHR   : chr [1:808554] "chr1" "chr1" "chr1" "chr1" ...
#$ pos   : int [1:808554] 10525 10848 10850 15865 18827 29407 29425 68849 68889 69591 ...
#$ betafc: num [1:808554] 0.00021 -0.000104 0.000222 -0.000228 0.000316 ...
#$ indfdr: num [1:808554] 0.944 0.898 0.947 0.944 0.914 ...
#$ is.sig: logi [1:808554] FALSE FALSE FALSE FALSE FALSE FALSE ...
#- attr(*, "row.names")= int [1:808554] 1 2 3 4 5 6 7 8 9 10 ...
#- attr(*, "class")= chr "annot"

##pcutoff= 0.05
# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")


results.ranges
#GRanges object with 1 range and 6 metadata columns:
#   seqnames            ranges strand |   no.cpgs              minfdr
#         <Rle>         <IRanges>  <Rle> | <integer>           <numeric>
#  [1]     chr6 33288180-33288372      * |         7 1.1813037202563e-10
#              Stouffer            maxbetafc           meanbetafc
#              <numeric>            <numeric>            <numeric>
#  [1] 0.560201486993845 -0.00551313927791463 -0.00371280179127199
#       overlapping.promoters
#                <character>
#  [1]    DAXX-010, DAXX-008
#-------
#  seqinfo: 1 sequence from an unspecified genome; no seqlengths

# visualization
dmr.table <- data.frame(results.ranges)

#seqnames    start      end width strand no.cpgs       minfdr  Stouffer
#1    chr17 30243828 30243860    33      *       2 7.350857e-09 0.3683492
#2     chr6 33288180 33288372   193      *       7 7.350857e-09 0.7347085
#     maxbetafc   meanbetafc overlapping.promoters#
#1 -0.002068872 -0.001869703                  <NA>
#  2 -0.005450478 -0.003694398    DAXX-010, DAXX-008



write.table(dmr.table, paste0(rdir, "/PDQ/PDQ_DMR1_FELi.csv"), sep=",", row.names=FALSE)
##############################################################

