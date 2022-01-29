#Code for FELICITy:- Phenotype: PDQ

#####Part 1#####################################################################
options(stringsAsFactors = FALSE)
rm(list=ls())

#packages
library(openxlsx)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
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
library(Glimma)
library(RColorBrewer)

library(Homo.sapiens)
library(missMethyl)
library(DT)
library(Gviz)
#install.packages("Gviz")

#source("https://bioconductor.org/biocLite.R")
#biocLite("AnnotationHub")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("AnnotationHub")
#BiocManager::install("Glimma")
#BiocManager::install("bacon")

#library(bacon)

#library(DMRcate)
#library(Glimma)


##paths###

idatdir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/infos/ScanData/"
#list.files(idatdir)
phenodir <-  "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/FELICITy_data_anal_Jim/Updated/"
rdir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/Results_models_FELICITy/"
list.files(phenodir)

basedir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"


###function for fixing one id###
fixIt <- function(tofix) {
  tmp <- strsplit(tofix, "_")
  tmp <- lapply(tmp, function(x) {x[2] <- sprintf("%03d", as.numeric(x[2])); return(x)})
  out <- do.call(c, lapply(tmp, paste, collapse = "_"))
  out
}

samps <- read.xlsx(paste0(phenodir, "EPIC_plate_FELICITy_21090919.xlsx"))
samps$Patient.Code <- fixIt(samps$Patient.Code)
#samplesheet
targets <- read.metharray.sheet(idatdir, "M00936_Pl1_3_n114_Felicity.csv") 


targets$Sample_Name <- fixIt(targets$Sample_Name)
samps <- samps[match(targets$Sample_Name, samps$Patient.Code),]
##Phenotype data
phenos <- read.xlsx(paste0(phenodir, "Felicity_phenotypes_114.xlsx"))

phenos <- phenos[match(targets$Sample_Name, phenos$Code),]
phenos$Group <- factor(phenos$Group)

triFun <- function(x)
  cut(x, quantile(x, seq(0,1,1/3), na.rm = T), include.lowest = TRUE)

phenos$TriCortisol <- triFun(phenos$Cortisol)
phenos$TriBMI <- triFun(phenos$BMI_study_entry)
phenos$TriAge <- triFun(phenos$Age_at_inclusion)
phenos$TriGestAge <- triFun(phenos$Gestational_age_birth)
phenos$TriBW <- triFun(phenos$Birthweight_gm)
samps <- cbind(samps[,-c(11,13)], phenos[,-1])
samps$Stress.group <- factor(samps$Stress.group)
samps$Illumina.Plate <- factor(samps$Illumina.Plate)
samps$AChE.BChE <- samps$fet_AChE/samps$fet_BChE
samps$Stress.group <- relevel(samps$Stress.group, "2")
samps$Gender <- factor(samps$Gender)
samps$Smoking <- factor(samps$Smoking)

load(paste0(idatdir, "methdataall.Rdata"))
dat 

annoEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annoEPIC)

annoEPICSub <- annoEPIC[match(rownames(getM(eset.nosex)),annoEPIC$Name),                    
                        c(1:4,12:19,22:ncol(annoEPIC))]  

#######DMPs analysis######################################################################

sum(is.na(samps$PDQ))#0

sum(is.na(samps$Autoimmunediseases))#0
sum(is.na(samps$Gestational_Diabetes)) #0

dim(model.matrix(~ PDQ + Gender  + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth + Illumina.Plate, samps))
#114 , 8

design <- model.matrix(~ PDQ + Gender  + Smoking + Autoimmunediseases + Gestational_Diabetes + Gestational_age_birth +Illumina.Plate, samps)
sv <- sva(getM(eset.nosex), design, vfilter = 5e4)

design <- cbind(design, sv$sv)
fit <- lmFit(getM(eset.nosex), design)
fit2 <- eBayes(fit)

colnames(fit2)

toptable(fit2, 2)
#               logFC   AveExpr         t      P.Value  adj.P.Val        B
#cg06542869  0.02023160 -3.676819  5.776614 7.459816e-08 0.04945990 6.783856
#cg25217535 -0.01786849  3.333815 -5.484353 2.761658e-07 0.04945990 5.502898
#cg07990198 -0.02419701  3.821696 -5.462523 3.041145e-07 0.04945990 5.408697
#cg22861369  0.02892472 -2.731930  5.458921 3.089839e-07 0.04945990 5.393177
#cg14284689 -0.01330927  1.835948 -5.448675 3.232585e-07 0.04945990 5.349055
#cg05491000 -0.02547266  5.994000 -5.392569 4.136384e-07 0.04945990 5.108278
#cg01629131  0.03047279 -0.987858  5.384674 4.281956e-07 0.04945990 5.074510
#cg20587970 -0.05568583  1.566467 -5.326899 5.510777e-07 0.05569701 4.828284
#cg03105159  0.01597937  3.080270  5.223497 8.625066e-07 0.07748701 4.391486
#cg03127222  0.02160031 -4.511320  5.165803 1.105159e-06 0.08935807 4.149990



DMPs <- toptable(fit2, coef=2,confint = TRUE, genelist=annoEPICSub, number= 808554) 
bonfc <- 0.05/nrow(fit2)
bonfc
#[1] 6.183879e-08
sum(DMPs$P.Value < bonfc) #0
sum(DMPs$adj.P.Val < 0.05) #7
DMPs10 <- topTable(fit2, 2,confint = TRUE, genelist = annoEPICSub)


#write.table(DMPs10, "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/sample script from Rory_Norm/codes_2020_Jan/Model_codes/Stress_group/PDQtop10mod2_FELICITy_EPIC_DMPs.csv", 
#            sep=",", row.names=FALSE)

bval <- getBeta(eset.nosex)
bval[1:5, 1:5]

pdf(paste0(basedir,"cpg_pdq_adj.pdf"),onefile=T)
par(mfrow=c(2,2))
plot1 <- sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bval, cpg=cpg, pheno= samps$PDQ, ylab = "Beta values")
})
print(plot1)
dev.off()



##qqplot and lambda correction ##
basedir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"
setwd(basedir)
#install.packages("qqman")
library(qqman)
vignette("qqman")

##main model##
png(paste0(basedir,"qqplot_feli_pdqmod2.png"))
print(qq(DMPs$P.Value))
dev.off()

chisq = qchisq(DMPs$P.Value,1,lower.tail=FALSE)
lamb <- median(chisq, na.rm = T) / qchisq(0.5,1)
lamb
# 1.264856


##BACON correction
library(bacon)

tstas <- DMPs[c(42)] ##change to name

bc <- bacon(DMPs$t)
bc
#Bacon-object containing 1 set(s) of 808554 test-statistics.
#...estimated bias: -0.16.
#...estimated inflation: 1.1.

estimates(bc) 
#p.0          p.1          p.2       mu.0     mu.1      mu.2  sigma.0
#[1,] 0.9992606 2.074652e-05 0.0007186992 -0.1629567 2.972664 -2.978265 1.086775
#sigma.1   sigma.2
#[1,] 0.4971337 0.2516315

inflation(bc)
#sigma.0
#[1] 1.08681

bias(bc)
#mu.0
#[1] -0.1629976

head(pval(bc))
#                     t
#cg06542869 4.624436e-08
#cg25217535 9.765551e-07
#cg07990198 1.081390e-06
#cg22861369 2.305371e-07
#cg14284689 1.153414e-06
#cg05491000 1.495377e-06


head(tstat(bc))
#                   t
#cg06542869  5.465178
#cg25217535 -4.896305
#cg07990198 -4.876218
#cg22861369  5.172861
#cg14284689 -4.863477
#cg05491000 -4.811853

#histogram#
png(paste0(basedir,"bc_pdq_hist.png"))
plot1 <- plot(bc, type="hist")
print(plot1)
dev.off()

#qqplot##
png(paste0(basedir,"qqplot_bc2.png")) ##bc plot
plot(bc, type="qq")
dev.off()

png(paste0(basedir,"qqplot_bc_pdqmod2.png")) ##p-values plot
print(qq(pval(bc)))
dev.off()

chisq = qchisq(pval(bc),1,lower.tail=FALSE)
lamb <- median(chisq, na.rm = T) / qchisq(0.5,1)
lamb
#[1] 1.047653

BC_pvals <- pval(bc)

DMPs1 <- (cbind(DMPs, BC_pvals)) 
dim(DMPs)
#[1] 808554     44
 dim(DMPs1)
#[1] 808554     45

head(DMPs1)
#               P.Value adj.P.Val        B    corrpvals
#cg06542869 7.459816e-08 0.0494599 6.783856 4.618781e-08
#cg25217535 2.761658e-07 0.0494599 5.502898 9.746568e-07
#cg07990198 3.041145e-07 0.0494599 5.408697 1.079302e-06
#cg22861369 3.089839e-07 0.0494599 5.393177 2.302902e-07
#cg14284689 3.232585e-07 0.0494599 5.349055 1.151196e-06
#cg05491000 4.136384e-07 0.0494599 5.108278 1.492553e-06

bonfc <- 0.05/nrow(DMPs1)
bonfc # 6.183879e-08

DMPs1 <- mutate(DMPs1,fdr = p.adjust(BC_pvals,method="fdr"), bonfp = pmin(BC_pvals*nrow(DMPs1),1))
sum(DMPs1$bonfp < 0.05,na.rm=T) #[1] 1
sum(DMPs1$fdr   < 0.05,na.rm=T) #[1] 1

write.table(DMPs1, "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/sample script from Rory_Norm/codes_2020_Jan/Model_codes/Stress_group/PDQcorrecp-vals_mod2_FELICITy_EPIC_DMPs.csv", 
            sep=",", row.names=FALSE)

save(DMPs1, samps,  file = paste0(rdir,"DMPS_pdq_adjustedmod.Rdata"))
DMPs1 <- DMPs1[order(DMPs1$ BC_pvals), ]


#####Manhattan plot: PDQ################################################
##Manhattan Plot#####
#package
library(ggplot2)
library(readxl)
library(limma)
library(dplyr)load(paste0(rdir, "DMPS_pdq_adjustedmod.Rdata"))

 ##loading the results datarame after it has been run by glm
dim(DMPs)
#[1] 808554      45
rownames(DMPs1) <- DMPs1$Name

resultsdf <- DMPs1 ##check the rownames carefully

resultsdf[1:6, 1:6]
#               chr       pos strand       Name Probe_rs Probe_maf
#cg06542869 chr11 102028101      + cg06542869       <NA>        NA
#cg22861369  chr5 131608620      + cg22861369       <NA>        NA
#cg01629131 chr20  10787894      + cg01629131       <NA>        NA
#cg03105159  chr2    291062      + cg03105159       <NA>        NA
#cg03127222  chr7   7221920      + cg03127222   rs758263  0.347315
#cg25217535 chr10  35027638      - cg25217535 rs11009905  0.012587

ppath <- "/mnt/nas/global/ame/users/Ritika/Projects/KORA_Age_Project/KORA_AGE_data/Full_data/Data/"

anno <- read.csv(paste0(ppath,"MethylationEPIC_v-1-0_B4.csv"),skip=7,fill=T , stringsAsFactors=F) 

anno <- anno[anno$IlmnID %in% rownames(resultsdf), c("IlmnID","CHR","MAPINFO", "UCSC_RefGene_Name")]
dim(anno)
#[1] 808554      4


#anno$CHR<-as.numeric(anno$CHR)
rownames(anno) <- anno$IlmnID

head(anno)
#CHR   MAPINFO       UCSC_RefGene_Name
#cg07881041  19   5236016 PTPRS;PTPRS;PTPRS;PTPRS
#cg18478105  20  61847650                  YTHDF1
#cg23229610   1   6841125
#cg03513874   2 198303466
#cg09835024  NA  24072640                  EIF2S3
#cg05451842  14  93581139       ITPK1;ITPK1;ITPK1


anno <- anno[rownames(resultsdf),]
colnames(anno)
#[1] "CHR"               "MAPINFO"           "UCSC_RefGene_Name"

max(-log10(resultsdf$P),na.rm = T) #checking the limits of y-axis, i.e the maximum value for P.
#[1]  7.447564

lin <- 0.05/nrow(resultsdf) ##setting the threshold for ewas ##should be -log10
#[1]  6.183879e-08

allannores <- cbind(anno,resultsdf) 
#rm(anno,resultsdf)

allannores[1:10, 1:10]
#               IlmnID CHR   MAPINFO                        UCSC_RefGene_Name
#cg11409463 cg11409463   5  77634319
#cg05306225 cg05306225   8   4848947                              CSMD1
#cg20905655 cg20905655  19  18520786
#cg22725303 cg22725303  19  36427617                              LRFN3
#cg26938014 cg26938014   5 122432628                              PRDM6
#cg20646907 cg20646907   5 154799021
#cg22700843 cg22700843   6 169053865                        SMOC2;SMOC2
#cg16992373 cg16992373   6 125117600 NKAIN2;NKAIN2;NKAIN2;NKAIN2;NKAIN2
#cg01097027 cg01097027  22  30598622
#cg01495737 cg01495737   1 163037823                     RGS4;RGS4;RGS4
#                 logFC    AveExpr         t      P.Value   adj.P.Val         B
#cg11409463  0.002862145  2.8029989  6.312115 7.076325e-09 0.005721591 6.5900298
#cg05306225 -0.002842984 -4.3913055 -5.574817 2.008067e-07 0.054931589 3.2991356
#cg20905655  0.004060955  2.6358929  5.571445 2.038142e-07 0.054931589 3.2845611
#cg22725303 -0.002183348 -4.7556816 -5.379382 4.717749e-07 0.095363880 2.4625221
#cg26938014 -0.002542892  3.1972674 -5.138596 1.322840e-06 0.191301192 1.4550844
#cg20646907  0.001512035  0.9293826  5.093733 1.598633e-06 0.191301192 1.2703670
#cg22700843 -0.002203494  4.1399653 -5.072794 1.745831e-06 0.191301192 1.1844841
#cg16992373  0.002045072  3.1500214  5.053543 1.892773e-06 0.191301192 1.1057089
#cg01097027 -0.001696773 -2.8284404 -4.997879 2.388756e-06 0.203959253 0.8789575
#cg01495737 -0.001962601  3.5159481 -4.939668 3.042410e-06 0.203959253 0.6434723

a<- resultsdf[1:10, 1:47]
write.table(a, "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/sample script from Rory_Norm/codes_2020_Jan/Model_codes/Stress_group/PDQcorrecp-vals_mod2_FELICITy_EPIC_DMPs10.csv", 
            sep=",", row.names=FALSE)

#ordering the dataframe according first to chromosome and then by base position
allannores <- allannores[order(allannores$CHR, allannores$MAPINFO), ] 
#%>% rename( c("Gene" = "UCSC_RefGene_Name", "BP" = "MAPINFO"))


sapply(allannores, class)

unique(allannores$CHR)

as.data.frame(table(allannores$CHR))
allannores <- (allannores[,-c(17,18)])

##---------------------------------------------------------------------------------------------------------------------------------------
#packages and paths again

basedir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
##--------------------------------------------------------------------------------------------
fac <- c("CHR", "IlmnID")

allannores_clean <- allannores
allannores_clean [fac] <- lapply(allannores[fac], factor)

glimpse(allannores_clean)

max(-log10(resultsdf$fdr),na.rm = T) #checking the limits of y-axis, i.e the maximum value for P.
#[1]  7.447564


sig = 6.183879e-08 # significant threshold line
sugg = 8.079453e-08  ##suggestive threshold limit
#------------------
mah_plot <- function(allannores, threshold){
  # format allannores
  plot_data <- allannores %>% 
    
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len= as.numeric(max(MAPINFO))) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(allannores_clean, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each IlmnID
    arrange(CHR, MAPINFO) %>%
    mutate( BPcum = as.numeric(MAPINFO +tot))  #
  
  # Add highlight and annotation information
  # mutate( is_highlight=ifelse(IlmnID %in% hlight, "yes", "no")) %>%
  # mutate( is_annotate=ifelse(P < threshold, "yes", "no"))
  
  # get chromosome center positions for x-axis
  axisdf <- plot_data %>% 
    group_by(CHR) %>% 
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  plot <- ggplot(plot_data, aes(x=BPcum, y=-log10(P.Value))) +
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(c("#4EAFAF", "#2C9696", "#0F8F8F", "#057272", "#005A5A"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR,  breaks= axisdf$center) + 
    scale_y_continuous(expand = c(0, 0)) + # expand=c(0,0)removes space between plot area and x axis 
    
    
    ylim(0,20) +
    theme_light() +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.line = element_line(color = "black")) +
    xlab("Chromosome") + 
    # add x label
    geom_label_repel( data=plot_data %>% filter(P.Value < sugg), # add annotation value
                      aes(label=IlmnID), size=3) + # add annotation
    geom_point(data= plot_data %>% filter(P.Value < threshold), # add annotation value
               color="orange", size=2) + # Add highlighted points 
    geom_hline(yintercept = -log10(threshold)) +
    geom_hline(yintercept = -log10(sugg), linetype="dashed")# threshold line
  
  return(plot) # return the final plot
}

png(paste0(basedir,"feli_pdq_mah.png")) #, width=1425, height=975
plot <- mah_plot(allannores_clean, threshold = 6.28402e-08) # run function
print(plot)
dev.off()


##DNA Methylated Regions for PDQ################################################################
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
#maxbetafc   meanbetafc overlapping.promoters#
#1 -0.002068872 -0.001869703                  <NA>
#  2 -0.005450478 -0.003694398    DAXX-010, DAXX-008


write.table(dmr.table, "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/sample script from Rory_Norm/codes_2020_Jan/Model_codes/Stress_group/PDQ_DMR1_FELi.csv", 
            sep=",", row.names=FALSE)

pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(samps$PDQ))]
#names(groups) <- levels(factor(samps$Group))


# draw the plot for the top DMR
basedir <- "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/EPIC data/"
setwd(basedir)

phen.col <- c(rep("orange", 38), rep("blue", 38))

pdf(paste0(basedir,"dmrpdq1_basic.pdf"), onefile=T)   
par(mfrow=c(1,1))
plot <- DMR.plot(ranges = results.ranges, dmr = 1, CpGs = getBeta(eset.nosex), phen.col=phen.col, 
         what = "Beta", arraytype = "EPIC", genome = "hg19")
print(plot)
dev.off()

##gviz plot#####

head(DMRs$results)
#          coord no.cpgs       minfdr  Stouffer    maxbetafc
#1 chr6:33288180-33288372       7 1.181304e-10 0.5602015 -0.005513139
#    meanbetafc
#1 -0.003712802

##indicate which genome is being used
#setting up the genomic region 
gen= "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 1

# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

cpgislanddir <- "/mnt/nas/global/ame/users/Ritika/Study in R/Data_methylation_paper/"

##cpg islands
islandHMM <- read.csv(paste(cpgislanddir, "Cpg_islands_chr6",sep="/"),
                      sep="\t", stringsAsFactors=FALSE, header=TRUE)

head(islandHMM)

islandData <- GRanges(seqnames=Rle(islandHMM[ ,'chrom']),
                      ranges=IRanges(start=islandHMM[,'chromStart'], end=islandHMM[,'chromEnd']),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData

# DNAseI hypersensitive sites
dnase <- read.csv(paste(cpgislanddir,"Dnashs_regions",
                        sep="/"), sep="\t",stringsAsFactors=FALSE,header=FALSE)
head(dnase)

dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])
dnaseData


#Now, set up the ideogram, genome and RefSeq tracks that will provide context for our methylation data.

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")

rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq",
                    from=minbase, to=maxbase, trackType="GeneRegionTrack",
                    rstarts="exonStarts", rends="exonEnds", gene="name",
                    symbol="name2", transcript="name", strand="strand",
                    fill="darkblue",stacking="squish", name="RefSeq",
                    showId=TRUE, geneSymbol=TRUE)


#Ensure that the methylation data is ordered by chromosome and base position.

annEPICkOrd <- annoEPIC[order(annoEPIC$chr,annoEPIC$pos),]
bvalOrd <- bval[match(annEPICkOrd$Name,rownames(bval)),]

#Create the data tracks:
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(annEPICkOrd$chr),
                   ranges=IRanges(start=annEPICkOrd$pos, end=annEPICkOrd$pos),
                   strand=Rle(rep("*",nrow(annEPICkOrd))),
                   betas=bvalOrd)

# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])


pal <- brewer.pal(8,"Dark2")
# methylation data track
methTrack <- DataTrack(range=cpgData, 
                       groups=samps$PDQ, 
                       genome = gen,
                       chromosome=chrom,
                       ylim=c(-0.05,1.05),
                       col=pal,
                       type=c("a","p"), 
                       name="DNA Meth.\n(beta value)",
                       background.panel="white", 
                       legend=TRUE, 
                       cex.title=0.8,
                       cex.axis=0.8, 
                       cex.legend=0.8)



# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")


tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks

pdf(paste0(basedir, filename = "dmrpdq.pdf") #, width = 15, height = 10, units = "in", res = 400))
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))

dev.off()


pdf(paste0(basedir,  "dmrpdq.pdf"))
    plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
               add35=TRUE, sizes = sizes,  length(tracks))
    
    dev.off()
    
    
#####extracting the cpgs present in the DMR#####################

RR <- as.data.frame(results.ranges)
RR$DMRID <-rownames(RR)
row.names(RR) = NULL
RR$DMRNO <- rownames (RR)
row.names(RR) = RR$DMRID
RR$DMRID = NULL
RR <-RR[order(RR$minfdr), , drop = FALSE]
RR

#Now you need to pull the CpG info from the dmr output that was used to make the results ranges
cgID <- as.data.frame(DMRs$input)

#Look at your RR file, choose a DMR that looks good, and take note of the DMRNO; Run
DMRNUM <- readline(prompt = "What is your DMR Number:") ##fill 1

#Enter the number into the console and hit enter, then run these lines 
#and it should spit out a table listing the probes as well as other useful info

assign(paste0("DMR_",DMRNUM), subset(subset(RR,DMRNO==DMRNUM)))

assign(paste0("DMR_",DMRNUM,"_probelist"), 
      subset(cgID, cgID$CHR==assign(paste0("DMR_",DMRNUM), 
      subset(subset(RR,DMRNO==DMRNUM)))$seqnames & cgID$pos>assign(paste0("DMR_",DMRNUM), 
        subset(subset(RR,DMRNO==DMRNUM)))$start-1 & cgID$pos<assign(paste0("DMR_",DMRNUM),
         subset(subset(RR,DMRNO==DMRNUM)))$end+1))

DMR_1_probelist
                ID  weights  CHR      pos        betafc    indfdr is.sig
#670766 cg03477252 1.734360 chr6 33288180 -0.0008179342 0.7299469  FALSE
#670767 cg22904406 2.548896 chr6 33288296 -0.0037361208 0.5764809  FALSE
#670768 cg26500914 2.970843 chr6 33288323 -0.0030850470 0.4901230  FALSE
#670769 cg09365002 2.764500 chr6 33288329 -0.0049693718 0.5335826  FALSE
#670770 cg09597022 3.412261 chr6 33288332 -0.0055131393 0.4189778  FALSE
#670771 cg07905975 3.749222 chr6 33288366 -0.0042468376 0.3719578  FALSE
#670772 cg24498636 2.807470 chr6 33288372 -0.0036211619 0.5266657  FALSE
#                raw          fdr  sig
#670766 2.805744e-15 3.240851e-10 TRUE
#670767 8.343604e-16 1.181304e-10 TRUE
#670768 7.887452e-16 1.181304e-10 TRUE
#670769 7.877290e-16 1.181304e-10 TRUE
#670770 7.884092e-16 1.181304e-10 TRUE
#670771 8.532850e-16 1.181304e-10 TRUE
#670772 8.766047e-16 1.181304e-10 TRUE

DMR_cpgs <- merge(DMR_1_probelist, DMPs1, by.x = "ID", by.y = "Name") 

write.table(DMR_cpgs, "/mnt/nas/global/ame/users/Ritika/Projects/FELICITy_Project/EPIC_DATA_norm/sample script from Rory_Norm/codes_2020_Jan/Model_codes/Stress_group/PDQ_DMRcpgs_FELi.csv", 
            sep=",", row.names=FALSE)
       
#####p cutoff:- 0.1###


DMRs1 <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = 0.1)

head(DMRs1$results)
results.ranges1 <- extractRanges(DMRs1, genome = "hg19")


#########################################pathway analysis##########################################

resultspdq <-  load(paste0(rdir, "DMPS_pdq_adjustedmod.Rdata"))

DMPs1 <- DMPs1[order(DMPs1$ BC_pvals), ] ###arrange according to corrected pvals

##top 100 cpgs
sigCpGs <- DMPs1$Name[1:100]
sigCpGs[1:10]
length(sigCpGs) #100

all <- DMPs1$Name
length(all)
#808554


png(paste0(basedir,"gometh_pdq.png"))
par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE, array.type = c("EPIC"))
print(gst)
dev.off()


topGSA(gsa, number=10)
#GO:0009986                                                                            cell surface
#GO:0060025                                                         regulation of synaptic activity
#GO:0030054                                                                           cell junction
#GO:0019842                                                                         vitamin binding
#GO:0003779                                                                           actin binding
#GO:0004494                                                       methylmalonyl-CoA mutase activity
#GO:1904990          regulation of adenylate cyclase-inhibiting dopamine receptor signaling pathway
#GO:1904992 positive regulation of adenylate cyclase-inhibiting dopamine receptor signaling pathway
#GO:0019481                                          L-alanine catabolic process, by transamination
#GO:0047305                           (R)-3-amino-2-methylpropionate-pyruvate transaminase activity
#Ont    N DE         P.DE FDR
#GO:0009986  CC  756 12 0.0001287742   1
#GO:0060025  BP    3  2 0.0001663392   1
#GO:0030054  CC 1213 15 0.0019655416   1
#GO:0019842  MF  134  4 0.0022674207   1
#GO:0003779  MF  400  8 0.0022886451   1
#GO:0004494  MF    1  1 0.0027469596   1
#GO:1904990  BP    1  1 0.0028101151   1
#GO:1904992  BP    1  1 0.0028101151   1
#GO:0019481  BP    1  1 0.0028852743   1
#GO:0047305  MF    1  1 0.0028852743   1

