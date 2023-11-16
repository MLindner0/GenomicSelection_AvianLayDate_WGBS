
### Detecting multilocus adaptation using Redundancy Analysis (RDA)
# following: https://popgen.nescent.org/2018-03-27_RDA_GEA.html

options(width=200)

# Load packages
# library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(adegenet) # Used to load data
library(dplyr)
library(ggplot2)
library(cowplot)

# read snps data:
snps.0 <- read.PLINK("plink/individuals_.9_QC_ok_RDA/ERC_wgbs.raw")
snps <- read.PLINK("plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs.raw")

# get matrix with individuals (rows) and snps (columns) and remove NAs:
snp.mat <- as.data.frame(t(as.matrix(snps)))
snp.mat.0 <- snp.mat[complete.cases(snp.mat),]
snp.dat <- t(as.matrix(snp.mat.0))

# get phenotype data:
SampleInfo <- read.table("src/SampleData.txt", sep="\t", header=T)
SampleData <- SampleInfo[SampleInfo$Pool==1 & SampleInfo$Run==1,c(2,4,5,8)]
names(SampleData)[1] <- "sample.id"
Phen_dat <- SampleData[match(snps@ind.names, SampleData$sample.id),]
pred <- data.frame(Line_factor=as.factor(Phen_dat$Line_numeric))

# run RAD
snps.rda <- rda(snp.dat ~ Line_factor, data=pred, scale=T)
snps.rda
save(snps.rda, file="temp/trial_snps.rda.RData")

RsquareAdj(snps.rda)
summary(eigenvals(snps.rda, model = "constrained"))

signif.full <- anova.cca(snps.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
write.table(as.data.frame(signif.full), "out/snps.signif.full", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# same here as only one predictor
signif.axis <- anova.cca(snps.rda, by="axis", parallel=getOption("mc.cores"))
write.table(as.data.frame(signif.axis), "out/snps.signif.axis", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

# check variance inflation factors (irrelevant for one explanatory variable)
vif.cca(snps.rda)

# extract the SNP loadings:
load.rda <- scores(snps.rda, choices=c(1), display="species")

# fix SNP identifier 
help <- str_split_fixed(rownames(load.rda), "_", 3)
help.2 <- paste(help[,1], help[,2], sep="_")
rownames(load.rda) <- help.2

# make plot 
hist(load.rda[,1], main="Loadings on RDA1")

load.rda.0 <- as.data.frame(load.rda)
col <- "darkslategrey"
up <- mean(load.rda.0$RDA1)+3*sd(load.rda.0$RDA1)
low <- mean(load.rda.0$RDA1)-3*sd(load.rda.0$RDA1)

Plot <- load.rda.0  %>% 
  ggplot() +
  geom_histogram(aes(x = RDA1), position="identity", binwidth=0.02, fill=col, color="white") +  
  ylab("Frequency") + xlab("Load on RAD1") +
  theme_classic() +
  geom_vline(xintercept=up, linetype=2, size=0.5, color="darkorange") +
  geom_vline(xintercept=low, linetype=2, size=0.5, color="darkorange") +
  geom_vline(xintercept=mean(load.rda.0$RDA1), linetype=1, size=0.5, color="darkorange") +
  scale_x_continuous(breaks=seq(-0.09, 0.09, length.out=7), expand=c(0,0)) +
  scale_y_continuous(breaks=seq(0, 30000, length.out=4), labels=c(paste(seq(0, 30, length.out=4), c("", rep("k",3)), sep="")), limits=c(0,40000), expand=c(0,0)) +
  theme(axis.text=element_text(size=8),
        axis.title.y=element_text(size=9, margin=margin(t = 0, r = 9, b = 0, l = 0)),
        axis.title.x=element_text(size=9, margin=margin(t = 9, r = 0, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin=margin(t = 10, r = 10, b = 4, l = 4)) 
save_plot("plots/RDA.SNPs.pdf", Plot, base_height = 3, base_width = 3.5)

# function to find outlier SNPs
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand.snps <- outliers(load.rda[,1],3)
cand.snp.names <- names(cand.snps)

# overlapp with bayescan?
load(file="out/Fst_all.RData")
snps.bayes <- Fst_all[Fst_all$bayeS_qval<0.05,]$snp
cand.snp.names[cand.snp.names %in% snps.bayes] 

# snps in data set after filter?
snps.0 <- colnames(snps.rda$Ybar)
snps.0[snps.0 %in% snps.bayes] 

# mark "significant" snps:
rda.load.snps <- as.data.frame(load.rda)
rda.load.snps$outlier <- ifelse(row.names(rda.load.snps) %in% cand.snp.names, 1,0)
save(rda.load.snps, file="final/rda.load.snps.RData")

# prepare table for manuscript
snp_annotation <- read.table("temp/annotated_Snps", sep="\t")
RDA.genes <- snp_annotation[snp_annotation$V4 %in% cand.snp.names,]
RDA.genes$V9 <- gsub("LOC107203824", "SOX3", RDA.genes$V9)
write.table(unique(RDA.genes$V9), "out/RDA.genes.txt", quote=F, col.names=F, row.names=F, sep="\n")

load.rda.out <- as.data.frame(load.rda); load.rda.out$V4 <- row.names(load.rda.out)
load.rda.out0 <- load.rda.out[row.names(load.rda.out) %in% cand.snp.names,]
RDA.out <- merge(load.rda.out0, snp_annotation[,c(4,11,9)], by="V4", all.x=TRUE)
RDA.out$V9 <- gsub("LOC107203824", "LOC107203824 (SOX3)", RDA.out$V9)

# fix chr name!
Chr_length_names <- read.table("../wgbs_snakemake_bismark/genome/chr_length_names", header=F, sep="\t") 
names(Chr_length_names) <- c("Chr_name", "Chr_name_short", "chr", "length")

RDA.out$chr <- str_split_fixed(RDA.out$V4, ":", 2)[,1]
help <- str_split_fixed(RDA.out$V4, ":", 2)[,2]
help.0 <- str_split_fixed(help, ":", 2)

RDA.out$help <- paste(help.0[,1], help.0[,2], sep="_")

RDA.out0 <- merge(RDA.out, Chr_length_names, by="chr")
RDA.out0 <- RDA.out0[order(abs(RDA.out0$RDA1), decreasing=TRUE),]
RDA.out0$SNP <- paste(RDA.out0$Chr_name, RDA.out0$help, sep="_")

RDA.out1 <- RDA.out0[,c(10,3:5)]
names(RDA.out1) <- c("SNP", "SNP loading", "Genomic region", "Gene symbol")
write.table(RDA.out1, "out/RDA.mainTable.txt", quote=F, col.names=TRUE, row.names=F, sep="\t")