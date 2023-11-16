### --------------------------------------------------------------------------------
###
### SNPs from wgbs data: filter
### Project: WGBS ERC selection lines
###
### Author: M.Lindner (M.Lindner@nioo.knaw.nl), v202107

options(width=200)
library(GenABEL)

library(ggplot2)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(dplyr)
library(stringr)



### --- Individuals ------------------------------------------------------

snps_raw <- read.table("plink/individuals_raw/ERC_wgbs_trans.tped")
samples_raw <- read.table("plink/individuals_raw/ERC_wgbs_trans.tfam")

# set up loop to format input data
sample_names <- samples_raw$V2
helper <- seq(5, 43, by=2)
snps_base <- data.frame(snp_id=snps_raw$V2)

for(i in 1:length(sample_names)){
  i2 <- helper[i]
  snps_base$new <- snps_raw[,i2]
  names(snps_base)[1+i] <- sample_names[i]
}
row.names(snps_base) <- snps_base$snp_id
snps_base$snp_id <- NULL

EarlyBird <- snps_base[,grep("F3_E", names(snps_base))]
EB_temp <- apply(EarlyBird, 1, function(x) length(x[x!=0]))
EB_temp2 <- data.frame(n_early=EB_temp)

LateBird <- snps_base[,grep("F3_L", names(snps_base))]
LB_temp <- apply(LateBird, 1, function(x) length(x[x!=0]))
LB_temp2 <- data.frame(n_late=LB_temp)

if(unique(row.names(EB_temp2)==row.names(LB_temp2))) ELB_temp <- cbind(EB_temp2, LB_temp2)
ELB_temp$passed <- ifelse(ELB_temp$n_early==10 & ELB_temp$n_late==10, 1, 0)
ELB_temp$passed_9 <- ifelse(ELB_temp$n_early>=9 & ELB_temp$n_late>=9, 1, 0)
ELB_temp$passed_8 <- ifelse(ELB_temp$n_early>=8 & ELB_temp$n_late>=8, 1, 0)

# get "good snps"
Good_snps_9 <- row.names(ELB_temp[ELB_temp$passed_9==1,])


# prepare for QC in GenABEL
SampleInfo <- read.table("src/SampleData.txt", sep="\t", header=T)
SampleData <- SampleInfo[SampleInfo$Pool==1 & SampleInfo$Run==1,c(2,4,5,8)]
names(SampleData)[1] <- "sample.id"
Phen_dat <- SampleData[match(samples_raw$V2, SampleData$sample.id),]
Phen_dat$sex <- rep(0,nrow(Phen_dat))
write.table(Phen_dat, "plink/individuals_raw/phenotype.dat", col.names=T, row.names=F, quote=F)

convert.snp.tped(tpedfile="plink/individuals_raw/ERC_wgbs_trans.tped",tfamfile="plink/individuals_raw/ERC_wgbs_trans.tfam",out="plink/individuals_raw/genotype.raw")
genotype.data0 <- load.gwaa.data(phe="plink/individuals_raw/phenotype.dat",
                                 gen="plink/individuals_raw/genotype.raw",
                                 id = "sample.id")

# remove snps with call rate<=0.9; 1,556,663 snps remain
genotype.data <- genotype.data0[, Good_snps_9]
save(genotype.data, file="plink/individuals_.9_QC_ok/gwaa.data.RData")

QC1 <- check.marker(genotype.data, p.level=0, maf=0.125, callrate=0, perid.call=0.9)
write.table(summary(QC1)$`Per-SNP fails statistics`, "plink/individuals_.9_QC_ok/QC_summary", col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

# 451,600 snps passed QC
write.table(QC1$snpok, "plink/individuals_.9_QC_ok/QC_good_snps", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
Indv.good_snps <- QC1$snpok


## population structure
Cdata <- genotype.data[, QC1$snpok]
Cdata.gkin <- ibs(Cdata, weight = "freq")
Cdata.dist <- as.dist(0.5-Cdata.gkin)
Cdata.mds <- as.data.frame(cmdscale(Cdata.dist))
Cdata.mds$sample.id <- row.names(Cdata.mds); row.names(Cdata.mds) <- 1:nrow(Cdata.mds)
Cdata.mds.phen <- merge(Phen_dat[,c(1,2,4)], Cdata.mds, by="sample.id")

PhenoInfo <- read.table("../WGBS_Snakemake_Bismark/out/PhenotypeInfo", sep="\t", header=T)
helper <- str_split_fixed(Phen_dat$sample.id, "_", 4)
helper2 <- data.frame(sample.id=SampleData$sample.id, Ring.number=paste(helper[,3], helper[,4], sep="..."))
PhenoInfo2 <- merge(PhenoInfo[,c(1,6:8)], helper2, by="Ring.number")
Cdata.mds.phen_v0 <- merge(Cdata.mds.phen, PhenoInfo2[,-1], by="sample.id")

write.table(Cdata.mds.phen_v0[,-c(3,6)], "plink/individuals_.9_QC_ok/mds.phen", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# get variance explained
eig <- cmdscale(Cdata.dist, eig = T)$eig
Var_explained <- round(eig*100/sum(eig),2)
write.table(Var_explained, "plink/individuals_.9_QC_ok/mds.Var_explained", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")

# make plot
col <- c("azure4", "darkslategrey")
order <- as.vector(unique(Cdata.mds.phen$Line))
name <- c("Early", "Late")

Plot_Line <- Cdata.mds.phen_v0 %>% mutate(Line = factor(Line, levels = order)) %>%
  ggplot(aes(x=V1, y=V2, col=Line)) +
  geom_point(shape=19, size=5, alpha=1)+
  scale_color_manual(name="", values=col, labels=name) +
  ylab("PC1") + xlab("PC2") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(0.8, "cm"))

col <- c("khaki1", "darkslategrey")
Plot_LD <- Cdata.mds.phen_v0  %>%
  ggplot(aes(x=V1, y=V2, col=LayDateApril)) +
  geom_point(shape=19, size=5, alpha=1) +
  scale_color_gradient(name = "Laying date\n(April days)", low=col[1], high=col[2], na.value="gray50") +
  ylab(paste("PC1")) + xlab("PC2") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.1, "cm"))

Plot_GEBV <- Cdata.mds.phen_v0  %>%
  ggplot(aes(x=V1, y=V2, col=GEBV)) +
  geom_point(shape=19, size=5, alpha=1) +
  scale_color_gradient(name = "GEBV", low=col[1], high=col[2], na.value="gray50") +
  ylab(paste("PC1")) + xlab("PC2") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.1, "cm"))


Plot <- plot_grid(Plot_Line, Plot_LD, Plot_GEBV, labels = "AUTO", label_size = 20, scale = 0.9, ncol=2)
save_plot("plots/mds_individuals_.9.pdf", Plot, base_height = 10, base_width = 14)