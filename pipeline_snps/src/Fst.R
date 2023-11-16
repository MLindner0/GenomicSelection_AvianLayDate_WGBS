
### --------------------------------------------------------------------------------
###
### SNPs from wgbs data: BayeScan, Arlequin and plink
### Project: WGBS ERC selection lines
###
### Author: M.Lindner (M.Lindner@nioo.knaw.nl), v202107

options(width=200)

library(adegenet)
library(dartR)
library(GenABEL)
library(coda)

library(stringr)

setwd("/Users/melanielindner/Documents/NIOO/Projects/ERC_molecular_data/WBGS_selection_lines/wgbs_snakemake_snps")

### Format data for BayeScan analysis -----------------------------------------
### ---------------------------------------------------------------------------

# read snps:
# filter: per indiv call rate = 0.9, maf=0.125 (2 indvs)
snps <- read.PLINK("plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs.raw")
groups <- read.table("src/groups")
if(unique(snps$ind.names==groups$V2)) pop(snps) <- groups$V3

# flags missing, so reset them:
source("src/utils.reset.flags.R")
snps.flags <- utils.reset.flags(snps, set = FALSE, value = 2)
!is.null(snps.flags@other$loc.metrics)
!is.null(snps.flags@other$loc.metrics.flags)

gl2bayescan(snps.flags, outfile="plink/individuals_.9_QC_ok_Bayescan/ERC_wgbs.bayescan.txt", outpath=".")


### Analyses ------------------------------------------------------------------

# analyses are run in bash


### Evaluate results from BayeScan, Arlequin and plink ------------------------
# based on 451,600 snps after QC (filter: perid=0.9 and maf=0.125)

## annotations
# get snps and some snp-specific details:
snp_names <- snps@loc.names
temp <- str_split_fixed(snp_names, ":", 4)

# fix snp name
help <- str_split_fixed(snp_names, "_", 3)
snp_names.0 <- paste(help[,1], "_", help[,2], sep="")

snp_locations <- data.frame(chr=temp[,1], start=as.numeric(temp[,2]), stop=as.numeric(temp[,2])+1, snp=snp_names.0)
write.table(snp_locations, "temp/snp_locations", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# get annotations in bash and load
snp_annotation <- read.table("temp/annotated_Snps", sep="\t")


## (1) BayeScan results 

# Evaluate convergence (from: http://evomics.org/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf)
# get chain:
chain <- read.table("plink/individuals_.9_QC_ok_Bayescan_pr_odds100/ERC_wgbs.bayescan.sel")
chain <- mcmc(chain,thin=10)

# trace plots and plot of posterior distribution
pdf(file = "plots/BayeScan_convergence_plots.pdf", width = 6, height = 8)
plot(chain) 
dev.off()
# summary (empirical mean, quantiles)
summary(chain)
# auto-correlation of samples
autocorr.diag(chain)
# effective sample size (smaller for LogL as quto-correlation degrades slower)
effectiveSize(chain) 

# Geweke’s convergence diagnostic (based on the comparison of the means of the first and last parts of a Markov chain)
# ok; convergence (equality of means) cannot be rejected as z-scores >-1.96 or <1.96
geweke.diag(chain, frac1=0.1, frac2=0.5) 
# Heidelberg and Welch’s convergence diagnostic
heidel.diag(chain, eps=0.1, pvalue=0.05) # ok
# even better; run >1 chain and compare them (compare between- and within-chain variances)

# Evaluate convergence again, but now when more samples are taken (number of outputted iterations = 25,000; default:5000) and higher thinning interval is specified (thinning interval = 50; default:10)
# get chain:
chain <- read.table("plink/individuals_.9_QC_ok_Bayescan_pr_odds100_more_samples/ERC_wgbs.bayescan.sel")
chain <- mcmc(chain,thin=50)

# trace plots and plot of posterior distribution
pdf(file = "plots/BayeScan_convergence_plots_more_samples.pdf", width = 6, height = 8)
plot(chain) 
dev.off()
# summary (empirical mean, quantiles)
summary(chain)
# auto-correlation of samples
autocorr.diag(chain)
# effective sample size (smaller for LogL as quto-correlation degrades slower)
effectiveSize(chain) 

# Geweke’s convergence diagnostic (based on the comparison of the means of the first and last parts of a Markov chain)
# ok; convergence (equality of means) cannot be rejected as z-scores >-1.96 or <1.96
geweke.diag(chain, frac1=0.1, frac2=0.5) 
# Heidelberg and Welch’s convergence diagnostic
heidel.diag(chain, eps=0.1, pvalue=0.05) # ok

# load results
bayeScan.out_raw <- read.table("plink/individuals_.9_QC_ok_Bayescan_pr_odds100_more_samples/ERC_wgbs.bayescan_fst.txt")
bayeScan.out_raw$sign <- ifelse(bayeScan.out_raw$qval<0.05,1,0)
bayeScan.out_raw$snp <- snp_names.0
# add "bayeS_" to columns names:
names(bayeScan.out_raw)[1:6] <- paste("bayeS", names(bayeScan.out_raw)[1:6], sep="_")
Fst_all <- bayeScan.out_raw[,c(7,1:6)]

## (2) Arlequin results 

# load results
Arlequin.out_raw <- read.table("plink/individuals_.9_QC_ok_Arlequin/ERC_wgbs.res/fdist2_ObsOut_no_header.txt", sep="\t", fill=TRUE, row.names=1)
names(Arlequin.out_raw) <- c("Obs_Het", "Obs_FST", "FST_pval", "1-FST_quantile")

Arlequin.out_raw$sign <- ifelse(Arlequin.out_raw$FST_pval<0.05/nrow(Arlequin.out_raw),1,0)
Arlequin.out_raw$sign_0.01 <- ifelse(Arlequin.out_raw$FST_pval<0.01,1,0)
Arlequin.out_raw$snp <- snp_names.0
names(Arlequin.out_raw)[1:6] <- paste("Arlequin", names(Arlequin.out_raw)[1:6], sep="_")

# combine BayeScan and Arlequin results
if(unique(Fst_all$snp==Arlequin.out_raw$snp)) Fst_all <- cbind(Fst_all, Arlequin.out_raw[,-7])


## (3) plink results and stats

# load data
plink.fst.raw <- read.table("plink/individuals_.9_QC_ok/Fst_ERC_wgbs.fst", header=TRUE)

# add to other results
if(unique(Fst_all$snp==plink.fst.raw$SNP)) Fst_all$plink_FST <- plink.fst.raw$FST

# add results from plink's Hardy-Weinberg equilibrium exact test statistic
# output file explanation:
# A text file with a header line, and one line per marker with the following nine fields:
#         
#         CHR	Chromosome code
# SNP	Variant identifier
# TEST	Type of test: one of {'ALL', 'AFF', 'UNAFF', 'ALL(QT)', 'ALL(NP)'}
# A1	Allele 1 (usually minor)
# A2	Allele 2 (usually major)
# GENO	'/'-separated genotype counts (A1 hom, het, A2 hom)
# O(HET)	Observed heterozygote frequency
# E(HET)	Expected heterozygote frequency
# P	Hardy-Weinberg equilibrium exact test p-value

plink.hwe.raw <- read.table("plink/individuals_.9_QC_ok/Hardy_ERC_wgbs.hwe", header=TRUE)
names(plink.hwe.raw)[4:9] <- paste("plink.hwe", names(plink.hwe.raw)[4:9], sep="_")

# add to other results
if(unique(Fst_all$snp==plink.hwe.raw$SNP)) Fst_all <- cbind(Fst_all, plink.hwe.raw[,4:9])

# save data:
save(Fst_all, file="out/Fst_all.RData")
write.table(Fst_all, "out/Fst_all.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="/t")


### Make Manhattan plots

## get snp positions on genome
# first, get chr length and order chromosomes
Chr_length_names <- read.table("../wgbs_snakemake_bismark/genome/chr_length_names", header=F, sep="\t") 
names(Chr_length_names) <- c("Chr_name", "Chr_name_short", "chr", "length")
# remove Z and MT (as not present in data)
Chr_length_names <- Chr_length_names[Chr_length_names$Chr_name_short!="Z" & Chr_length_names$Chr_name_short!="MT",]

order <- c("1", "1A", "2", "3", "4", "4A", as.character(5:15), as.character(17:24), "25LG1", "25LG2", "26", "27", "28", "LGE22", "Sc")
col_code <- rep(c(1,2), length(order)/2)
helper <- data.frame("Chr_name_short"=order, col=as.factor(col_code), index=1:length(order))
Chr_length_names.col <- merge(Chr_length_names, helper, by="Chr_name_short"); Chr_length_ordered <- Chr_length_names.col[order(Chr_length_names.col$index),]
Chr_length_ordered$index <- 1:nrow(Chr_length_ordered)
Chr_length_ordered$CumSum <- cumsum(Chr_length_ordered$length)

# prepare snp positions
help <- data.frame(Fst_all$snp, str_split_fixed(Fst_all$snp, ":", 4)); names(help) <- c("snp", "chr", "pos", "A1", "A2")
help$pos <- as.numeric(help$pos)
snps.positions.0 <- merge(help, Chr_length_ordered, by="chr")
snps.positions <- snps.positions.0[order(snps.positions.0$index, snps.positions.0$pos),]
snps.positions$genome.pos <- snps.positions$CumSum-snps.positions$length+snps.positions$pos

# get axis label
# non-sc part
axis.help.0 <- Chr_length_ordered[Chr_length_ordered$Chr_name_short!="Sc",]
axis.help.0 $lab.pos <- axis.help.0$CumSum-(axis.help.0$length/2)
# sc part
axis.help.1 <- Chr_length_ordered[Chr_length_ordered$Chr_name_short=="Sc",]
sum.sc <- sum(axis.help.1$length)
axis.help.2 <- axis.help.1[nrow(axis.help.1),]
axis.help.2$lab.pos <- axis.help.2$CumSum-(sum.sc/2)
# combine (MT in end)
axis.labs <- rbind(axis.help.0, axis.help.2)

# save snp positions
save(snps.positions, axis.labs, file="out/For.Manahttan.RData")

# prepare plot and set plot parameter
dat.plot <- Fst_all[match(snps.positions$snp, Fst_all$snp),]
if(unique(dat.plot$snp==snps.positions$snp)) dat.plot <- cbind(dat.plot, snps.positions[,c(9,12)])

row.names(axis.labs) <- 1:nrow(axis.labs)
help <- axis.labs$Chr_name_short
axis.labs$x.lab <- c(help[1:10], "", help[12], rep("", 4), help[17], rep("", 7), help[25], rep("", 6), help[32])

x.lim <- c(range(dat.plot$genome.pos)[1]-range(dat.plot$genome.pos)[1], range(dat.plot$genome.pos)[2]+range(dat.plot$genome.pos)[1])

# make plot: (1) Arlequin
pdf(file = "plots/Fst_Arlequin.pdf", width = 12, height = 5)
par(mfrow = c(1,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0), xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(-0.06,0.9), xlab="",
     ylab=expression(paste("Arlequin F"["st"])), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,0.8,0.2), labels=seq(0,0.8,0.2), lwd = 0, lwd.ticks = 1)
points(dat.plot$genome.pos, dat.plot$Arlequin_Obs_FST, pch=19, cex=ifelse(dat.plot$Arlequin_sign==1,1,0.5), col=ifelse(dat.plot$Arlequin_sign==1, "darkorange", ifelse(dat.plot$col==1, "darkslategrey", "azure4")))
dev.off()

png(file = "plots/Fst_Arlequin.png", width = 800, height = 300)
par(mfrow = c(1,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0), xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(-0.06,0.9), xlab="",
     ylab=expression(paste("Arlequin F"["st"])), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,0.8,0.2), labels=seq(0,0.8,0.2), lwd = 0, lwd.ticks = 1)
points(dat.plot$genome.pos, dat.plot$Arlequin_Obs_FST, pch=19, cex=ifelse(dat.plot$Arlequin_sign==1,1,0.5), col=ifelse(dat.plot$Arlequin_sign==1, "darkorange", ifelse(dat.plot$col==1, "darkslategrey", "azure4")))
dev.off()


# make plot: (2) BayeScan
pdf(file = "plots/Fst_bayeScan.pdf", width = 12, height = 5)
par(mfrow = c(1,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0), xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(0.0235,0.25), xlab="",
     ylab=expression(paste("BayeScan F"["st"])), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0.06,0.24,0.06), labels=seq(0.06,0.24,0.06), lwd = 0, lwd.ticks = 1)
points(dat.plot$genome.pos, dat.plot$bayeS_fst, pch=19, cex=ifelse(dat.plot$bayeS_sign==1,1,0.5), col=ifelse(dat.plot$bayeS_sign==1, "darkorange", ifelse(dat.plot$col==1, "darkslategrey", "azure4")))
dev.off()

png(file = "plots/Fst_bayeScan.png", width = 800, height = 300)
par(mfrow = c(1,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0), xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(0.0235,0.25), xlab="",
     ylab=expression(paste("BayeScan F"["st"])), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0.06,0.24,0.06), labels=seq(0.06,0.24,0.06), lwd = 0, lwd.ticks = 1)
points(dat.plot$genome.pos, dat.plot$bayeS_fst, pch=19, cex=ifelse(dat.plot$bayeS_sign==1,1,0.5), col=ifelse(dat.plot$bayeS_sign==1, "darkorange", ifelse(dat.plot$col==1, "darkslategrey", "azure4")))
dev.off()

# make plot like for DMA result:

# helper
dat.plot$bayeS_col_code <- as.factor(ifelse(dat.plot$bayeS_qval<0.05, 3, ifelse(dat.plot$bayeS_qval>0.05 & dat.plot$col==1,1,2)))
base_col <- c("darkslategrey", "azure4", "darkorange")
axis.labs.0 <- axis.labs[c(1:11, 13, 15, 18, 21, 25, 32),]

# make plot
Plot.Manhattan.Fst <- dat.plot %>% 
        ggplot() +
        geom_point(aes(x=genome.pos, y=bayeS_fst, color=bayeS_col_code, size=bayeS_col_code), shape=19) +  
        scale_x_continuous(breaks=axis.labs.0$lab.pos, expand=c(0,0), labels=axis.labs.0$Chr_name_short) +
        scale_y_continuous(expand=c(0,0), breaks=seq(0.06,0.24,0.06), limits=c(0.023,0.25)) +
        scale_color_manual(name="", values=base_col) +
        scale_size_manual(name="", values=c(2,2,3)) +
        ylab(expression(paste("BayeScan F"["st"]))) + xlab("Chromosome position") +
        theme_classic() +
        theme(axis.text.x=element_text(size=20, margin=margin(t = 5, r = 0, b = 0, l = 0)),
              axis.text.y=element_text(size=20),
              axis.title.x=element_text(size=22, margin=margin(t = 10, r = 0, b = 0, l = 0)),
              axis.title.y=element_text(size=22, margin=margin(t = 0, r = 7, b = 0, l = 0)),
              axis.line=element_line(colour="black"),
              plot.margin=unit(c(0.2,0.5,0.2,0.3),"cm"),
              legend.position = "none") 
save_plot("plots/manhattan.Fst_plot_Final_v1.png", Plot.Manhattan.Fst, base_height = 4.5, base_width = 15)
save_plot("plots/manhattan.Fst_plot_Final_v1.pdf", Plot.Manhattan.Fst, base_height = 4.5, base_width = 15)

# make plot: (3) plink
pdf(file = "plots/Fst_plink.pdf", width = 12, height = 5)
par(mfrow = c(1,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0),  xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(-0.13,0.9), xlab="",
     ylab=expression(paste("plink F"["st"])), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,0.8,0.2), labels=seq(0,0.8,0.2), lwd = 0, lwd.ticks = 1)
points(dat.plot$genome.pos, dat.plot$plink_FST, pch=19, cex=0.5, col=ifelse(dat.plot$col==1, "darkslategrey", "azure4"))
dev.off()

png(file = "plots/Fst_plink.png", width = 800, height = 300)
par(mfrow = c(1,1))
par(mar = c(4,4,2,2), mgp=c(2, 0.5, 0),  xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(-0.15,0.9), xlab="",
     ylab=expression(paste("plink F"["st"])), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,0.8,0.2), labels=seq(0,0.8,0.2), lwd = 0, lwd.ticks = 1)
points(dat.plot$genome.pos, dat.plot$plink_FST, pch=19, cex=0.5, col=ifelse(dat.plot$col==1, "darkslategrey", "azure4"))
dev.off()

## prepare tables for manuscript

# get outlier SNPs:
Sign_SNPs <- bayeScan.out_raw[bayeScan.out_raw$bayeS_sign==1,]
Annotations.Fst.out <- snp_annotation[snp_annotation$V4 %in% Sign_SNPs$snp,]

# for main:
main.fst <- merge(Sign_SNPs[,c(7,5,3,1)], Annotations.Fst.out[,c(4,11,9)], by.x="snp", by.y="V4", all.x=TRUE)
names(main.fst)[5:6] <- c("Genomic region", "Gene symbol")
main.fst$`Gene symbol`<- gsub("LOC107203824", "LOC107203824 (SOX3)", main.fst$`Gene symbol`)
main.fst$snp <- gsub("NC_031771.1", "chr4", main.fst$snp)
main.fst$snp <- gsub("NC_031772.1", "chr4A", main.fst$snp)
write.table(main.fst, "out/out_Fst.txt", quote=F, col.names=TRUE, row.names=F, sep="\t")

IN <- load(file="plink/individuals_.9_QC_ok/gwaa.data.RData")
genotype.data.sign <- genotype.data[,Sign_SNPs]
genotypes <- as.genotype(genotype.data.sign)
genotypes$Line <- genotype.data.sign@phdata$Line

# which genes?
snp_locations[snp_locations$snp %in% Sign_SNPs,]
snp_annotation[snp_annotation$V4 %in% Sign_SNPs,]


## overlap between tools:

# make plot
pdf(file = "plots/Fst_correlation_scross_tools.pdf", width = 6, height = 6)
par(mfrow=c(2,2), mar=c(4.1,5.1,4.1,2.1), xaxs = "i", yaxs = "i")
plot(dat.plot$bayeS_fst, dat.plot$Arlequin_Obs_FST, pch=20, cex.lab=1.5, cex.axis=1.5, bty = "l", 
     xlim=c(-0.1255, 0.9), ylim=c(-0.1255, 0.9),
     xlab=expression("F"["st"]*" bayeScan"), ylab=expression("F"["st"]*" Arlequin"),
     col=ifelse(dat.plot$bayeS_qval<0.05, "darkorange", "darkslategrey"),
     cex=ifelse(dat.plot$bayeS_qval<0.05, 1.5, 1))
abline(0,1, lty=2)
plot(dat.plot$bayeS_fst, dat.plot$plink_FST, pch=20, cex.lab=1.5, cex.axis=1.5, bty = "l",
     xlim=c(-0.1255, 0.9), ylim=c(-0.1255, 0.9),
     xlab=expression("F"["st"]*" bayeScan"), ylab=expression("F"["st"]*" plink"),
     col=ifelse(dat.plot$bayeS_qval<0.05, "darkorange", "darkslategrey"),
     cex=ifelse(dat.plot$bayeS_qval<0.05, 1.5, 1))
abline(0,1, lty=2)
plot(dat.plot$Arlequin_Obs_FST, dat.plot$plink_FST, pch=20, cex.lab=1.5, cex.axis=1.5, bty = "l",
     xlim=c(-0.1255, 0.9), ylim=c(-0.1255, 0.9),
     xlab=expression("F"["st"]*" Arlequin"), ylab=expression("F"["st"]*" plink"),
     col=ifelse(dat.plot$bayeS_qval<0.05, "darkorange", "darkslategrey"),
     cex=ifelse(dat.plot$bayeS_qval<0.05, 1.5, 1))
abline(0,1, lty=2)
dev.off()


## compare Arlequin output to Arlequin analysis in previous study
# Verhagen, I, Gienapp, P, Laine, VN, et al. Genetic and phenotypic responses to genomic selection for timing of breeding in a wild songbird. Funct Ecol. 2019; 33: 1708–1721. https://doi.org/10.1111/1365-2435.13360

# get genes from vero's analysis:
SNPchip <- read.csv("temp/Verhagen_et_al_TableS3.csv") # from previous study
Chr_length_names <- read.table("../wgbs_snakemake_bismark/genome/chr_length_names", header=F, sep="\t") 
names(Chr_length_names) <- c("Chr_name", "Chr_name_short", "chr", "length")

# all Chr names of SNP chip data present in translation file?
SNPchip.chr <- unique(SNPchip$CHR)
Chr <- unique(Chr_length_names$Chr_name_short)
length(Chr[Chr %in% SNPchip.chr])==length(SNPchip.chr) # Yes

change.chr.names <- Chr_length_names[,c(2,3)]; names(change.chr.names)[1] <- "CHR"
SNPchip.new <- merge(SNPchip[c(1,3)], change.chr.names, by="CHR")
SNPchip.out <- SNPchip.new[,c(3,2)]
SNPchip.out$POS2 <- SNPchip.out$POS+1
write.table(SNPchip.out, "temp/SNP.chip_SNP.positions.txt", quote=F, col.names=F, row.names=F, sep="\t")

# get overlap with annotations in bash and load annotations
SNPchip.annotations <- read.table("temp/SNPchip.annotated_Snps")
SNPchip.Genes <- unique(SNPchip.annotations$V8)

# get genes from my Arlequin analysis:
Arlequin.sig0.01 <- Arlequin.out_raw[Arlequin.out_raw$Arlequin_FST_pval<0.01,]
Arlequin.annotations <- snp_annotation[snp_annotation$V4 %in% Arlequin.sig0.01$snp, c(4,9)]
Arlequin.genes <- unique(Arlequin.annotations$V9)

overlapp <- intersect(Arlequin.genes, SNPchip.Genes)
Only.in.SNPchip <- SNPchip.Genes[!SNPchip.Genes %in% overlapp]
Only.in.wgbs <- Arlequin.genes[!Arlequin.genes %in% overlapp]
write.table(overlapp, "out/Genes.overlapp.SNPchip.txt", quote=F, col.names=F, row.names=F, sep="\n")

# check for SNPs that are significant in both
help <- str_split_fixed(Arlequin.sig0.01$snp, ":", 3)
Arlequin.sig0.01_comparison <- as.data.frame(help[,1:2])
wgbs_sign <- paste(Arlequin.sig0.01_comparison$V1, Arlequin.sig0.01_comparison$V2, sep="_")
SNPchip_sign <- paste(SNPchip.out$chr, SNPchip.out$POS, sep="_")
shared_outlier <- intersect(wgbs_sign, SNPchip_sign)

an.0 <- snp_annotation
an.0$help <- paste(an.0$V1, an.0$V2, sep="_")
an.0.shared <- an.0[an.0$help %in% shared_outlier,]
write.table(unique(an.0.shared$V9), "out/shared_Fst_outliers_v2.genes.txt", quote=F, col.names=F, row.names=F, sep="\n")

# check how many snps are in both data sets (irrespective of significance)
SNPchip_all <- read.table("temp/gt_NLallSept2018_nobadsnps_noZUns_SELECTION_forFST.bim") # input data from previous analysis (which is not publicly avaialable)
SNPchip_all <- SNPchip_all[,c(1,4)]
names(SNPchip_all)[1] <- "CHR"
SNPchip_all.new <- merge(SNPchip_all, change.chr.names, by="CHR")

SNPchip_all <- paste(SNPchip_all.new$chr, SNPchip_all.new$V4, sep="_")

help <- str_split_fixed(Arlequin.out_raw$snp, ":", 3)
Arlequin_comparison <- as.data.frame(help[,1:2])
wgbs_all <- paste(Arlequin_comparison$V1, Arlequin_comparison$V2, sep="_")

shared_snp_positions <- intersect(wgbs_all, SNPchip_all)

shared_snp_positions_sign_SNPchip <- intersect(shared_snp_positions, SNPchip_sign)
shared_snp_positions_sign_wgbs <- intersect(shared_snp_positions, wgbs_sign)

shared_snp_positions_sign_SNPchip_check <- intersect(shared_snp_positions_sign_SNPchip, shared_outlier)
shared_snp_positions_sign_wgbs_check <- intersect(shared_snp_positions_sign_wgbs, shared_outlier)

