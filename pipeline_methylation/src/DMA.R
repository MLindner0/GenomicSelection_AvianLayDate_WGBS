## ----------------------------------------------------------------------
###
###
### Project: Molecular signatures of selection
###
### Author: M.Lindner (M.Lindner@nioo.knaw.nl), v2021-02


### -------------- load required packages & set up environment -----------

options(width=200)

# analysis of methylation data:
library(methylKit)
library(lme4)

# data formatting and visualization:
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# genome annotation:
library(GenomicFeatures)
library(rtracklayer)
library(genomation)

# set working directory
base <- "/home/NIOO.INT/melaniel/projects/WGBS_Snakemake_Bismark"
setwd(base)


### -------------- load and filter methylation data -----------------------

# get file names
FileLocation <- "out/meth"
FileNames <- list.files(path = FileLocation, pattern = "*CpG_report.txt")

# select methylation data of individual samples (we also sequenced pool which are ignored here)
FileNames_ind <- c(FileNames[grep("F3_E_BD", FileNames)], FileNames[grep("F3_L_BD", FileNames)])

# get sample data
SampleData <- read.table("src/SampleData.txt", header=T, sep="\t")

# remove tissue sample
Out <- as.character(SampleData[SampleData$SampleType=="Animal","sample.id"])
FileNames_ind <- FileNames_ind[-grep(Out, FileNames_ind)]

# get file location
FileNames_Location <- paste(FileLocation, FileNames_ind, sep="/")
file_list <- as.list(FileNames_Location)

## Make methylkit object
sample.id <- str_split_fixed(FileNames_ind, "\\.", 4)[,1] # sample id
sample_id <- as.list(sample.id)
if(unique(sample.id==SampleData$sample.id)) Line_numeric <- SampleData$Line_numeric2 # selection line (i.e. treatment)

myobj <- methRead(file_list, pipeline="bismarkCytosineReport", sample.id=sample_id, assembly="p.major1.1", treatment=as.numeric(Line_numeric), context = "CpG", mincov=1)

# apply coverage and percentile filter
filtered.myobj <- filterByCoverage(myobj, lo.count=1, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

# unite samples and combine cytosines of both strands into CpG sites. For this, destrand set to TRUE as there is strand info to combine Cs. Minimum coverage per x samples within treatment group (i.e. selectio line) is set below, hence min.per.group=1L to not filter any CpG sites
meth <- unite(filtered.myobj, destrand=TRUE, min.per.group=1L)
save(meth, file="temp/Methylkit_meth_v0.RData") # save data

## Format data
## to (1) only keep sites with at least 10x per CpG site in all samples, (2) get all sites with 0/100% methylation across samples, (3) calculate mean coverage and mean difference, and (4) prepare PCA

# get data frames with methylation proportion, total counts, and methylated counts
methData <- getData(meth)
sample_id <- meth@sample.ids
site <- paste(methData$chr, methData$start, sep="_")
methDataNew <- data.frame(site=site)
methDataNew_Cov <- data.frame(site=site)
methDataNew_Meth <- data.frame(site=site)

for(f in seq(5,59,3)){
  cov <- methData[,f]
  methCount <- methData[,(f+1)]
  methDataNew$methProp <- methCount/cov
  methDataNew_Cov$cov <- cov
  methDataNew_Meth$meth <- methCount
  no <- (f-2)/3
  names(methDataNew)[no+1] <- sample_id[no]
  names(methDataNew_Cov)[no+1] <- sample_id[no]
  names(methDataNew_Meth)[no+1] <- sample_id[no]
}
row.names(methDataNew) <- methDataNew$site
row.names(methDataNew_Cov) <- methDataNew_Cov$site
row.names(methDataNew_Meth) <- methDataNew_Meth$site

# (1) only keep sites with at least 10x per CpG site in all samples
temp <- methDataNew_Cov
temp$site <- NULL
temp[temp<10]=NA

EarlyBird <- temp[,grep("F3_E", names(temp))]
EB_temp <- apply(EarlyBird, 1, function(x) length(x[!is.na(x)]))
EB_temp2 <- data.frame(n_early=EB_temp)

LateBird <- temp[,grep("F3_L", names(temp))]
LB_temp <- apply(LateBird, 1, function(x) length(x[!is.na(x)]))
LB_temp2 <- data.frame(n_late=LB_temp)

if(unique(row.names(EB_temp2)==row.names(LB_temp2))) ELB_temp <- cbind(EB_temp2, LB_temp2)
ELB_temp$passed <- ifelse(ELB_temp$n_early==10 & ELB_temp$n_late==9, 1, 0)
GoodSites <- row.names(ELB_temp[ELB_temp$passed==1,])

# update data:
methDataNew <- methDataNew[methDataNew$site %in% GoodSites,]
methDataNew_Cov <- methDataNew_Cov[methDataNew_Cov$site %in% GoodSites,]
methDataNew_Meth <- methDataNew_Meth[methDataNew_Meth$site %in% GoodSites,]

# (2) get and remove all sites with 0/100% methylation across samples
Cov_temp <- methDataNew_Cov[, -1]
Cov_sum <- apply(Cov_temp, 1, sum)
Cov_sum_df <- as.data.frame(Cov_sum)

Meth_temp <- methDataNew_Meth[, -1]
Meth_sum <- apply(Meth_temp, 1, sum)
Meth_sum_df <- as.data.frame(Meth_sum)

if(unique(row.names(Cov_sum_df)==row.names(Meth_sum_df))) Out_temp <- cbind(Cov_sum_df, Meth_sum_df)
Out_temp$Diff_sum <- Out_temp$Cov_sum-Out_temp$Meth_sum
Out_temp$out <- ifelse(Out_temp$Meth_sum==0 | Out_temp$Diff_sum==0, 1, 0)
GoodSites2 <- row.names(Out_temp[Out_temp$out==0,])

# update data:
methDataNew <- methDataNew[methDataNew$site %in% GoodSites2,]
methDataNew_Cov <- methDataNew_Cov[methDataNew_Cov$site %in% GoodSites2,]
methDataNew_Meth <- methDataNew_Meth[methDataNew_Meth$site %in% GoodSites2,]
save(methDataNew, methDataNew_Cov, methDataNew_Meth, file="temp/Meth_Dat_filtered_v0.RData") # save updated data

# (3) calculate mean difference and mean coverage
EarlyBird <- methDataNew[,grep("F3_E", names(methDataNew))]
EB_temp <- apply(EarlyBird, 1, mean)
EB_temp2 <- data.frame(Early_mean=EB_temp)

LateBird <- methDataNew[,grep("F3_L", names(methDataNew))]
LB_temp <- apply(LateBird, 1, mean)
LB_temp2 <- data.frame(Late_mean=LB_temp)

if(unique(row.names(EB_temp2)==row.names(LB_temp2))) ELB_temp <- cbind(EB_temp2, LB_temp2)
ELB_temp$Diff_Meth_E_L <- ELB_temp$Early_mean-ELB_temp$Late_mean

Cov_temp <- methDataNew_Cov[, -1]
Cov_sum <- apply(Cov_temp, 1, sum)
Cov_sum_df <- as.data.frame(Cov_sum)
Cov_sum_df$Cov_mean <- Cov_sum_df$Cov_sum/19

if(unique(row.names(ELB_temp)==row.names(Cov_sum_df))) Meth_Dat_Info <- cbind(ELB_temp, Cov_sum_df)
save(Meth_Dat_Info, file="temp/Meth_Dat_Info_v0.RData")

# (4) PCA
# run PCA
Dat_PCA <- methDataNew
Dat_PCA$site <- NULL
PCA <- prcomp(Dat_PCA, center=F, scale.=F)
pcaData <- data.frame(PCA$rotation)
if(unique(row.names(pcaData)==SampleData$sample.id)) PCA_Data <- cbind(pcaData, SampleData)
save(PCA_Data, file="temp/PCA_Dat_v0.RData")

# calculate var explaines
eigs <- PCA$sdev^2
var <- eigs/sum(eigs)
VarExplained <- data.frame(PC=paste("PC", 1:19, sep=""), Var=var)
write.table(VarExplained, "plots/PCA_Var_explained.txt", quote=F, sep="\t", row.names=F, col.names=T)

# plot PCA output differentiated by
base_col <- c("darkorange", "darkslategrey")
colfunc <- colorRampPalette(base_col)

# (1) tissue type
order <- as.vector(unique(PCA_Data$SampleType))
name <- c("Blood\n(no\nplasma)", "Blood")
col <- base_col

Plot1 <- PCA_Data %>% mutate(SampleType = factor(SampleType, levels = order)) %>% 
  ggplot(aes(x=PC2, y=PC1, color=SampleType)) +
  geom_point(shape=19, size=5, alpha=1)+
  scale_color_manual(name="", values=col, labels=name) +
  ylab(paste("PC1 (95.67%)")) + xlab("PC2 (0.30%)") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"))

# (2) selection line
order <- as.vector(unique(PCA_Data$Line))
name <- c("Early", "Late")
col <- base_col

Plot2 <- PCA_Data %>% mutate(Line = factor(Line, levels = order)) %>% 
  ggplot(aes(x=PC2, y=PC1, color=Line)) +
  geom_point(shape=19, size=5, alpha=1)+
  scale_color_manual(name="", values=col, labels=name) +
  ylab(paste("PC1 (95.67%)")) + xlab("PC2 (0.30%)") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.1, "cm"))

# (3) reproductive status
order <- as.vector(unique(PCA_Data$Status))
name <- c("Laying", "Post-\nlaying", "Not\nlaying")
col <- c(base_col, "gray70")

Plot3 <- PCA_Data %>% mutate(Status = factor(Status, levels = order)) %>% 
  ggplot(aes(x=PC2, y=PC1, color=Status)) +
  geom_point(shape=19, size=5, alpha=1)+
  scale_color_manual(name="", values=col, labels=name) +
  ylab(paste("PC1 (95.67%)")) + xlab("PC2 (0.30%)") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.key.size = unit(1.5, "cm"))

# (4) lay date
col <- base_col
Plot4 <- PCA_Data %>% 
  ggplot(aes(x=PC2, y=PC1, color=LayDateApril)) +
  geom_point(shape=19, size=5, alpha=1)+
  scale_color_gradient(name = "Lay date", low=col[1], high=col[2], na.value="gray70") +
  ylab(paste("PC1 (95.67%)")) + xlab("PC2 (0.30%)") +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1.1, "cm"))

Plot <- plot_grid(Plot2, Plot1, Plot3, Plot4, labels = "AUTO", label_size = 20, scale = 0.9, ncol=2, nrow=2)
save_plot("plots/PCA.pdf", Plot, base_height = 8, base_width = 11)



### -------------- differential methylation analysis ---------------------

## differential methylation analysis
# run methylkit:
diffMeth <- calculateDiffMeth(meth_clean, adjust="fdr", overdispersion="MN", mc.cores = 20)

# get results:
Dat_diffMeth <- getData(diffMeth)
Dat_diffMeth$site <- paste(Dat_diffMeth$chr, Dat_diffMeth$start, sep="_")

# get differentially methylated sites (DMS):
DMS_individuals <- Dat_diffMeth[Dat_diffMeth$qvalue<0.05,]$site

save(Dat_diffMeth, DMS_individuals, file="temp/MethylKit_dma_v0.RData") # save data


## make validation plots
# get sign thresholds (fdr and bonferroni)
dummy <- Dat_diffMeth[Dat_diffMeth$qvalue>0.05,]
sign_threshold_fdr <- -log10(min(dummy$pvalue)) # fdr
sign_threshold <- -log10(0.05/nrow(Dat_diffMeth)) # bf

# get functions for making qq-plots
gg_qqplot <- function(ps, ci = 0.95, fill=fill, col=col, al=al) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], "(p-value)"))
  log10Po <- expression(paste("Observed -log"[10], "(p-value)"))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 19, size = 3, color=col, alpha=al) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}
inflation <- function(pvalue) {
  chisq <- qchisq(1 - pvalue, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}

# Make plots:
# (1) p-value distribution
col <- "darkslategrey"
col <- "gray30"
# p-values:
Plot1 <- Dat_diffMeth %>% 
  ggplot() +
  geom_histogram(aes(x = pvalue), position="identity", binwidth=0.02, fill=col, color="white") +  
  ylab("Frequency") + xlab("p-value") +
  theme_classic() +
  theme(axis.text=element_text(size=18),
        axis.title.y=element_text(size=20, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x=element_text(size=20, margin=margin(t = 10, r = 0, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none") 

# (2) coverage vs. -log10(p-values)
if(unique(Dat_diffMeth$site==row.names(Meth_Dat_Info))) Dat_diffMeth$Cov_mean <- Meth_Dat_Info$Cov_mean
Dat_diffMeth$log10_pvalue <- -log10(Dat_diffMeth$pvalue)

Plot2 <- Dat_diffMeth %>% 
  ggplot(aes(x=Cov_mean, y=log10_pvalue)) +
  geom_point(shape=19, size=3, alpha=1, color=col) + 
  ylab(expression(paste("-log"["10"],"(p-value)"))) + xlab("Mean coverage across samples") +
  theme_classic() +
  theme(axis.text=element_text(size=18),
        axis.title.y=element_text(size=20, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x=element_text(size=20, margin=margin(t = 10, r = 0, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

# (2) coverage vs. p-values
Plot2b <- Dat_diffMeth %>% 
  ggplot(aes(x=Cov_mean, y=pvalue)) +
  geom_point(shape=19, size=3, alpha=1, color=col) + 
  ylab(expression(paste("p-value"))) + xlab("Mean coverage across samples") +
  theme_classic() +
  theme(axis.text=element_text(size=18),
        axis.title.y=element_text(size=20, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x=element_text(size=20, margin=margin(t = 10, r = 0, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

# (3) qq plot
Plot3 <- gg_qqplot(Dat_diffMeth$pvalue, col=col, al=1) +
  scale_y_continuous(limits=c(0,10)) +
  theme_classic() +
  theme(axis.text=element_text(size=18),
        axis.title.y=element_text(size=20, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x=element_text(size=20, margin=margin(t = 10, r = 0, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none")  
(lambda <- inflation(Dat_diffMeth$pvalue)) # 1.14619

# (4) difference in methylation between selection lines vs. -log10(p-values):
Plot4 <- Dat_diffMeth %>% 
  ggplot(aes(x=meth.diff, y=log10_pvalue)) +
  geom_hline(yintercept=sign_threshold, linetype=2, size=0.5, color="black") +
  annotate(geom="text", x=-61, y=sign_threshold+0.3, label="BF", size=6) +
  geom_hline(yintercept=sign_threshold_fdr, linetype=2, size=0.5, color="black") +
  annotate(geom="text", x=-59, y=sign_threshold_fdr+0.3, label="FDR", size=6) +
  geom_vline(xintercept=-10, linetype=2, size=0.5, color="black") +
  geom_vline(xintercept=10, linetype=2, size=0.5, color="black") +
  geom_point(shape=19, size=3, alpha=1, color=col) + 
  ylab(expression(paste("-log"["10"],"(p-value)"))) + xlab(expression(paste(Delta, " Methylation level", sep=" "))) +
  theme_classic() +
  theme(axis.text=element_text(size=18),
        axis.title.y=element_text(size=20, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x=element_text(size=20, margin=margin(t = 10, r = 0, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none") 

# save plot (.png and .pdf)
Plot <- plot_grid(Plot1, Plot2, Plot3, Plot4, labels = "AUTO", label_size = 20, scale = 0.9, ncol=2)
save_plot("plots/DMA_Validation_MethyLkit_Final_v1.png", Plot, base_height = 14, base_width = 14)
Plot <- plot_grid(Plot1, Plot2, Plot3, Plot4, labels = "", label_size = 20, scale = 0.9, ncol=2)
save_plot("plots/DMA_Validation_MethyLkit_Final_v1.pdf", Plot, base_height = 14, base_width = 14)


# add more details to (2) coverage vs. -log10(p-values)
test <- Dat_diffMeth[,8:10]
range(test$Cov_mean) # 10-190

# split data into coverage bins of length 10 starting
bin.borders <- seq(10,180,10)

test.split <- list()
for(i in 1:length(bin.borders)) {
  part <- test[test$Cov_mean>bin.borders[i] & test$Cov_mean<=(bin.borders[i]+10),]
  test.split[[i]] <- part
}

# get mean -log10(p-value) and CI for each coverage bins
window.median.cov <- bin.borders + 5
mean.log.p <- unlist(lapply(test.split, function(x) mean(x$log10_pvalue)))
library(rethinking)
CI.log.p <- unlist(lapply(test.split, function(x) PI(x$log10_pvalue, prob=0.95)))
n.log.p <- unlist(lapply(test.split, function(x) length(x$log10_pvalue)))
# 16th, 17th and 18th window only have 4, 0, 1 CPGs, espectively --> remove

# combine data:
cov.windows.dat <- data.frame(window.median.cov=window.median.cov, window.mean.logP=mean.log.p, window.low95CI.logP=CI.log.p[seq(1,length(CI.log.p),2)], window.up95CI.logP=CI.log.p[seq(2,length(CI.log.p),2)], window.n.CpGs=n.log.p)
cov.windows.dat_sub <- cov.windows.dat[cov.windows.dat$window.n.CpGs>4,]

# make a quick test:
mod <- lm(Dat_diffMeth$log10_pvalue ~ Dat_diffMeth$Cov_mean)
summary(mod)

# plot data
help.n.scale <- max(cov.windows.dat$window.n.CpGs)/max(Dat_diffMeth$log10_pvalue)

pdf(file = "plots/Revision_DMA_Validation_COVERAGE.pdf", width = 5, height = 5)
par(mfrow=c(1,1), mar = c(4,4,2,4), mgp=c(2, 0.5, 0))

# make plot base
plot(NULL, xlim=range(Dat_diffMeth$Cov_mean), ylim=range(Dat_diffMeth$log10_pvalue), xlab="Mean coverage across samples",
     ylab=expression(paste("-log"["10"],"(p-value)")))
axis(4, at=seq(100,600,100)/(help.n.scale/1000), labels=paste(seq(100,600,100),"k", sep=""), lwd = 0, lwd.ticks = 1)
mtext("n CpGs", side=4, line=2)
# add data
points(Dat_diffMeth$Cov_mean, Dat_diffMeth$log10_pvalue, pch=21, col=rgb(47/255, 79/255, 79/255, 0.25))
# add arrow for Cis
arrows(x0=cov.windows.dat_sub$window.median.cov, y0=cov.windows.dat_sub$window.low95CI.logP, x1=cov.windows.dat_sub$window.median.cov, y1=cov.windows.dat_sub$window.up95CI.logP, code=3, angle=90, length=0.05, col="black", lwd=1)
# add mean
points(cov.windows.dat_sub$window.median.cov, cov.windows.dat_sub$window.mean.logP, pch=19, cex=1.5, col="black")
# add line for n of CpGs within window
lines(cov.windows.dat_sub$window.median.cov, cov.windows.dat_sub$window.n.CpGs/(help.n.scale), lwd=1, lty=2, col=lines(cov.windows.dat_sub$window.median.cov, cov.windows.dat_sub$window.n.CpGs/(help.n.scale), lwd=2, lty=2, col="darkorange")
)
dev.off()

# same but save as .png
png(file = "plots/Revision_DMA_Validation_COVERAGE.png", width = 300, height = 300)
par(mfrow=c(1,1), mar = c(4,4,2,4), mgp=c(2, 0.5, 0))

plot(NULL, xlim=range(Dat_diffMeth$Cov_mean), ylim=range(Dat_diffMeth$log10_pvalue), xlab="Mean coverage across samples",
     ylab=expression(paste("-log"["10"],"(p-value)")))
axis(4, at=seq(100,600,100)/(help.n.scale/1000), labels=paste(seq(100,600,100),"k", sep=""), lwd = 0, lwd.ticks = 1)
mtext("n CpGs", side=4, line=2)
points(Dat_diffMeth$Cov_mean, Dat_diffMeth$log10_pvalue, pch=21, col=rgb(47/255, 79/255, 79/255, 0.25))
arrows(x0=cov.windows.dat_sub$window.median.cov, y0=cov.windows.dat_sub$window.low95CI.logP, x1=cov.windows.dat_sub$window.median.cov, y1=cov.windows.dat_sub$window.up95CI.logP, code=3, angle=90, length=0.05, col="black", lwd=1)
points(cov.windows.dat_sub$window.median.cov, cov.windows.dat_sub$window.mean.logP, pch=19, cex=1.5, col="black")
lines(cov.windows.dat_sub$window.median.cov, cov.windows.dat_sub$window.n.CpGs/(help.n.scale), lwd=1, lty=2, col=lines(cov.windows.dat_sub$window.median.cov, cov.windows.dat_sub$window.n.CpGs/(help.n.scale), lwd=2, lty=2, col="darkorange")
)
dev.off()


## test for difference in methylation level of DMS between selection lineslines
# format data:
new.Cov <- gather(methDataNew_Cov, names(methDataNew_Cov)[-1],
                  key = "Sample", value = "Cov")
new.Meth <- gather(methDataNew_Meth, names(methDataNew_Meth)[-1],
                   key = "Sample", value = "Meth")
if(unique(new.Cov$site==new.Meth$site & new.Cov$Sample==new.Meth$Sample)) new.Cov$Meth <- new.Meth$Meth
new.Cov$Line <- as.factor(ifelse(grepl("F3_E", new.Cov$Sample), 1, 2))
new.Cov$unMeth <- new.Cov$Cov-new.Cov$Meth
new.Cov$Sample_f <- as.factor(new.Cov$Sample)

# filter for significant sites:
DMS.fdr <- Dat_diffMeth[Dat_diffMeth$qvalue<0.05,]$site
new.Cov.DMS <- new.Cov[new.Cov$site %in% DMS.fdr,]

# fit model and null model:
glmm.DMS <- glmer(cbind(Meth, unMeth) ~ Line + (1|site) + (1|Sample_f), new.Cov.DMS, family=binomial)
glmm0.DMS <- glmer(cbind(Meth, unMeth) ~ (1|site) + (1|Sample_f), new.Cov.DMS, family=binomial)

# compare models; p-val=9.429e-10
DMS.comp <- anova(glmm0.DMS, glmm.DMS, test="LRT")
DMS.comp[,8][2]

# prepare result for plotting
new.Cov.DMS$methLevel <- (new.Cov.DMS$Meth/new.Cov.DMS$Cov)*100
mean <- aggregate(methLevel ~ Line, data=new.Cov.DMS, FUN="mean")
sd <- aggregate(methLevel ~ Line, data=new.Cov.DMS, FUN="sd")
summary.DMS <- data.frame(Line=c("early", "late"), mean=mean$methLevel, sd=sd$methLevel)
summary.DMS$error <- (qnorm(0.95)*summary.DMS$sd)/sqrt(nrow(new.Cov.DMS))
summary.DMS$CI_lower <- summary.DMS$mean-summary.DMS$error
summary.DMS$CI_upper <- summary.DMS$mean+summary.DMS$error

# make plot
base_col <- c("darkorange", "darkslategrey")
Plot.glmm <- summary.DMS %>% 
  ggplot(aes(x=Line, y=mean)) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), color="black", width=.2, size=1.1) +
  geom_point(aes(color=Line), shape=19, size=8, alpha=1) + 
  scale_color_manual(name="", values=base_col) +
  ylab("Methylation level (%)") + xlab("") +
  scale_y_continuous(limits=c(40,54), breaks=seq(41,53,4)) +
  theme_classic() +
  theme(axis.text=element_text(size=18),
        axis.title.y=element_text(size=20, margin=margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title.x=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 5, l = 10))
save_plot("plots/Diff.meth.lines.pdf", Plot.glmm, base_height = 4, base_width = 4)

dummy <- new.Cov.DMS
dummy$Line <- ifelse(dummy$Line==1,"early", "late")
Plot.glmm.data <- summary.DMS %>% 
  ggplot(aes(x=Line, y=mean)) +
  geom_violin(data=dummy, aes(x=Line, y=methLevel, color=Line, fill=Line), alpha=0.3) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), color="black", width=.5) +
  geom_point(aes(color=Line), shape=19, size=5, alpha=1) + 
  scale_color_manual(name="", values=base_col) +
  scale_fill_manual(name="", values=base_col) +
  ylab("Methylation level (%)") + xlab("") +
  scale_y_continuous(limits=c(0,100)) +
  theme_classic() +
  theme(axis.text=element_text(size=8),
        axis.title.y=element_text(size=9, margin=margin(t = 0, r = 4, b = 0, l = 0)),
        axis.title.x=element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = margin(t = 10, r = 10, b = 5, l = 5))
save_plot("plots/Diff.meth.lines.2.pdf", Plot.glmm.data, base_height = 2, base_width = 2)


### -------------- genome annotation -------------------------------------

# get information on length of chr/scaffolds
Chr_length <- read.table("genome/chr_length", header=F, sep="\t") 
names(Chr_length) <- c("chrom","length")

# load gff-file
gff <- makeTxDbFromGFF("genome/GCF_001522545.3_Parus_major1.1_genomic.gff", format="gff3", organism="Parus major", chrominfo=Chr_length) 

# get annotations
PROmoters_2k_new <- promoters(gff, upstream=2000, downstream=200, columns=c("tx_name", "gene_id")) 
PROmoters_2k.t <- trim(PROmoters_2k_new, use.names=TRUE)
TSS.laine <- promoters(gff, upstream=300, downstream=50, columns=c("tx_name", "gene_id")) 
TSS.laine.t <- trim(TSS.laine, use.names=TRUE)
downstream.laine <- flank(genes(gff), 10000, start=FALSE, both=FALSE, use.names=TRUE)
downstream.laine.t <- trim(downstream.laine, use.names = TRUE)
upstream.laine <- promoters(genes(gff), upstream=10000, downstream=0)
upstream.laine.t <- trim(upstream.laine, use.names = TRUE)

## export files 
export(genes(gff), "annotation/genes.gff3")
export(PROmoters_2k.t, "annotation/promoters_2k.t.gff3")
export(TSS.laine.t, "annotation/TSS.laine.t.gff3")
export(downstream.laine.t, "annotation/downstream.laine.t.gff3")
export(upstream.laine.t, "annotation/upstream.laine.t.gff3")
## save as GRanges object
Genes <- genes(gff)
save(Genes, PROmoters_2k.t, TSS.laine.t, downstream.laine.t, upstream.laine.t, file="annotation/annotations-as-GRanges.RData")

### get overlap between CpG sites and annotations in bash! ### see Genome_annotation.sh
# for which you need CpG site locations:
CpGs_all <- Dat_diffMeth[,1:2]
CpGs_all$end <- CpGs_all$start+1
write.table(CpGs_all, "temp/CpGs_all", quote=F, sep="\t", row.names=F, col.names=F)

# load overlap between CpG sites and annotations
Annotations.0 <- read.table("temp/annotated_CpGs.out", sep="\t")
names(Annotations.0) <- c("site", "chr", "pos", "annotation", "gene_id", "id", "strand_ann")

# remove duplicated genes (multiple transcripts) !!!
Annotations <- Annotations.0[!duplicated(Annotations.0[,c(1,4,5)]), ]
length(unique(Annotations$gene_id))
n_intergenic <- data.frame(Annotation="intergenic", n=length(row.names(Dat_diffMeth))-length(unique(Annotations$site)))

# summarize data:
n_annotation_temp <- Annotations %>% group_by(annotation) %>%
  summarise(
    count = n())
names(n_annotation_temp) <- names(n_intergenic)
n_annotation <- as.data.frame(rbind(n_annotation_temp, n_intergenic))
write.table(n_annotation, "final/n_annotation.txt", quote=F, sep="\t", row.names=F, col.names=T)

# Prepare plot
n_annotation$n_k <- n_annotation$n/1000
annotation_order <- n_annotation$Annotation[order(n_annotation$n)]
annotation_name <- annotation_order

# make plot 
col <- "darkslategrey"
Plot <- n_annotation %>% mutate(Annotation = factor(Annotation, levels = annotation_order)) %>% 
  ggplot() +
  geom_bar(aes(x=Annotation, y=n_k), fill=col, color="white", stat="identity", width=0.8) +
  ylab("Number of CpG sites") + xlab("") +
  scale_y_continuous(breaks=seq(100,700,200), limits=c(0,780), expand=c(0,0), labels=paste(seq(100,700,200), "k", sep="")) +
  scale_x_discrete(labels=c("transcription start site", annotation_name[2], "10k downstream", "10k upstream", annotation_name[5], "gene body")) +
  theme_classic() +
  theme(axis.text.x=element_text(size=8, angle=90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=8),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9, margin=margin(t = 0, r = 10, b = 0, l = 0)),
        axis.line = element_line(colour = "black"),
        legend.position = "none") 
save_plot("plots/Annotations_CpGs_InData.pdf", Plot, base_height = 3, base_width = 2)

# save annotation file:
diffMeth_Annotation <- merge(Dat_diffMeth[,c(8,5:7)], Annotations[,c(1,4:7)], by="site", all.x=T)
save(diffMeth_Annotation, file="final/DMA_out_Annotation_double.transcripts.removed.RData")

# get gene lists for GO analysis
# get genes with significant DMS
Genes_DMS <- diffMeth_Annotation[diffMeth_Annotation$qvalue<0.05,]
table(Genes_DMS$annotation); length(Genes_DMS$annotation[is.na(Genes_DMS$annotation)])

# exchange LOC genes for gene symbol when good prediction available:
Genes.0 <- unique(Genes_DMS$gene_id[!is.na(Genes_DMS$gene_id)])
Genes.0 <- gsub("LOC107211503", "SETDB2", Genes.0); Genes.0 <- gsub("LOC117244202", "SNORA14", Genes.0)
write.table(Genes.0, "final/Genes_DMS.txt", quote=F, col.names=F, row.names=F, sep="\n")

Genes_DMS_bf <- diffMeth_Annotation[diffMeth_Annotation$pvalue<0.05/nrow(Dat_diffMeth),]
table(Genes_DMS_bf$annotation); length(Genes_DMS_bf$annotation[is.na(Genes_DMS_bf$annotation)])

write.table(Genes_DMS_bf, "final/DMS.BF.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(unique(Genes_DMS_bf$gene_id[!is.na(Genes_DMS_bf$gene_id)]), "final/Genes_DMS_bf.txt", quote=F, col.names=F, row.names=F, sep="\n")

# prepare output files, but use readable chromosome name
# get good chr name
Chr_length_names <- read.table("genome/chr_length_names", header=F, sep="\t") 
names(Chr_length_names) <- c("Chr_name", "Chr_name_short", "chr", "length")
help <- str_split_fixed(Genes_DMS$site, "_", 3)
Genes_DMS$chr <- paste(help[,1], help[,2], sep="_")
Genes_DMS$site2 <- help[,3]
Genes_DMS_clean <- merge(Genes_DMS, Chr_length_names[,c(1,3)], by="chr", all.y=FALSE)
Genes_DMS_clean$pos <- paste(Genes_DMS_clean$Chr_name, Genes_DMS_clean$site2, sep="_")

# add gene symbol to LOC genes when good prediction available:
Genes_DMS_clean$gene_id <- gsub("LOC107211503", "LOC107211503 (SETDB2)", Genes_DMS_clean$gene_id); Genes_DMS_clean$gene_id <- gsub("LOC117244202", "LOC117244202 (SNORA14)", Genes_DMS_clean$gene_id)
Genes_DMS_clean$annotation <- gsub("gene", "gene body", Genes_DMS_clean$annotation); Genes_DMS_clean$annotation <- gsub("downstream10k", "10k downstream", Genes_DMS_clean$annotation)
Genes_DMS_clean$annotation <- gsub("upstream10k", "10k upstream", Genes_DMS_clean$annotation); Genes_DMS_clean$annotation[is.na(Genes_DMS_clean$annotation)] <- "intergenic"
write.table(Genes_DMS_clean[,c(12,3,5:7)], "final/DMS.FDR_good.chr.names.txt", quote=F, col.names=T, row.names=F, sep="\t")

# get annotations in wich significant DMS are located

# get background gene list
BackGenes <- unique(diffMeth_Annotation$gene_id[!is.na(diffMeth_Annotation$gene_id)])
BackGenes <- gsub("LOC107211503", "LOC107211503 (SETDB2)", BackGenes); BackGenes <- gsub("LOC117244202", "LOC117244202 (SNORA14)", BackGenes)
write.table(BackGenes, "final/GO_Background.txt", quote=F, col.names=F, row.names=F, sep="\n")


## get meth diff from pools
IN <- load(file="temp/MethylKit_dma_POOLS_v0.RData")
DMS_out <- NULL
for (i in 1:length(diffMeth_POOLS)) {
  dat <- diffMeth_POOLS[[i]]
  dat$site <- paste(dat$chr, dat$start, sep="_")
  dat_dms <- dat[dat$site %in% Genes_DMS_bf$site,]
  DMS_out <- rbind(DMS_out, dat_dms)
}


### -------------- Manhattan plot ----------------------------------------

# prepare CpG site location across chromosomes
Chr_length_names <- read.table("genome/chr_length_names", header=F, sep="\t") 
names(Chr_length_names) <- c("Chr_name", "Chr_name_short", "chr", "length")

order <- c("1", "1A", "2", "3", "4", "4A", as.character(5:15), as.character(17:24), "25LG1", "25LG2", "26", "27", "28", "LGE22", "Z", "Sc", "MT")
col_code <- rep(c(1,2), length(order)/2)
helper <- data.frame(order, col_code)

Chr_length_ordered <- Chr_length_names[c(1,6,2:5,7:nrow(Chr_length_names)),]
Chr_length_ordered$index <- 1:nrow(Chr_length_ordered)
Chr_length_ordered$CumSum <- cumsum(Chr_length_ordered$length)

# make plot data 
prep_Manh_temp <- merge(Dat_diffMeth, Chr_length_ordered, by="chr")
prep_Manh <- prep_Manh_temp[order(prep_Manh_temp$index, prep_Manh_temp$start),]
prep_Manh$Position <- prep_Manh$CumSum-prep_Manh$length+prep_Manh$start
prep_Manh_col <- merge(prep_Manh, helper, by.x="Chr_name_short", by.y="order")

# get chr axis label positions
# for chr + MT
Axis_lab_help_temp <- Chr_length_ordered[Chr_length_ordered$Chr_name_short!="Sc",]
Axis_lab_help_temp$axis_pos <- Axis_lab_help_temp$CumSum-(Axis_lab_help_temp$length/2)
# for scaffolds
Axis_lab_help_sc <- Chr_length_ordered[Chr_length_ordered$Chr_name_short=="Sc",]
Length_all_sc <- sum(Axis_lab_help_sc$length)
Axis_lab_help_sc2 <- Axis_lab_help_sc[nrow(Axis_lab_help_sc),]
Axis_lab_help_sc2$axis_pos <- Axis_lab_help_sc2$CumSum-(Length_all_sc/2)
# combine (MT in end)
Axis_lab_help <- rbind(Axis_lab_help_temp[-nrow(Axis_lab_help_temp),], Axis_lab_help_sc2, Axis_lab_help_temp[nrow(Axis_lab_help_temp),])

# remove MT
Data_Manh <- prep_Manh_col[prep_Manh_col$Chr_name_short!="MT",]

# prepare plot and set plot parameter
Data_Manh$col_code <- as.factor(Data_Manh$col_code)

Axis_lab_help_short <- Axis_lab_help[c(1:11,13,15,18,21,28,32,33),]
help <- Axis_lab_help$Chr_name_short
Axis_lab_help$short_lab <- c(help[1:11], "", help[13], "", help[15], "", "", help[18],
                            "", "", help[21], "", "", "", "", "", "", help[28], "", "", "", as.character(help[32:33]), "")
Axis_lab_help <- Axis_lab_help[-34,]
base_col <- c("darkslategrey", "azure4", "darkorange")

# make highlight for significant CpGs:
Data_Manh$col_code_highlights <- as.factor(ifelse(Data_Manh$qvalue<0.05, 3, ifelse(Data_Manh$qvalue>0.05 & Data_Manh$col_code==1,1,2)))

# make plot
Plot.Manhattan.DMA <- Data_Manh %>% 
  ggplot() +
  geom_hline(yintercept=sign_threshold, linetype=2, size=0.5, color="black") +
  annotate(geom="text", x=19000000, y=sign_threshold+0.3, label="BF", size=7) +
  geom_hline(yintercept=sign_threshold_fdr, linetype=2, size=0.5, color="black") +
  annotate(geom="text", x=25000000, y=sign_threshold_fdr+0.3, label="FDR", size=7) +
  geom_point(aes(x=Position, y=log10_pvalue, color=col_code_highlights, size=col_code_highlights), shape=19) +  
  scale_x_continuous(breaks=Axis_lab_help$axis_pos, expand=c(0,0), labels=Axis_lab_help$short_lab) +
  scale_y_continuous(expand=c(0,0), breaks=c(seq(2,10,length.out=5)), limits=c(0,10)) +
  scale_color_manual(name="", values=base_col) +
  scale_size_manual(name="", values=c(2,2,3)) +
  ylab(expression(paste("-log"["10"],"(p-value)"))) + xlab("Chromosome position") +
  theme_classic() +
  theme(axis.text.x=element_text(size=20, margin=margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.y=element_text(size=20),
        axis.title.x=element_text(size=22, margin=margin(t = 10, r = 0, b = 0, l = 0)),
        axis.title.y=element_text(size=22, margin=margin(t = 0, r = 7, b = 0, l = 0)),
        axis.line=element_line(colour="black"),
        plot.margin=unit(c(0.7,0.5,0.2,0.3),"cm"),
        legend.position = "none") 
save_plot("plots/manhattan_plot_Final_v1.png", Plot.Manhattan.DMA, base_height = 7, base_width = 15)
save_plot("plots/manhattan_plot_Final_v1.pdf", Plot.Manhattan.DMA, base_height = 7, base_width = 15)


