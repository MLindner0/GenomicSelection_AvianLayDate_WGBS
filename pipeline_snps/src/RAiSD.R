

### RAiSD
options(width=200)
library(stringr)

setwd("/home/nioo/melaniel/projects/wgbs_snakemake_snps")

###############
####### LOAD DATA

# get file names / location
Line <- c("Early", "Late")
FileNames_Location <- NULL

# make loop to list output files // Report
for(f in 1:length(Line)){
  
  # get file names for all runs:
  FileLocation <- paste("plink/individuals_.9_QC_ok_RAiSD", Line[f], sep="_")
  FileNames = list.files(path = FileLocation, pattern = "RAiSD_Report*")
  FileNames_Location.temp <- paste(FileLocation, FileNames, sep="/")
  
  # combine runs
  FileNames_Location <- append(FileNames_Location, FileNames_Location.temp, after = length(FileNames_Location))
}

# make loop to list output files // Info
Line <- c("Early", "Late")
FileNames_Location_info <- NULL

# make loop to list files
for(f in 1:length(Line)){

  # get file names for all runs:
  FileLocation <- paste("plink/individuals_.9_QC_ok_RAiSD", Line[f], sep="_")
  FileNames = list.files(path = FileLocation, pattern = "RAiSD_Info_formatted*")
  FileNames_Location.temp <- paste(FileLocation, FileNames, sep="/")

  # combine runs
  FileNames_Location_info <- append(FileNames_Location_info, FileNames_Location.temp, after = length(FileNames_Location_info))
}

# read data
RAiDS_reports <- NULL
for(i in 1:length(FileNames_Location)) {
  if(file.size(FileNames_Location[[i]]) > 0) {
  dat <- read.table(FileNames_Location[[i]], header=FALSE, row.names=NULL, sep="\t")
  names(dat) <- c("GenomicLocation", "WindowStart", "WindowEnd", "VAR", "SFS", "LD", "mu")
  chr <- str_split_fixed(FileNames_Location[[i]], "\\.", 4)[,4]
  line.0 <- str_split_fixed(FileNames_Location[[i]], "/", 3)[,2]
  line <- str_split_fixed(line.0, "_", 6)[,6]
  dat.chr <- data.frame(SelectionLine=rep(line, nrow(dat)), Chr=rep(chr, nrow(dat)), dat)
  RAiDS_reports <- rbind(RAiDS_reports, dat.chr)
  }
}
save(RAiDS_reports, file="out/RAiDS_reports.RData")

# get number of snps considered:
RAiDS_info <- NULL
for(i in 1:length(FileNames_Location_info)) {
    dat <- read.table(FileNames_Location_info[[i]], header=FALSE, row.names=NULL, sep="\n")
    temp <- str_split_fixed(dat$V1, " ", 15)
    sites <- sum(as.numeric(temp[,7]))
    snps <- sum(as.numeric(temp[,10]))
    line.0 <- str_split_fixed(FileNames_Location_info[[i]], "/", 3)[,2]
    line <- str_split_fixed(line.0, "_", 6)[,6]
    dat.out <- data.frame(SelectionLine=line, 2, snp_type=c("all", "considered"), snp_number=c(sites, snps))
    RAiDS_info <- rbind(RAiDS_info, dat.out)
}


###############
####### FORMAT DATA

setwd("/Users/melanielindner/Documents/NIOO/Projects/ERC_molecular_data/WBGS_selection_lines/wgbs_snakemake_snps")

table(RAiDS_reports$SelectionLine, RAiDS_reports$Chr)
# 37 Chr/Scaffolds of 274 Chr/Scaffolds have windows (same Chr/Scaffolds in both lines)

# split reports by selection line (as data output will be processed separately)
line <- unique(RAiDS_reports$SelectionLine)
RAiDS_reports_split <- list()

for(i in 1: length(line)) {
  RAiDS_reports_split[[i]] <- RAiDS_reports[RAiDS_reports$SelectionLine==line[i],]
}

lapply(RAiDS_reports_split, dim) # more windows for late selection line (E:238304, L:261924)

# define threshold for significance (top 0.05%)
mu_sign_threshold <- lapply(RAiDS_reports_split, function(x) quantile(x$mu, probs=(1-0.0005)))
# now many regions?
n_regions_sign_threshold <- lapply(RAiDS_reports_split, function(x) nrow(x)*0.0005)

# get significance label
RAiDS.out <- list()
for(i in 1:length(RAiDS_reports_split)) {
  dat <- RAiDS_reports_split[[i]]
  dat$sign <- ifelse(dat$mu>as.numeric(mu_sign_threshold[i]),1,0)
  RAiDS.out[[i]] <- dat
}

# filter for "significant" results 
RAiDS.out_sign <- lapply(RAiDS.out, function(x) x[x$sign==1,])
lapply(RAiDS.out_sign, function(x) table(x$Chr))



###############
####### Evaluate DATA

### Make "Manhattan plots"

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
help <- RAiDS_reports[!duplicated(RAiDS_reports[,2:3]),2:3]; names(help) <- c("chr", "pos") # all genomic locations for both lines
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
# combine
axis.labs <- rbind(axis.help.0, axis.help.2)

# save snp positions
save(snps.positions, axis.labs, file="out/RAiDS_for.Manahttan.RData")

# prepare plot and set plot parameter 
dat.early <- RAiDS.out[[1]]
dat.plot.e <- merge(dat.early, snps.positions[,c(1,2,6,9)], by.x=c("Chr", "GenomicLocation"), by.y=c("chr", "pos"))
dat.late <- RAiDS.out[[2]]
dat.plot.l <- merge(dat.late, snps.positions[,c(1,2,6,9)], by.x=c("Chr", "GenomicLocation"), by.y=c("chr", "pos"))

row.names(axis.labs) <- 1:nrow(axis.labs)
help <- axis.labs$Chr_name_short
axis.labs$x.lab <- c(help[1:10], "", help[12], rep("", 4), help[17], rep("", 7), help[25], rep("", 6), help[32])

help.x.lim <- unique(c(dat.plot.e$genome.pos, dat.plot.l$genome.pos))
x.lim <- c(range(help.x.lim)[1]-range(help.x.lim)[1], range(help.x.lim)[2]+range(help.x.lim)[1])

range(c(dat.plot.e$mu, dat.plot.l$mu))

# make plot: 
pdf(file = "plots/RAiDS.pdf", width = 9, height = 7)
par(mfrow = c(2,1))
par(mar = c(2,4,0.5,2), mgp=c(2, 0.5, 0), xaxs = "i", yaxs = "i")

plot(NULL, xlim=x.lim, ylim=c(0,8.5), xlab="",
     ylab=expression(mu), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,8,2), labels=seq(0,8,2), lwd = 0, lwd.ticks = 1)
text(40000000, 8.3, "Early", cex=1, pos=1, offset=0, font=2)
non_sig <- dat.plot.e[dat.plot.e$sign==0,]
points(non_sig$genome.pos, non_sig$mu, pch=19, cex=0.5, col=ifelse(non_sig$col==1, "darkslategrey", "azure4"))
sig <- dat.plot.e[dat.plot.e$sign==1,]
points(sig$genome.pos, sig$mu, pch=19, cex=1, col="darkorange")

plot(NULL, xlim=x.lim, ylim=c(0,8.5), xlab="",
     ylab=expression(mu), xaxt="n", yaxt="n", bty = "l")
axis(1, at=axis.labs$lab.pos, labels=axis.labs$Chr_name_short, lwd = 0, lwd.ticks = 1)
axis(2, at=seq(0,8,2), labels=seq(0,8,2), lwd = 0, lwd.ticks = 1)
text(40000000, 8.3, "Late", cex=1, pos=1, offset=0, font=2)
non_sig <- dat.plot.l[dat.plot.l$sign==0,]
points(non_sig$genome.pos, non_sig$mu, pch=19, cex=0.5, col=ifelse(non_sig$col==1, "darkslategrey", "azure4"))
sig <- dat.plot.l[dat.plot.l$sign==1,]
points(sig$genome.pos, sig$mu, pch=19, cex=1, col="darkorange")

dev.off()


## Get genome annotations
# prepare file for merging with annotations:

for_anotation <- RAiDS_reports[!duplicated(RAiDS_reports[,2:5]),c(2,4,5,3)]# all genomic locations for both lines
write.table(for_anotation, "temp/RAiDS_regions.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

# do intersection on hpc
annotations.RAiDS <- read.table("temp/annotated_RAiDS_regions")
annotations.RAiDS$V3 <- NULL
annotations.RAiDS$V4 <- NULL
annotations.RAiDS$V9 <- NULL
annotations.RAiDS$V11 <- NULL
names(annotations.RAiDS) <- c("Chr", "GenomicLocation", "AnnStart", "AnnStop", "Strand", "GeneID", "Feature")

annotations.RAiDS <- annotations.RAiDS[!duplicated(annotations.RAiDS),]
background.list.GO.0 <- unique(annotations.RAiDS$GeneID)

RAiDS.out_sign_annotated <- NULL
for(i in 1:length(RAiDS.out_sign)) {
  dat <- RAiDS.out_sign[[i]]
  annotated <- merge(RAiDS.out_sign[[i]], annotations.RAiDS, by=c("Chr", "GenomicLocation"), all.x=TRUE)
  annotated_goodchrname <- merge(Chr_length_names[,c(1,3)], annotated, by.x="chr", by.y="Chr")
  RAiDS.out_sign_annotated <- rbind(RAiDS.out_sign_annotated, annotated_goodchrname)
} 

LOC <- unique(RAiDS.out_sign_annotated$GeneID[grepl("LOC", RAiDS.out_sign_annotated$GeneID)])
write.table(LOC, "temp/LOCs_RAiDS.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
# read LOC replacements
LOC_replacements <- read.table("temp/LOC_replacements.txt", header=FALSE, sep="\t")
names(LOC_replacements) <- c("Loc", "GeneID")
LOC_replacements$LocAndGeneID <- paste(LOC_replacements$Loc, " (", LOC_replacements$GeneID, ")", sep="")

# get GO lists
# 1 BG
background.list.GO <- background.list.GO.0[!is.na(background.list.GO.0)]
for(i in 1:nrow(LOC_replacements)) {
  background.list.GO <- gsub(LOC_replacements$Loc[i], LOC_replacements$GeneID[i], background.list.GO)
}
write.table(sort(background.list.GO), "temp/RAiDS.background.list.GO.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
length(background.list.GO) # 15307

# 2 early
genes.early.0 <- unique(RAiDS.out_sign_annotated[RAiDS.out_sign_annotated$SelectionLine=="Early", "Chr_name"])
genes.early <- genes.early.0[!is.na(genes.early.0)]
for(i in 1:nrow(LOC_replacements)) {
  genes.early <- gsub(LOC_replacements$Loc[i], LOC_replacements$GeneID[i], genes.early)
}
write.table(sort(genes.early), "temp/RAiDS.early.GO.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
length(genes.early) # 44

genes.late.0 <- unique(RAiDS.out_sign_annotated[RAiDS.out_sign_annotated$SelectionLine=="Late", "GeneID"])
genes.late <- genes.late.0[!is.na(genes.late.0)]
for(i in 1:nrow(LOC_replacements)) {
  genes.late <- gsub(LOC_replacements$Loc[i], LOC_replacements$GeneID[i], genes.late)
}
write.table(sort(genes.late), "temp/RAiDS.late.GO.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\n")
length(genes.late) # 60

genes.early[genes.early %in% genes.late]
genes.late[genes.late %in% genes.early]

# get table for supplementary material:
RAiDS.supp <- RAiDS.out_sign_annotated[,c(4,2,3,5:10,16,15)]
for(i in 1:nrow(LOC_replacements)) {
  RAiDS.supp$GeneID <- gsub(LOC_replacements$Loc[i], LOC_replacements$LocAndGeneID[i], RAiDS.supp$GeneID)
}
RAiDS.supp$Feature[is.na(RAiDS.supp$Feature)] <- "intergenic"
RAiDS.supp$Feature[RAiDS.supp$Feature=="tss"] <- "transcription start site"
RAiDS.supp$Feature[RAiDS.supp$Feature=="gene"] <- "gene body"
RAiDS.supp$Feature[RAiDS.supp$Feature=="downstream10k"] <- "10k downstream"
RAiDS.supp$Feature[RAiDS.supp$Feature=="upstream10k"] <- "10k upstream"

write.table(RAiDS.supp, "out/RAiDS.supp.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")






