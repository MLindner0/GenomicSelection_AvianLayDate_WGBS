
# R script adjusted from: https://evodify.com/gatk-in-non-model-organism/
library(ggplot2)
library(cowplot)
library(stringr)
library(dplyr)


## ----- autosomes -----
FileLocation <- "snps"
FileNames <- paste(FileLocation, list.files(path = FileLocation, pattern = "*.table"), sep="/")

FileList <- lapply(FileNames, function(x) read.csv(x, header=T, na.strings=c("","NA"), sep="\t"))

lapply(FileList, function(x) dim(x))
lapply(FileList, function(x) head(x))
Percentiles_5_99 <- lapply(FileList, function(x) quantile(x$NA00001.DP, probs=c(0.05,0.99)))
Percentiles_f <- data.frame(Percentile_5=unlist(sapply(Percentiles_5_99, "[", 1)),
                            Percentile_99=unlist(sapply(Percentiles_5_99, "[", 2)))

sample_temp <- str_split_fixed(FileNames, "/", 2)[,2]
sample_help <- str_split_fixed(sample_temp, "\\.", 3)[,1]
Percentiles_f$Sample <- sample_help
Percentiles_f$Pool_yes <- ifelse(grepl("Pool", Percentiles_f$Sample),1,0)
write.table(Percentiles_f, "temp/Percentiles_Auto.txt", quote = F, row.names = F, col.names = T, sep = "\t")

Mean_Percentiles_f <- Percentiles_f %>% group_by(Pool_yes) %>%
  summarise(
    Mean_Percentile_99=mean(Percentile_99),
    Mean_Percentile_5=mean(Percentile_5),
    SD_Percentile_99=sd(Percentile_99),
    SD_Percentile_5=sd(Percentile_5)
    )
Mean_Percentiles_f <- as.data.frame(Mean_Percentiles_f)
write.table(Mean_Percentiles_f , "temp/Mean_Percentiles_Auto.txt", quote = F, row.names = F, col.names = T, sep = "\t")

Data_plots <- NULL
for(i in 1:length(FileList)) {
  data <- FileList[[i]]
  data$DP <- NULL
  sample <- sample_help[i]

  # remove sample name from headers
  names(data)[grep("NA00001", names(data))] <- str_split_fixed(names(data)[grep("NA00001", names(data))], "\\.", 2)[,2]
  data$SAMPLE <- rep(sample, nrow(data))
  data$Pool_yes <- ifelse(grepl("Pool", data$SAMPLE),1,0)
  Data_plots <- rbind(Data_plots, data)
}

# Make plots: here DP, GQ, QUAL
base_col <- c("darkgoldenrod1", "darkslategrey")

# All individuals
VCF_f_temp <- Data_plots[Data_plots$Pool_yes==0,]
# remove snps with DP higher than 99th percentile to aid visualization
VCF_f <- VCF_f_temp[VCF_f_temp$DP<as.numeric(Mean_Percentiles_f[Mean_Percentiles_f$Pool_yes==0,2]),]

DP <- ggplot(VCF_f, aes(x=DP)) + geom_density(alpha=.5, size=0.1, fill=base_col[2]) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank())

GQ <- ggplot(VCF_f, aes(x=GQ)) + geom_density(alpha=.5, size=0.1, fill=base_col[2]) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank())

QUAL <- ggplot(VCF_f, aes(x=QUAL)) + geom_density(alpha=.5, size=0.1, fill=base_col[2]) +
  scale_fill_manual(name="", values=col) + scale_y_continuous(expand=c(0,0)) +
  theme_classic() + theme(legend.position = "none", axis.text.y=element_blank(), axis.ticks.y=element_blank())

Plot <- plot_grid(DP, GQ, QUAL, labels = "AUTO", label_size = 10, scale = 0.9, ncol=3)
save_plot("plots/Diagnostics_Auto.pdf", Plot, base_height = 4, base_width = 9)
