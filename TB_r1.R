
## Sulfur gene detection
"1133136", "1132874", "1132770"

remove_SAGs <- c("1132750", "1132751", "1132752", "1132753", "1132754", "1132755", "1132756", "1132758", "1132759", "1132760", "1132762", "1132764", "1132766", "1132767", "1132768", "1132769", "1132771", "1132772", "1132870", "1132871", "1132878", "1132879", "1132880", "1132881", "1132882", "1132883", "1132887", "1132889", "1132891", "1132892", "1133008", "1133010", "1133012", "1133015", "1133016", "1133017", "1133019", "1133024", "1133025", "1133026", "1133027", "1133028", "1133127", "1133128", "1133130", "1133131", "1133132", "1133133", "1133134", "1133135", "1133141", "1133142", "1133144", "1133146")

metadata <- read.table("Bin1_coverages_R.txt", header = TRUE, row.names = 1)[, c("sample", "location", "date", "fraction")]
metadata <- metadata[order(metadata$sample),]

sulfur_detect <- read.table("IMG_sulfur_genes-DETECTION.txt", header = TRUE, row.names = 1)
gene_calls_sulfur <- read.table("genecalls.txt", header = TRUE, row.names = 1)
gene_calls_sulfur <- gene_calls_sulfur[rownames(sulfur_detect),]
colnames(sulfur_detect) <- gsub("a", "", colnames(sulfur_detect))
rownames(sulfur_detect) <- paste(gene_calls_sulfur$contig, rownames(gene_calls_sulfur), sep = "_")
sulfur_detect <- sulfur_detect[,setdiff(colnames(sulfur_detect), remove_SAGs)]
sulfur_detect <- t(sulfur_detect)
metadata1 <- metadata[rownames(sulfur_detect),]
sulfur_detect <- cbind(sulfur_detect, metadata1)
sulfur_detect_pos_m <- data_clean2(melt(remove_controls(sulfur_detect)))
sulfur_detect_pos_m$date <- factor(sulfur_detect_pos_m$date)
sulfur_detect_pos_m$variable <- factor(sulfur_detect_pos_m$variable, levels = unique(sort(as.character(sulfur_detect_pos_m$variable))))
sulfur_plot_data <- aggregate(value ~ date + variable, data = sulfur_detect_pos_m, FUN= "median")

sulfur_plot_data$Bin <- str_split_fixed(sulfur_plot_data$variable, "_", 5)[,4]
sulfur_plot_data$Bin <- factor(sulfur_plot_data$Bin, levels = c("00003", "00002", "00001", "00005", "00004"))

p <- ggplot(sulfur_plot_data, aes(x = date, y = variable)) + 
  geom_tile(aes(fill = value)) + scale_fill_viridis(name="Proportion detection") +
  ylab("gene") +
  ggtitle("Detection of sulfur genes") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1))

p + geom_hline(yintercept = (c(max(grep("Bin_00004", levels(sulfur_plot_data$variable))),
                                  max(grep("Bin_00005", levels(sulfur_plot_data$variable))),
                                  max(grep("MAG_00001", levels(sulfur_plot_data$variable))),
                                  max(grep("MAG_00002", levels(sulfur_plot_data$variable))),
                                  max(grep("MAG_00003", levels(sulfur_plot_data$variable))))+0.5), col = "red") +
  geom_vline(xintercept = (c(max(grep("2005", levels(sulfur_plot_data$date))),
                             max(grep("2007", levels(sulfur_plot_data$date))),
                             max(grep("2008", levels(sulfur_plot_data$date))),
                             max(grep("2009", levels(sulfur_plot_data$date))),
                             max(grep("2010", levels(sulfur_plot_data$date))),
                             max(grep("2012", levels(sulfur_plot_data$date))),
                             max(grep("2013", levels(sulfur_plot_data$date))))+0.5), col = "red")



sulfur_cov <- read.table("IMG_sulfur_genes-COVERAGES.txt", header = TRUE, row.names = 1)
gene_calls_sulfurcov <- read.table("genecalls.txt", header = TRUE, row.names = 1)
gene_calls_sulfurcov <- gene_calls_sulfur[rownames(sulfur_cov),]
colnames(sulfur_cov) <- gsub("a", "", colnames(sulfur_cov))
rownames(sulfur_cov) <- paste(gene_calls_sulfurcov$contig, rownames(gene_calls_sulfurcov), sep = "_")
sulfur_cov <- sulfur_cov[,setdiff(colnames(sulfur_cov), remove_SAGs)]
sulfur_cov <- t(sulfur_cov)
metadata1 <- metadata[rownames(sulfur_cov),]
sulfur_cov <- cbind(sulfur_cov, metadata1)
num_reads2 <- data.frame(num_reads[rownames(remove_controls(sulfur_cov)),])
sulfur_cov <- remove_controls(sulfur_cov)
sulfur_cov$readCount <- factor(num_reads2$num_reads.rownames.remove_controls.sulfur_cov.....)
sulfur_cov_pos_m <- data_clean2(melt(sulfur_cov))
sulfur_cov_pos_m$date <- factor(sulfur_cov_pos_m$date)
sulfur_cov_pos_m$variable <- factor(sulfur_cov_pos_m$variable, levels = unique(sort(as.character(sulfur_cov_pos_m$variable))))

sulfur_cov_pos_m$norm_value <- (sulfur_cov_pos_m$value/as.numeric(sulfur_cov_pos_m$readCount))*1000000

sulfur_plot_datacov <- aggregate(norm_value ~ date + variable, data = sulfur_cov_pos_m, FUN= "median")

p <- ggplot(sulfur_plot_datacov, aes(x = date, y = variable)) + 
  geom_tile(aes(fill = log10(sulfur_plot_datacov$norm_value+1))) + 
  scale_fill_gradientn(name="Median normalized coverage", colors = c(viridis_pal(begin = 0, end = 0.3)(12),
                                                                      viridis_pal(begin = 0.3, end = 0.5)(5),
                                                                      viridis_pal(begin = 0.5, end = 1)(15))) +
  ylab("gene") +
  ggtitle("Coverage of sulfur genes") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1))

p + geom_hline(yintercept = (c(max(grep("Bin_00004", levels(sulfur_plot_data$variable))),
                               max(grep("Bin_00005", levels(sulfur_plot_data$variable))),
                               max(grep("MAG_00001", levels(sulfur_plot_data$variable))),
                               max(grep("MAG_00002", levels(sulfur_plot_data$variable))),
                               max(grep("MAG_00003", levels(sulfur_plot_data$variable))))+0.5), col = "red") +
  geom_vline(xintercept = (c(max(grep("2005", levels(sulfur_plot_data$date))),
                             max(grep("2007", levels(sulfur_plot_data$date))),
                             max(grep("2008", levels(sulfur_plot_data$date))),
                             max(grep("2009", levels(sulfur_plot_data$date))),
                             max(grep("2010", levels(sulfur_plot_data$date))),
                             max(grep("2012", levels(sulfur_plot_data$date))),
                             max(grep("2013", levels(sulfur_plot_data$date))))+0.5), col = "red")

sulfur_KEGG <- read.table("sulfur_kegg.txt", header = FALSE, row.names = 1)
rownames(sulfur_KEGG) <- colnames(sulfur_detect)[1:72]

```