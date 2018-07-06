
viral_bin_files <- c("viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_10/Bin_10-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_14/Bin_14-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_18/Bin_18-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_11/Bin_11-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_15/Bin_15-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_8/Bin_8-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_12/Bin_12-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_16/Bin_16-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_9/Bin_9-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_13/Bin_13-gene_coverages.txt", 
                     "viral_solo_summary/viral_solo_summary/bin_by_bin/Bin_17/Bin_17-gene_coverages.txt")
metadata <- read.table("Bin1_coverages_R.txt", header = TRUE, row.names = 1)[, c("sample", "location", "date", "fraction")]
metadata <- metadata[order(metadata$sample),]

MAG1_15_genes_bin11 <- as.character(c(1955, 1954, 1953, 1952, 1951, 1950, 1949, 1948, 1947, 1946, 1945, 1944, 1943, 1942, 1941, 1940, 1939, 1938, 1937, 1936, 1935, 1934))
Bin4_22_genes_bin6 <- as.character(c(219, 218, 217, 216, 215, 214, 213, 212, 211, 210, 209, 208, 207, 206, 205))
MAG2_83_genes_bin9 <- as.character(c(5963, 5962, 5961, 5960, 5959, 5958, 5957, 5956, 5955, 5954, 5953, 5952, 5951, 5950, 5949, 5948, 5947, 5946, 5945, 5944, 5943, 5942, 5941, 5940, 5939, 5938, 5937, 5936, 5935, 5934, 5933, 5932, 5931, 5930, 5929, 5928, 5927, 5926, 5925, 5924, 5923, 5922, 5921, 5920, 5919, 5918, 5917))

####################
### USING SCMG DOESN'T SEEM TO WORK SINCE MOST READS DON'T MAP TO MY CONTIGS
all_gene_cov <- read.table("/Users/mberg/gene_cov_detect-GENE-COVERAGES.txt", header = TRUE, row.names = 1)
colnames(all_gene_cov) <- gsub("a", "", colnames(all_gene_cov))

SCMG_IDs <- read.table("/Users/mberg/SCMG_IDs.txt", header = FALSE)

all_gene_cov_SCMG <- all_gene_cov[SCMG_IDs$V1,]
keep <- rownames(subset(metadata, fraction == "Hypolimnion"))
all_gene_cov_SCMG_noMDA <- all_gene_cov_SCMG[,keep]
### remove samples that have 0% or greater zeros in the SCMGs
all_gene_cov_SCMG_noMDA_cleaned <- all_gene_cov_SCMG_noMDA[, apply(all_gene_cov_SCMG_noMDA, 2, function(x) sum(x == 0)) < 174]
all_gene_cov_SCMG_noMDA_cleaned_ave <- apply(all_gene_cov_SCMG_noMDA_cleaned, 2, function(x) mean(x))

viral_bin_cov_6_clean <- viral_bin_cov_6[names(all_gene_cov_SCMG_noMDA_cleaned_ave), ]
vBin6_cov_norm <- apply(viral_bin_cov_6_clean, 2, function(x) x/(all_gene_cov_SCMG_noMDA_cleaned_ave))
vBin6_cov_norm <- vBin6_cov_norm[complete.cases(vBin6_cov_norm), ]
vBin6_cov_norm_cleaned <- vBin6_cov_norm[rowSums(vBin6_cov_norm != 0) > 0, ]
####################

####################
## USE TOTAL READ COUNT (FROM SEQUENCER) - it's better?
num_reads <- read.table("/Users/mberg/num_reads.txt", header = TRUE, row.names = 1)
#num_reads <- data.frame(num_reads[rownames(viral_bin_cov_6),])
#rownames(num_reads) <- rownames(viral_bin_cov_6)
#colnames(num_reads) <- c("num_reads")
#vBin6_cov_norm <- apply(viral_bin_cov_6, 2, function(x) x/(num_reads$num_reads))
#vBin6_cov_norm_cleaned <- vBin6_cov_norm[rowSums(vBin6_cov_norm != 0) > 0, ]
####################


##### FUNCTION ###############
## remove outliers; top and bottom 10% from each sample*date
######################
data_clean2 <- function(x) {
  mDate <- format(as.Date(x$date, format = "%m/%d/%Y"), "20%y/%m/%d")
  x$date <- as.Date(mDate)
  NCx <- ncol(x)
  x[,(NCx+1):(NCx+3)] <- str_split_fixed(x$date, "-", 3)
  colnames(x)[(NCx+1):(NCx+3)] <- c("year", "month", "day")
  
  spring <- c("03", "04", "05")
  summer <- c("06", "07", "08")
  fall <- c("09", "10", "11")
  
  x$season <- x$month
  x[x$season %in% spring, "season"] <- "spring"
  x[x$season %in% summer, "season"] <- "summer"
  x[x$season %in% fall, "season"] <- "fall"
  
  x$seasonYear <- paste(x$season, x$year, sep ="")
  x$seasonYear <- factor(x$seasonYear, levels = c("summer2005", "spring2007", "summer2007", "fall2007",
                                                  "spring2008", "summer2008", "fall2008",
                                                  "spring2009", "summer2009", "summer2010",
                                                  "summer2012", "summer2013", "summer2017", 
                                                  "fall2017"))
  return(x)
}
##############################

for(i in 1:length(viral_bin_files)){
  temp <- t(read.table(viral_bin_files[i], header = TRUE, row.names = 1))
  colnames(temp) <- paste(colnames(temp), "a", sep="")
  rownames(temp) <- gsub("a", "", rownames(temp))
  keep <- rownames(subset(metadata, fraction == "Hypolimnion"))
  temp <- temp[keep,]
  assign(paste("viral_bin_cov", i, sep="_"), temp)
  num_reads2 <- data.frame(num_reads[rownames(temp),])
  rownames(num_reads2) <- rownames(temp)
  temp_norm <- apply(temp, 2, function(x) x/(num_reads2$num_reads))
  temp_norm_cleaned <- temp_norm[rowSums(temp_norm) > 0, ]
  if(class(temp_norm_cleaned) == "numeric"){
    print(paste("Too many zeros in", i, sep =" "))
  }else{
    temp_norm_cleaned <- temp_norm_cleaned[, colSums(temp_norm_cleaned) > 0]
    assign(paste("viral_bin_cov_norm_cleaned", i, sep="_"), temp_norm_cleaned)
  }
}

# [1] "Too many zeros in 5"
# [1] "Too many zeros in 7"
# [1] "Too many zeros in 8"

colnames(viral_bin_cov_norm_cleaned_6) <- gsub("a", "", colnames(viral_bin_cov_norm_cleaned_6))
viral_bin_cov_norm_cleaned_6[, Bin4_22_genes_bin6]
colnames(viral_bin_cov_norm_cleaned_9) <- gsub("a", "", colnames(viral_bin_cov_norm_cleaned_9))
colnames(viral_bin_cov_norm_cleaned_11) <- gsub("a", "", colnames(viral_bin_cov_norm_cleaned_11))

cor_matrix = cor(data.matrix(viral_bin_cov_norm_cleaned_1), use = "pairwise")
cor_matrix = cor(data.matrix(viral_bin_cov_norm_cleaned_2), use = "pairwise")
cor_matrix = cor(data.matrix(viral_bin_cov_norm_cleaned_3), use = "pairwise")
cor_matrix = cor(data.matrix(viral_bin_cov_norm_cleaned_4), use = "pairwise")
cor_matrix = cor(data.matrix(viral_bin_cov_norm_cleaned_5), use = "pairwise")
cor_matrix = cor(data.matrix(viral_bin_cov_norm_cleaned_6[, Bin4_22_genes_bin6]), use = "pairwise")

pheatmap(
  mat               = cor_matrix,
  color             = brewer.pal(9, "RdBu"),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  #annotation_col    = data.frame(ID = colnames(viral_bin_cov_norm_cleaned_6[, Bin4_22_genes_bin6]), row.names = colnames(viral_bin_cov_norm_cleaned_6[, Bin4_22_genes_bin6])),
  #annotation_row    = data.frame(ID = colnames(viral_bin_cov_norm_cleaned_6[, Bin4_22_genes_bin6]), row.names = colnames(viral_bin_cov_norm_cleaned_6[, Bin4_22_genes_bin6])),
  #annotation_colors = mat_colors,
  #cluster_rows = FALSE,
  #cluster_cols = FALSE,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Correlation within contig",
  breaks = rev(seq(1, -1, length.out = 10))
)

#cbind(viral_bin_cov_norm_cleaned_1, metadata[rownames(viral_bin_cov_norm_cleaned_1),])
p <- ggplot(melt(viral_bin_cov_norm_cleaned_1 + 0.0000001), 
            aes(x = Var2, y = value, group = Var1)) + geom_line() + scale_y_log10()
  
  geom_point(size = 0.8) +
  ylab("Coverage") +
  ggtitle(paste("Contig coverage w/ median coverage line", m, sep = "-")) +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1)) +
  #scale_color_brewer(palette = "Set2")
  scale_color_manual(values = c("#fdae61", "#ca0020", "#2c7bb6")) +
  scale_shape_manual(values = c(1, 16, 16))



ggplot(aggregate(value ~ date + variable, data = viral_bin_cov_1, FUN= "median"), aes(x = date, y = variable)) + 
  geom_tile(aes(fill = value)) + scale_fill_viridis(name="Proportion detection") +
  ylab("gene") +
  ggtitle("Detection of viral genes - ChPeak_MDA_Bin_00004_000000000020") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 45, hjust = 1))


