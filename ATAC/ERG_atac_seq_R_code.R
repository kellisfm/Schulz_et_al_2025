### ERG ATAC-seq analysis


### figure 1
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("Rsubread", quietly = TRUE)){
  BiocManager::install("Rsubread")
}

if (!require("DESeq2", quietly = TRUE)){
  BiocManager::install("DESeq2")
}

# libraries
library(DESeq2)
library(Rsubread)
library(stringr)


BAM1 <- "~/mnt/largeprojects/Kai/collab/ERG/ATAC/WL2801_none_na_TeloHAEC_none_hsap_NA_WT-ERG_TCAG1_R_roughuniq_sortedByPosition.bam"
BAM2 <- "~/mnt/largeprojects/Kai/collab/ERG/ATAC/WL2803_none_na_TeloHAEC_none_hsap_NA_Del9-ERG-Clone_TCAG1_R_roughuniq_sortedByPosition.bam"
BAM3 <- "~/mnt/largeprojects/Kai/collab/ERG/ATAC/WL2802_none_na_TeloHAEC_none_hsap_NA_WT-ERG_TCAG1_R_roughuniq_sortedByPosition.bam"
BAM4 <- "~/mnt/largeprojects/Kai/collab/ERG/ATAC/WL2804_none_na_TeloHAEC_none_hsap_NA_Del15-ERG_TCAG1_R_roughuniq_sortedByPosition.bam"

saf.file <- "~/TOBIAS/TOBIAS_tmp_ERG_run/ATAC/merged_peaks_chronic.saf"
fc_PE = featureCounts(files=c(BAM1, BAM2,BAM3,BAM4,), annot.ext=saf.file, isPairedEnd=T)

saveRDS(fc_PE, "~/mnt/largeprojects/Kai/collab/ERG/ATAC/genrich_peaks/KO_wt_chronic_pseudo_diffbind.RDS")

fc_PE <- readRDS("~/mnt/largeprojects/Kai/collab/ERG/ATAC/genrich_peaks/KO_wt_chronic_pseudo_diffbind.RDS")
# read in sample metadata. 
meta.data <- read.csv("~/mnt/largeprojects/Kai/collab/ERG/ATAC/diffbind_guidesheet_CRISPR.csv",row.names = 1)

all(colnames(fc_PE$counts)==rownames(meta.data))
colnames(fc_PE$counts)==rownames(meta.data)

dds <- DESeqDataSetFromMatrix(countData = fc_PE$counts,
                              colData = meta.data,
                              design = ~ Treatment)
dds <- DESeq(dds)
res <- results(dds)
resName <- resultsNames(dds)


res <- results(dds, contrast=c("Treatment","WT","ERGKO"))
res[c('chromosome', 'peak_start',"peak_end")] <- str_split_fixed(rownames(res24), "\\.",3)
write.csv(res, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/deseq2_results_ERGKO.csv")



### figure 3
## differential accessibility
# set the working directory 
setwd("/home/kellis/TOBIAS/TOBIAS_tmp_ERG_run")

## generate variables for the path to each of your bams (not merged), along with the path to the merged saf file you made.
### note: Feature counts can read in any number of BAMS, so feel free to add more BAM variables as needed.
BAM1 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3519_none_aortic_HAEC_none_hsap_4_72hr-ERG-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM2 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3518_none_aortic_HAEC_none_hsap_4_24hr-ERG-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM3 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3521_none_aortic_HAEC_none_hsap_4_72hr-Ctrl-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM4 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3520_none_aortic_HAEC_none_hsap_4_24hr-Ctrl-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM5 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3515_none_aortic_HAEC_none_hsap_2_72hr-ERG-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM6 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3514_none_aortic_HAEC_none_hsap_2_24hr-ERG-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM7 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3517_none_aortic_HAEC_none_hsap_2_72hr-Ctrl-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"
BAM8 <- "/hpf/largeprojects/mdwilson/Kai/collab/ERG/ATAC/2024_timecourse/WL3516_none_aortic_HAEC_none_hsap_2_24hr-Ctrl-SiRNA_TCAG1_R_roughuniq_sortedByPosition_dedup.bam"

saf.file <- "~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/inputs/Peaks/merged_peaks.saf"

# run feature counts to read all your counts into R. The fc_PE variable will be a list of 4 objects - fc_PE$counts is the one we will use in the future
fc_PE = featureCounts(files=c(BAM1, BAM2,BAM3,BAM4,BAM5, BAM6,BAM7,BAM8), annot.ext=saf.file, isPairedEnd=T)
saveRDS(fc_PE, "~/mnt/largeprojects/Kai/collab/ERG/ATAC/genrich_peaks/KO_wt_chronic_pseudo_diffbind.RDS")

fc_PE <- readRDS("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/genrich_peaks/counts_matrix_for_deseq2.RDS")
# read in sample metadata. 
meta.data <- read.csv("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/diffbind_guidesheet_Timecourse_2_4_only.csv",row.names = 1)

### make sure that the colnames of your counts matrix have the exact same sample names and order as your metadata rownames.
# it is critical that the col and rownames are identical and in the same order. DEseq will not be happy if they are not
all(colnames(fc_PE$counts)==rownames(meta.data))
colnames(fc_PE$counts)==rownames(meta.data)

# using your metadata and counts matrix run deseq
### NOTE: you can easily use any diffrential expression package for this step. I am most familiar with DESeq, so that is what I have used here
# this paper seems to imply that for low quality samples RUVseq or edger may perform better: https://www.nature.com/articles/s41598-020-66998-4 

dds <- DESeqDataSetFromMatrix(countData = fc_PE$counts,
                              colData = meta.data,
                              design = ~ Treatment)
dds <- DESeq(dds)
res <- results(dds)
resName <- resultsNames(dds)


res24 <- results(dds, contrast=c("Treatment","ERG_24","Ctrl_24"))
res24[c('chromosome', 'peak_start',"peak_end")] <- str_split_fixed(rownames(res24), "\\.",3)
res72 <- results(dds, contrast=c("Treatment","ERG_72","Ctrl_72"))
res72[c('chromosome', 'peak_start',"peak_end")] <- str_split_fixed(rownames(res72), "\\.",3)
write.csv(res24, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/deseq2_results_24h.csv")
write.csv(res72, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/deseq2_results_72h.csv")


library(ggplot2)

res24 <- read.csv("/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/deseq2_results_24h.csv")


res24$diffexpressed <- "Non-DE"
res24$diffexpressed[res24$log2FoldChange>0.1 & res24$padj<0.05 ] <- "UP"
res24$diffexpressed[res24$log2FoldChange<0.1 & res24$padj<0.05 ] <- "DOWN"
res24$keep <- sample(0:7, length(res24$diffexpressed), replace = T)==1 
res24$keep[res24$diffexpressed %in%  c("UP", "DOWN")] <- T

res24 <- res24[res24$keep,]
library(gridExtra)
p24 <- ggplot(res24, aes(x = log2FoldChange, y = -log10(padj), colour =  diffexpressed)) +
  geom_point(, size = 5, alpha = 0.8) + 
  scale_color_manual(values=c("skyblue", "grey", "salmon")) +
  labs(
    title = "ErgKO 24 hours vs Ctrl siRNA",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") + # Significance threshold
  xlim(-7.5,7.5) + 
  ylim(0,150)
library(gg)
#geom_vline(xintercept = c(1), linetype = "dashed", color = "salmon")  
ggsave(p24, filename="/home/kellis/wilson_lab_documents/Figures/ERG_project/Figure_3/uniformly_scaled_volcano_24.pdf")
res72 <- read.csv("/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/deseq2_results_72h.csv")

res72$diffexpressed <- "Non-DE"
res72$diffexpressed[res72$log2FoldChange>0.1 & res72$padj<0.05 ] <- "UP"
res72$diffexpressed[res72$log2FoldChange<0.1 & res72$padj<0.05 ] <- "DOWN"

res72$keep <- sample(0:7, length(res72$diffexpressed), replace = T)==1 
res72$keep[res72$diffexpressed %in%  c("UP", "DOWN")] <- T

res72 <- res72[res72$keep,]


p72 <- ggplot(res72, aes(x = log2FoldChange, y = -log10(padj), colour =  diffexpressed)) +
  geom_point(, size = 5, alpha = 0.8) + 
  scale_color_manual(values=c("#00a7fffb", "grey", "#ff2a00ff")) +
  labs(
    title = "ErgKO 72 hours vs Ctrl siRNA",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  )+
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlim(-7.5,7.5)+ 
  ylim(0,150)

ggsave(p72, filename="/home/kellis/wilson_lab_documents/Figures/ERG_project/Figure_3/uniformly_scaled_volcano_72.pdf")

sum(na.omit(res72[,"padj"]) < 0.05)
sum(na.omit(res24[,"padj"]) < 0.05)

vsd <- vst(dds)
plotPCA(vsd, intgroup = "Treatment")

if (!require("stringr", quietly = TRUE)){
  install.packages("stringr")
}
library(stringr)
res[c('chromosome', 'peak_start',"peak_end")] <- str_split_fixed(rownames(res), "\\.",3)

## save as a CSV or an RDS for future analysis
write.csv(res, "~/mnt/largeprojects/Kai/collab/ERG/ATAC/genrich_peaks/KO_wt_chronic_pseudo_diffbind.csv")
saveRDS(dds, "~/mnt/largeprojects/Kai/collab/ERG/ATAC/genrich_peaks/KO_wt_chronic_pseudo_diffbind.RDS")


## dotplot
library(tidyverse)

preprocess_GREAT_File <- function(great_global_export_path, row_skip = 3 , condition_name=""){
  name <- gsub(".*/(.*)\\..*","\\1",great_global_export_path)
  gprofiler_output_red <- vroom::vroom(great_global_export_path, skip = row_skip)
  gprofiler_output_red <- gprofiler_output_red[,c(1,3,7,8)]
  
  colnames(gprofiler_output_red) <- c("source","term_name","negative_log10_of_adjusted_p_value", "Fold_enrich")
  colnames(gprofiler_output_red)[c(3,4)] <- paste(colnames(gprofiler_output_red)[c(3,4)],condition_name,sep="-")
  gprofiler_output_red[,3] <- -log10(gprofiler_output_red[,3])
  assign(name,gprofiler_output_red, envir = .GlobalEnv)
  
  
}

### preprocess all input files from GREAT (currently set up to use global export), merge into a final dataset
preprocess_GREAT_File("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/ROSE/Diffbind_rose/GREAT/72hr_super_down_ERGko.tsv") #, skip = 3
preprocess_GREAT_File("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/ROSE/Diffbind_rose/GREAT/Ch_super_down_ergko.tsv") #, skip = 3

preprocess_GREAT_File("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/ROSE/Diffbind_rose/GREAT/24hr_super_up_ergko.tsv") #, skip = 3
preprocess_GREAT_File("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/ROSE/Diffbind_rose/GREAT/72hr_super_up_ergko.tsv") #, skip = 3
preprocess_GREAT_File("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/ROSE/Diffbind_rose/GREAT/Ch_super_up_ergko.tsv") #, skip = 3

tc_downreg <- inner_join(`24h_KO_downreg`,`72h_KO_downreg`, by = c("source","term_name"))
all_downreg <- inner_join(tc_downreg,Chronic_down, by = c("source","term_name")) 


tc_upreg <- inner_join(`24h_up`,`72h_up`, by = c("source","term_name"))
all_upreg <- inner_join(tc_upreg,Chronic_up, by = c("source","term_name")) 


all_downreg <- inner_join(`72h_KO_downreg`,Chronic_down, by = c("source","term_name"))
all_downreg <- inner_join(tc_downreg,Chronic_down, by = c("source","term_name")) 


tc_upreg <- inner_join(`24h_up`,`72h_up`, by = c("source","term_name"))
all_upreg <- inner_join(tc_upreg,Chronic_up, by = c("source","term_name")) 


# Convert to tall format
tall_up_data <- all_upreg %>%
  pivot_longer(
    cols = starts_with("negative_log10_of_adjusted_p_value") | starts_with("Fold_enrich") | starts_with("Observed_hits") ,
    names_to = c(".value", "Timepoint"),
    names_sep = " "
  )
#tall_up_data <-  tall_up_data %>% replace_na(list(negative_log10_of_adjusted_p_value = 1, Fold_enrich = 0, Observed_hits=0))


tall_down_data <- all_downreg %>%
  pivot_longer(
    cols = starts_with("negative_log10_of_adjusted_p_value") | starts_with("Fold_enrich") | starts_with("Observed_hits") ,
    names_to = c(".value", "Timepoint"),
    names_sep = " "
  )
# tall_down_data <-  tall_down_data %>% replace_na(list(negative_log10_of_adjusted_p_value = 1, Fold_enrich = 0, Observed_hits=0))
library(dplyr)

# Filter for the "chronic" timepoint and select top 10 hits
top_10_chronic <- tall_down_data %>%
  filter(Timepoint == "Ch_super_down_ergko") %>%  # Replace "chronic" with the actual timepoint name
  arrange(negative_log10_of_adjusted_p_value) %>%               # Sort by p-value (lowest to highest)
  slice_head(n = 20)                 # Select the top 10 rows
top_10_data <- tall_down_data %>%
  filter(term_name %in% top_10_chronic$term_name)

# Create the dot plot
dot_plot <- ggplot(top_10_data, aes(x = Timepoint, y = term_name)) +
  geom_point(aes(size = Fold_enrich, color = -log10(negative_log10_of_adjusted_p_value))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
  scale_size_continuous(range = c(3, 10), name = "Fold Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "GO Term Enrichment Across Timepoints",
    color = "-log10(P-value)",
    size = "Gene Count"
  )
ggsave("/home/kellis/wilson_lab_documents/Figures/ERG_project/Figure_3/downreg_dotplot.pdf",dot_plot )
# Display the plot
print(dot_plot)

# Filter for the "chronic" timepoint and select top 10 hits
top_10_chronic <- tall_up_data %>%
  filter(Timepoint == "72hr_super_up_ergko" & source != "Mouse Phenotype|Human Phenotype") %>%  # Replace "chronic" with the actual timepoint name
  arrange(negative_log10_of_adjusted_p_value) %>%               # Sort by p-value (lowest to highest)
  slice_head(n = 10)                 # Select the top 10 rows
top_10_data <- tall_up_data %>%
  filter(term_name %in% top_10_chronic$term_name)

# Create the dot plot
dot_plot <- ggplot(top_10_data, aes(x = Timepoint, y = term_name)) +
  geom_point(aes(size = Fold_enrich, color = -log10(negative_log10_of_adjusted_p_value))) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P-value)") +
  scale_size_continuous(range = c(3, 10), name = "Fold Enrichment") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(
    title = "GO Term Enrichment Across Timepoints",
    color = "-log10(P-value)",
    size = "Gene Count"
  )


# Display the plot
print(dot_plot)
ggsave("/home/kellis/wilson_lab_documents/Figures/ERG_project/Figure_3/upreg_dotplot.pdf",dot_plot )




### figure 3 heatmaps:
d1_down <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/ERGKO_downregulated_24h_homer/knownResults.txt")
d1_up <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/ERGKO_upregulated_24h_homer/knownResults.txt")

d3_up <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/ERGKO_upregulated_72h_homer/knownResults.txt")
d3_down <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/ERGKO_downregulated_72h_homer/knownResults.txt")

chronic_down <- vroom::vroom("~/R_work/collab/ERG/ATAC/homer_down/knownResults.txt")
chronic_up <- vroom::vroom("~/R_work/collab/ERG/ATAC/homer_up/knownResults.txt")


chronic_down$`Motif Family` <- gsub(".*\\((.*)\\)","\\1",chronic_down$`Motif Name`)
chronic_down$`Motif Name` <- gsub("\\(.*","",chronic_down$`Motif Name`)


chronic_up$`Motif Family` <- gsub(".*\\((.*)\\)","\\1",chronic_up$`Motif Name`)
chronic_up$`Motif Name` <- gsub("\\(.*","",chronic_up$`Motif Name`)


d3_up$`Motif Family` <- gsub(".*\\((.*)\\)","\\1",d3_up$`Motif Name`)
d3_up$`Motif Name` <- gsub("\\(.*","",d3_up$`Motif Name`)


d3_down$`Motif Family` <- gsub(".*\\((.*)\\)","\\1",d3_down$`Motif Name`)
d3_down$`Motif Name` <- gsub("\\(.*","",d3_down$`Motif Name`)


d1_up$`Motif Family` <- gsub(".*\\((.*)\\)","\\1",d1_up$`Motif Name`)
d1_up$`Motif Name` <- gsub("\\(.*","",d1_up$`Motif Name`)


d1_down$`Motif Family` <- gsub(".*\\((.*)\\)","\\1",d1_down$`Motif Name`)
d1_down$`Motif Name` <- gsub("\\(.*","",d1_down$`Motif Name`)


d1_down[,4] <- -log10(d1_down[,3])
d1_down <- d1_down[,c(1,10,4)]
colnames(d1_down)[c(2,3)] <- paste(colnames(d1_down)[c(2,3)], "d1_down", sep = " ")

d1_up[,4] <- -log10(d1_up[,3])
d1_up <- d1_up[,c(1,10,4)]
colnames(d1_up)[c(2,3)] <- paste(colnames(d1_up)[c(2,3)], "d1_up", sep = " ")

d3_down[,4] <- -log10(d3_down[,3])
d3_down <- d3_down[,c(1,10,4)]
colnames(d3_down)[c(2,3)] <- paste(colnames(d3_down)[c(2,3)], "d3_down", sep = " ")

d3_up[,4] <- -log10(d3_up[,3])
d3_up <- d3_up[,c(1,10,4)]
colnames(d3_up)[c(2,3)] <- paste(colnames(d3_up)[c(2,3)], "d3_up", sep = " ")


chronic_down[,4] <- -log10(chronic_down[,3])
chronic_down <- chronic_down[,c(1,10,4)]
colnames(chronic_down)[c(2,3)] <- paste(colnames(chronic_down)[c(2,3)], "chronic_down", sep = " ")


chronic_up[,4] <- -log10(chronic_up[,3])
chronic_up <- chronic_up[,c(1,10,4)]
colnames(chronic_up)[c(2,3)] <- paste(colnames(chronic_up)[c(2,3)], "chronic_up", sep = " ")


library(dplyr)
library(ComplexHeatmap)
library(circlize)

up <- full_join(d1_up,d3_up, by = "Motif Name")
up <- full_join(up,chronic_up, by = "Motif Name")
up <- up[up$`Log P-value d3_up` > 50,]

up.mat <- data.matrix(up[,c(3,5,7)])[c(1:36),]
rownames(up.mat) <- up$`Motif Name`[c(1:36)]


col_up <- colorRamp2(c(0, 300), c("white", "salmon"))
col_up(seq(0,300))

a <- Heatmap(up.mat, col = col_up,show_row_dend = F, show_column_dend = F)



###
down<- full_join(d1_down,d3_down, by = "Motif Name")
down <- full_join(down,chronic_down, by = "Motif Name")
down <- down[down$`Log P-value d3_down` > 50,]

down.mat <- data.matrix(down[,c(3,5,7)])
rownames(down.mat) <- down$`Motif Name`


col_down <- colorRamp2(c(0, 300), c("white", "skyblue"))
col_down(seq(0,300))

b <- Heatmap(down.mat, col = col_down,show_row_dend = F, show_column_dend = F)
ggsave("/home/kellis/wilson_lab_documents/Figures/ERG_project/Figure_3/Heatmaps_upreg.pdf",a )
ggsave("/home/kellis/wilson_lab_documents/Figures/ERG_project/Figure_3/Heatmaps_downreg.pdf",b )



