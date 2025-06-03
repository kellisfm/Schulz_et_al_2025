### generating and plotting HEGS
# Kai Ellis 2025
setwd("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/")
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
db <- "org.Hs.eg.db" 
txd <- TxDb.Hsapiens.UCSC.hg38.knownGene  

annotate_tf_txt <- function(TF_binding_path){
  db <- "org.Hs.eg.db" 
  txd <- TxDb.Hsapiens.UCSC.hg38.knownGene  

  jun_data <- vroom::vroom(TF_binding_path)
  jun_data_red <- jun_data %>% dplyr::select(TFBS_chr, TFBS_start, TFBS_end,Ctrl_24hr_KO_24hr_log2fc, Ctrl_72hr_KO_72hr_log2fc,KO_24hr_bound,KO_72hr_bound,Ctrl_24hr_bound,Ctrl_72hr_bound)
  jun_grange <- makeGRangesFromDataFrame(jun_data_red,
                                         keep.extra.columns=T)
  
  dbReport_empa_anno <- annotatePeak(jun_grange, tssRegion=c(-3000, 3000),
                                     TxDb=txd, annoDb= db)
  dbReport_empa_anno <- data.frame(dbReport_empa_anno)


  dbReport_empa_summed <- dbReport_empa_anno |> 
    group_by(across(SYMBOL)) |>
    summarise(Ctrl_72hr_KO_72hr_log2fc = sum(Ctrl_72hr_KO_72hr_log2fc),
              Ctrl_24hr_KO_24hr_log2fc = sum(Ctrl_24hr_KO_24hr_log2fc),
              KO_bound = sum(KO_24hr_bound|KO_72hr_bound == T)>1,
              Ctrl_bound = sum(Ctrl_24hr_bound|Ctrl_72hr_bound == T)>1,
              d1_bound =  sum(KO_24hr_bound == T)>1,
              d3_bound =  sum(KO_72hr_bound == T)>1
    )|>
    ungroup()  
  #df1 = dbReport_empa_promoters[order(dbReport_empa_promoters[,'SYMBOL'],-dbReport_empa_promoters[,'Ctrl_72hr_KO_72hr_log2fc']),]
  #df = df1[!duplicated(df1$SYMBOL),]
  return(dbReport_empa_summed)
}



### processing randi data
Randi_h3K27_down <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ChIP/Randi_2019/siRNA_ERG/Peaks/diffrential_peaks_ERGKO_DOWN_MACS.txt")
Randi_h3K27_down$`Fold Change vs. Background` <- log10(Randi_h3K27_down$`Fold Change vs. Background`+1)
Randi_h3K27_up <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ChIP/Randi_2019/siRNA_ERG/Peaks/diffrential_peaks_ERGKO_UP_MACS.txt")
Randi_h3K27_up$`Fold Change vs. Background` <- log10(Randi_h3K27_up$`Fold Change vs. Background`)*-1
Randi_h3K27_down <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ChIP/Randi_2019/siRNA_ERG/Peaks/diffrential_peaks_ERGKO_DOWN_MACS.txt")
Randi_h3K27_down$`Fold Change vs. Background` <- log10(Randi_h3K27_down$`Fold Change vs. Background`+1)
Randi_h3K27_nc <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ChIP/Randi_2019/siRNA_ERG/Peaks/same_peaks_ERGKO_MACS.txt")
Randi_h3K27_nc <- data.frame(Randi_h3K27_nc)

up_in_erg <- Randi_h3K27_nc$Total.Tags < Randi_h3K27_nc$Background.Tags
Randi_h3K27_nc$`Fold.Change.vs..Background` <- log10(Randi_h3K27_nc$`Fold.Change.vs..Background`+1)
Randi_h3K27_nc$Fold.Change.vs..Background[up_in_erg] <- Randi_h3K27_nc$Fold.Change.vs..Background[up_in_erg]*-1
colnames(Randi_h3K27_nc) <- colnames(Randi_h3K27_down)




randi_h3k27 <- rbind(Randi_h3K27_up,Randi_h3K27_down)
randi_h3k27 <- rbind(randi_h3k27,Randi_h3K27_nc)
#write.csv(randi_h3k27,"~/mnt/largeprojects/Kai/collab/ERG/ChIP/Randi_2019/siRNA_ERG/Peaks/all_randi_peaks.csv")
ggplot(data = randi_h3k27, aes(x = `Fold Change vs. Background`, y= -log10(`p-value`))) +
  geom_point() +
  theme_classic()

randi_grange <- makeGRangesFromDataFrame(randi_h3k27,
                                       keep.extra.columns=T)




randi_anno <- annotatePeak(randi_grange, tssRegion=c(-3000, 3000),
                                   TxDb=txd, annoDb= db)
randi_anno <- data.frame(randi_anno)
randi_anno_summed <- randi_anno |> 
  group_by(across(SYMBOL)) |>
  summarise(l2fc_ChIP = sum(Fold.Change.vs..Background),
            min_padj_ChIP = min(p.value)) |>
  ungroup()

# tmp <- vroom::vroom("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/inputs/Peaks/ERG_merged_peaks.bed")
# colnames(tmp) <- c("chr","start","end")
# tmp <- makeGRangesFromDataFrame(tmp,keep.extra.columns=T)
# 
# tmp <- annotatePeak(tmp, tssRegion=c(-3000, 3000),
#                            TxDb=txd, annoDb= db)
# write.table(tmp, "~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/inputs/Peaks/merged_anno_peaks.bed", sep = "\t", row.names = F, col.names = F, quote = F)

randi_anno <- data.frame(randi_anno_summed)

#### read in RNA24 hr
rna_de_24 <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/RNA/ERG_siRNA/01RawData/res_Ctrl24_Erg24.csv")
colnames(rna_de_24)[1] <- "SYMBOL"
colnames(rna_de_24)[3] <- "RNA24hr_log2fc"
colnames(rna_de_24)[7] <- "RNA24hr_padj"



rna_de_24 <- rna_de_24[,c(1,3,7)]
rna_de_24$RNA24hr_log2fc <- rna_de_24$RNA24hr_log2fc*-1
rna_de_24$diffexpressed <- "Non-DE"
rna_de_24$diffexpressed[rna_de_24$RNA24hr_log2fc>0.1 & rna_de_24$RNA24hr_padj<0.05 ] <- "UP"
rna_de_24$diffexpressed[rna_de_24$RNA24hr_log2fc<0.1 & rna_de_24$RNA24hr_padj<0.05 ] <- "DOWN"

ggplot(data = rna_de_24, aes(x = RNA24hr_log2fc*-1, y= -log10(RNA24hr_padj), colour = diffexpressed)) +
  geom_point() +
  xlim(c(-10,10)) + 
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  ggrepel::geom_label_repel(data = utils::head(rna_de_24, 30), aes(label = SYMBOL), size = 3, label.size = NA, 
                            max.overlaps = 30, box.padding = 0.4) +
  theme_classic()



#### read in RNA72 hr
rna_de_72 <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/RNA/ERG_siRNA/01RawData/res_Ctrl72_Erg72.csv")
colnames(rna_de_72)[1] <- "SYMBOL"
colnames(rna_de_72)[3] <- "RNA72hr_log2fc"
colnames(rna_de_72)[7] <- "RNA72hr_padj"
rna_de_72 <- rna_de_72[,c(1,3,7)]
rna_de_72$RNA72hr_log2fc <- rna_de_72$RNA72hr_log2fc*-1
rna_de_72$diffexpressed <- "Non-DE"
rna_de_72$diffexpressed[rna_de_72$RNA72hr_log2fc>0.1 & rna_de_72$RNA72hr_padj<0.05 ] <- "UP"
rna_de_72$diffexpressed[rna_de_72$RNA72hr_log2fc<0.1 & rna_de_72$RNA72hr_padj<0.05 ] <- "DOWN"

ggplot(data = rna_de_72, aes(x = RNA72hr_log2fc*-1, y= -log10(RNA72hr_padj), colour =  diffexpressed)) +
  geom_point() +
  xlim(c(-10,10)) + 
 # geom_text_repel() +
 scale_color_manual(values=c("salmon", "grey", "skyblue")) +
  geom_vline(xintercept=c(-1, 1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  ggrepel::geom_label_repel(data = utils::head(rna_de_72, 30), aes(label = SYMBOL), size = 3, label.size = NA, 
                            max.overlaps = 30, box.padding = 0.4) +
  theme_classic()


#### read in ATAC
ATAC <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/ERGKO_v_Ctrl_24h.csv")

ATAC <- ATAC[,c(8:10,1:7)] 

ATAC24_grange <- makeGRangesFromDataFrame(ATAC,
                                         keep.extra.columns=T)

ATAC24_anno <- annotatePeak(ATAC24_grange, tssRegion=c(-3000, 3000),
                           TxDb=txd, annoDb= db)
ATAC24_anno <- data.frame(ATAC24_anno)
ATAC_24_summed <- ATAC24_anno |> 
  group_by(across(SYMBOL)) |>
  summarise(l2fc_ATAC_24h = sum(log2FoldChange),
            min_padj_ATAC_24h = min(padj)) |>
  ungroup()

merged_24 <- dplyr::inner_join(rna_de_24,ATAC_24_summed, by= "SYMBOL")
merged_24$directionally_uniform <- rowSums(merged_24[,c(2,4)]>0) == ncol(merged_24[,c(2,4)]) | rowSums(merged_24[,c(2,4)]<0) == ncol(merged_24[,c(2,4)]) 
merged_24$all_significant <-  rowSums(merged_24[,c(3,5)]<0.05) == ncol(merged_24[,c(3,5)])
high_evidence_genes_24 <- merged_24[merged_24$directionally_uniform & merged_24$all_significant,]
write.csv(high_evidence_genes_24, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_24hr.csv")
high_evidence_genes_24 <- read.csv("/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_24hr.csv")

high_evidence_genes_24_nearby_peaks <- dplyr::inner_join(high_evidence_genes_24,ATAC24_anno, by = "SYMBOL")
write.csv(high_evidence_genes_24_nearby_peaks, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_nearby_atac_peaks_24hr.csv")
high_evidence_genes_24_nearby_peaks_significant <- high_evidence_genes_24_nearby_peaks[high_evidence_genes_24_nearby_peaks$padj<0.05,]
high_evidence_genes_24_nearby_peaks_significant_uniform <- high_evidence_genes_24_nearby_peaks_significant[rowSums(high_evidence_genes_24_nearby_peaks_significant[,c(2,15)]>0) == ncol(high_evidence_genes_24_nearby_peaks_significant[,c(2,15)]) | rowSums(high_evidence_genes_24_nearby_peaks_significant[,c(2,15)]<0) == ncol(high_evidence_genes_24_nearby_peaks_significant[,c(2,15)]) ,]
write.csv(high_evidence_genes_24_nearby_peaks_significant_uniform, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_nearby_sig_uniform_atac_peaks_24hr.csv")

ATAC_72 <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/ERGKO_v_Ctrl_72h.csv")
#readr::write_delim(ATAC[,-1], "~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/Diffbind/RiP_normalized_run/ERG_24_vs_Ctrl_24.bed", col_names = F, quote = "none", delim = "\t")
ATAC_72 <- ATAC_72[,c(8:10,1:7)] 

ATAC72_grange <- makeGRangesFromDataFrame(ATAC_72,
                                          keep.extra.columns=T)

ATAC72_anno <- annotatePeak(ATAC72_grange, tssRegion=c(-3000, 3000),
                            TxDb=txd, annoDb= db)
ATAC72_anno <- data.frame(ATAC72_anno)
ATAC_72_summed <- ATAC72_anno |> 
  group_by(across(SYMBOL)) |>
  summarise(l2fc_ATAC_24h = sum(log2FoldChange),
            min_padj_ATAC_24h = min(padj)) |>
  ungroup()



merged_72 <- dplyr::inner_join(rna_de_72,ATAC_72_summed, by= "SYMBOL")
merged_72$directionally_uniform <- rowSums(merged_72[,c(2,4)]>0) == ncol(merged_72[,c(2,4)]) | rowSums(merged_72[,c(2,4)]<0) == ncol(merged_72[,c(2,4)]) 
merged_72$all_significant <-  rowSums(merged_72[,c(3,5)]<0.05) == ncol(merged_72[,c(3,5)])
high_evidence_genes_72 <- merged_72[merged_72$directionally_uniform & merged_72$all_significant,]
write.csv(high_evidence_genes_72, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_72hr.csv")


high_evidence_genes_72_nearby_peaks <- dplyr::inner_join(high_evidence_genes_72,ATAC72_anno, by = "SYMBOL")
high_evidence_genes_72_nearby_peaks_significant <- high_evidence_genes_72_nearby_peaks[high_evidence_genes_72_nearby_peaks$padj<0.05,]
high_evidence_genes_72_nearby_peaks_significant_uniform <- high_evidence_genes_72_nearby_peaks_significant[rowSums(high_evidence_genes_72_nearby_peaks_significant[,c(2,15)]>0) == ncol(high_evidence_genes_72_nearby_peaks_significant[,c(2,15)]) | rowSums(high_evidence_genes_72_nearby_peaks_significant[,c(2,15)]<0) == ncol(high_evidence_genes_72_nearby_peaks_significant[,c(2,15)]) ,]


write.csv(high_evidence_genes_72_nearby_peaks, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_nearby_atac_peaks_72hr.csv")
write.csv(high_evidence_genes_72_nearby_peaks_significant_uniform, "/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/pseudoDiffbind/integrated_RNA/HEG_nearby_sig_uniform_atac_peaks_72hr.csv")

####
mat <- matrix(merged_72$RNA72hr_log2fc*-1)
mat2 <- matrix(merged_72$l2fc_ATAC_72*-1)
mat3 <- matrix(as.numeric(merged_72$directionally_uniform))
mat4 <- matrix(as.numeric(merged_72$all_significant))

a <- Heatmap(mat) #add col=col_fun to alter color palette
b <- Heatmap(mat2)
c <- Heatmap(mat3)
d <- Heatmap(mat4)
a + b + c + d 


####Erg ChIP
ERG <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/Diffbind/RiP_normalized_run/ERG_bound_changing_peaks.bed", col_names = F)
ERG <- data.frame(ERG)
ERG <- ERG[,c(32,33)]
colnames(ERG) <- c("SYMBOL","Full_name")
ERG_summed <- ERG |> 
  group_by(across(SYMBOL)) |>
  summarise(n_ERG_bound = nrow(Full_name),
  ) |>
  ungroup()
ERG <- data.frame(table(ERG$SYMBOL))
colnames(ERG) <- c("SYMBOL","N_erg_bound")
####


Fosl1 <- annotate_tf_txt("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/BINDetect_output_all_motifs/FOSL1_MA0477.2/FOSL1_MA0477.2_overview.txt")
colnames(Fosl1)[c(2:6)] <- paste0("FOSL1_",colnames(Fosl1)[c(2:6)]) 
Jun <- annotate_tf_txt("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/BINDetect_output_all_motifs/Jun_MA0489.2/Jun_MA0489.2_overview.txt")
colnames(Jun)[c(2:6)] <- paste0("Jun_",colnames(Jun)[c(2:6)]) 
Erg <- annotate_tf_txt("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/BINDetect_output_all_motifs/Erg_MA0474.3/Erg_MA0474.3_overview.txt")
colnames(Erg)[c(2:6)] <- paste0("ERG_",colnames(Erg)[c(2:6)]) 
Erg

merged <- dplyr::inner_join(data.frame(Fosl1),rna_de_24, by= "SYMBOL")
#merged <- dplyr::right_join(data.frame(Erg),merged, by= "SYMBOL")  
merged <- dplyr::inner_join(merged,rna_de_72, by= "SYMBOL")  
merged <- dplyr::inner_join(merged,randi_anno, by= "SYMBOL")  
merged <- dplyr::inner_join(merged,ATAC_summed, by= "SYMBOL")  


merged <- data.frame(merged)
rownames(merged) <- merged$SYMBOL
merged <- merged[,-1]
merged$sum_l2fc <- rowSums(merged)
hist(merged$sum_l2fc)

merged<- merged*-1
merged$FOSL1_KO_bound <- merged$FOSL1_KO_bound*-1
plot(merged$sum_l2fc,merged$RNA72hr_log2fc)
plot(merged$l2fc_ATAC,merged$l2fc_ChIP)
merged$directionally_uniform <- rowSums(merged>0) == ncol(merged)  | rowSums(merged<0) == ncol(merged) 


merged$directionally_uniform[is.na(merged$directionally_uniform)] <- 0
merged$sum_l2fc_modified <- merged$sum_l2fc
merged$sum_l2fc_modified[!merged$directionally_uniform] <- merged$sum_l2fc_modified[!merged$directionally_uniform]/2
#write.csv(merged,"ERG_FOSL1_ATAC_CHIP_RNA_MERGE.csv")




library(circlize)
col_fun = colorRamp2(c(-10, 0, 10), c("blue", "white", "red"))
col_fun(seq(-3, 3))


mat <- matrix(merged$FOSL1_Ctrl_24hr_KO_24hr_log2fc)
mat2 <- matrix(merged$RNA24hr_log2fc)
mat3 <- matrix(merged$l2fc_ChIP)
mat4 <- matrix(merged$l2fc_ATAC)
mat5 <- matrix(as.numeric(merged$directionally_uniform))
mat6 <- matrix(merged$sum_l2fc)
mat7 <- matrix(merged$FOSL1_KO_bound)
library(ComplexHeatmap)
a <- Heatmap(mat, col = col_fun) #add col=col_fun to alter color palette
b <- Heatmap(mat2)
c <- Heatmap(mat3)
d <- Heatmap(mat4)
e <- Heatmap(mat5)
f <- Heatmap(mat6)
g <- Heatmap(mat7)
a + b +  d + f + e + g
plot(mat3,mat4)

sum(merged$Ctrl_24hr_KO_24hr_log2fc > 0)


#### setup a csv for export:


top_0.05 <- head(erg_merged,90)
bot_0.05 <- tail(erg_merged,90)
write.csv(top_0.05,"ap1_top_evidence.csv")
write.csv(bot_0.05,"ap1_top_evidence_down.csv")



write.csv(merged,"combined_evidence.csv")



#### no jun:


merged_all <- dplyr::inner_join(rna_de_24,rna_de_72, by= "SYMBOL")  
merged_all <- dplyr::inner_join(merged_all,randi_anno, by= "SYMBOL")  
merged_all <- dplyr::inner_join(merged_all,ATAC_24_summed, by= "SYMBOL") 
merged_all <- dplyr::left_join(merged_all,Erg, by = "SYMBOL",)
merged_all[is.na(merged_all$ERG_KO_bound),c(10:13)] <- 0
merged_all[is.na(merged_all$ERG_Ctrl_bound),c(10:13)] <- 0
merged_all <- data.frame(merged_all)
rownames(merged_all) <- merged_all$SYMBOL
merged_all <- merged_all[,-1]
merged_all[,c(1,3,5,7)]<- merged_all[,c(1,3,5,7)]*-1
merged_all$Ctrl_24hr_KO_24hr_log2fc <- merged_all$ERG_Ctrl_24hr_KO_24hr_log2fc*-1
merged_all$Ctrl_72hr_KO_72hr_log2fc <- merged_all$ERG_Ctrl_72hr_KO_72hr_log2fc*-1 
merged_all$sum_l2fc <- rowSums(merged_all[,c(1,3,5,7)])
merged_all$directionally_uniform <- rowSums(merged_all[,c(1,3,5,7)]>0) == ncol(merged_all[,c(1,3,5,7)]) | rowSums(merged_all[,c(1,3,5,7)]<0) == ncol(merged_all[,c(1,3,5,7)]) 
merged_all$sum_l2fc_modified <- merged_all$sum_l2fc
merged_all$sum_l2fc_modified[!merged_all$directionally_uniform] <- merged_all$sum_l2fc_modified[!merged_all$directionally_uniform]/2
merged_all$all_significant <-  rowSums(merged_all[,c(2,4,6,8)]<0.05) == ncol(merged_all[,c(2,4,6,8)])
#write.csv(merged_all,"combined_RNA_ChIP_ATAC.csv")

high_evidence_genes <- merged_all[merged_all$directionally_uniform & merged_all$all_significant,]
high_evidence_genes$gain_of_Jun_binding <- high_evidence_genes$Jun_KO_bound>high_evidence_genes$Jun_Ctrl_bound
high_evidence_genes$loss_of_Erg_binding <- high_evidence_genes$ERG_KO_bound<high_evidence_genes$ERG_Ctrl_bound
sum(high_evidence_genes[high_evidence_genes$RNA24hr_log2fc<0,12])
sum(high_evidence_genes[high_evidence_genes$RNA24hr_log2fc<0,19])
#write.csv(high_evidence_genes, "high_evidence_genes_Jun.csv")
mat2 <- matrix(high_evidence_genes$RNA24hr_log2fc)
mat3 <- matrix(high_evidence_genes$l2fc_ChIP)
mat4 <- matrix(high_evidence_genes$l2fc_ATAC)
mat5 <- matrix(high_evidence_genes$sum_l2fc)
mat6 <- matrix(as.numeric(merged_all$directionally_uniform))
mat7 <- matrix(as.numeric(merged_all$all_significant))
mat8 <- matrix(high_evidence_genes$Ctrl_24hr_KO_24hr_log2fc)
mat9 <- matrix(as.numeric(high_evidence_genes$Jun_Ctrl_bound))
mat10 <- matrix(as.numeric(high_evidence_genes$Jun_KO_bound))
mat11 <- matrix(as.numeric(high_evidence_genes$Jun_KO_bound>high_evidence_genes$Jun_Ctrl_bound))
b <- Heatmap(mat2)
c <- Heatmap(mat3)
d <- Heatmap(mat4)
e <- Heatmap(mat5)
f <- Heatmap(mat6)
g <- Heatmap(mat7)
h <- Heatmap(mat8)
i <- Heatmap(mat9)
j <- Heatmap(mat10)
k <- Heatmap(mat11)
b + c + d + h + i + j + k

col_fun = colorRamp2(c(0, 5, 10), c("blue", "white", "red"))

mat2 <- matrix(high_evidence_genes$RNA24hr_log2fc)
mat3 <- matrix(high_evidence_genes$l2fc_ChIP)
mat4 <- matrix(high_evidence_genes$l2fc_ATAC)
mat5 <- matrix(high_evidence_genes$Ctrl_24hr_KO_24hr_log2fc)
mat6 <- matrix(high_evidence_genes$Ctrl_72hr_KO_72hr_log2fc)
mat7 <- matrix(as.numeric(high_evidence_genes$KO_bound))

b + c + d #+ e + f + g

top_0.05 <- head(erg_merged_all,90)
bot_0.05 <- tail(erg_merged_all,90)
write.csv(top_0.05,"ap1_top_evidence.csv")
write.csv(bot_0.05,"ap1_top_evidence_down.csv")

hist(merged_all$sum_l2fc)

plot(merged_all$sum_l2fc,merged_all$RNA72hr_log2fc)


library(ggplot2)
library(ggpubr)

setwd("/home/kellis/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/Figures")
pdf("ATAC_CHIP_CORRELATION.pdf")
ggscatter(merged_all, x = "l2fc_ATAC_24h", y = "l2fc_ChIP", add = "reg.line") +
  stat_cor(label.x = -10, label.y = 20) +
  stat_regline_equation(label.x = -10, label.y = 15) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
dev.off()

pdf("ATAC_RNA_CORRELATION.pdf")
ggscatter(merged_all, x = "l2fc_ATAC_24h", y = "RNA24hr_log2fc", add = "reg.line") +
  stat_cor(label.x = -10, label.y = 5) +
  stat_regline_equation(label.x = -10, label.y = 4) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
dev.off()
merged_all$l2fc_ChIP <-merged_all$l2fc_ChIP*-1

setwd("/home/kellis/mnt/largeprojects/Kai/collab/ERG/Integrated/24hour/Figures/")
pdf("ChIP_RNA_CORRELATION.pdf")
ggscatter(merged_all, x = "l2fc_ChIP", y = "RNA24hr_log2fc", add = "reg.line") +
  stat_cor(label.x = -15, label.y = 3.5) +
  stat_regline_equation(label.x = -15, label.y = 4) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
dev.off()

endmt_genes <- read.csv("~/Downloads/composite.gene.lists.all(1).csv")

library(dplyr)
Jun_endMt <- Jun[Jun$SYMBOL%in%endmt_genes$EndMT.genes[1:242],]

Jun_endMt$Ctrl_72hr_KO_72hr_log2fc <- Jun_endMt$Ctrl_72hr_KO_72hr_log2fc*-1
Jun_endMt$Ctrl_24hr_KO_24hr_log2fc <- Jun_endMt$Ctrl_24hr_KO_24hr_log2fc*-1

Jun_endMt$KO_bound <- as.numeric(Jun_endMt$KO_bound)*10
Jun_endMt$Ctrl_bound <- as.numeric(Jun_endMt$Ctrl_bound)*10
Jun_endMt$d1_bound <- as.numeric(Jun_endMt$d1_bound)*10
Jun_endMt$d3_bound <- as.numeric(Jun_endMt$d3_bound)*10

rna_de_24$`RNA24hr_-log10(padj)` <- -log10(rna_de_24$RNA24hr_padj)
Jun_endMt <- inner_join(Jun_endMt, rna_de_24[,-3], by = "SYMBOL")

mat3 <- matrix(as.numeric(Jun_endMt$Ctrl_bound))
mat4 <- matrix(as.numeric(Jun_endMt$d1_bound))
mat5 <- matrix(as.numeric(Jun_endMt$d3_bound))
mat2 <- matrix(Jun_endMt$Ctrl_24hr_KO_24hr_log2fc)
mat <- matrix(Jun_endMt$Ctrl_72hr_KO_72hr_log2fc)



mat <- as.matrix(data.frame(Jun_endMt[,-1]))
rownames(mat) <- Jun_endMt$SYMBOL

library(circlize)
col_fun = colorRamp2(c(, 0,), c("skyblue", "white", "salmon"))
col_fun(seq(-3, 3))

library(ComplexHeatmap)
Heatmap(mat,  row_names_gp = grid::gpar(fontsize = 6)) #add col=col_fun to alter color palette
b <- Heatmap(mat2)
c <- Heatmap(mat3)
d <- Heatmap(mat4)
e <- Heatmap(mat5)

a + b + c + d + e

scaled_mat = t(scale(t(mat)))



write.csv(Jun_endMt, "Jun_endMT.csv")

