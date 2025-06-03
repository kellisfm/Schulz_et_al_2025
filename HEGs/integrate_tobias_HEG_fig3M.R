### combine tobias with high evidence genes

library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ComplexHeatmap)
library(dplyr)

setwd("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/BINDetect_output_all_motifs/")
### read in the HEG list
genes <- read.csv("high_evidence_genes.csv")
## split into up and down
HE_up <- genes[genes$RNA24hr_log2fc>0,1]
HE_down <- genes[genes$RNA24hr_log2fc<0,1]
ATAC <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/ATAC/2024_timecourse/Diffbind/RiP_normalized_run/ERG_24_vs_Ctrl_24.csv")

### select only ATAC-peaks assigned to genes consistantly up or down, extract l2fc
atac_heUP <- ATAC[ATAC$SYMBOL%in%HE_up,]
atac_heUP_summed <- atac_heUP |> 
  group_by(across(SYMBOL)) |>
  summarise(ctrl_24_atac_l2fc = sum(Fold)
  )|>
  ungroup() 

colnames(atac_heUP)

atac_heDown <- ATAC[ATAC$SYMBOL%in%HE_down,]
atac_heDown_summed <- atac_heDown |> 
  group_by(across(SYMBOL)) |>
  summarise(ctrl_24_atac_l2fc = sum(Fold)
  )|>
  ungroup() 

### read in and annotate the TOBIAS Factor
dirs <- list.dirs( "~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/BINDetect_output_all_motifs/")
dirs <- dirs[!grepl("beds|plots",dirs)]
dirs <- dirs[-1]
correlation_frame <- data.frame(TF = as.character(),
                                up_correlation = as.numeric(),
                                down_correlation = as.numeric())
for (i in dirs){   #### for each directory in the bindectect directory
  print(i)
  TF <- gsub(".*//","",i)
  tf_path <- paste0(i,"/",TF,"_overview.txt") ## extract the overview file
  target <- vroom::vroom(tf_path) ### read it in 
  
  ### annotate it
  target_red <- target %>% dplyr::select(TFBS_chr, TFBS_start, TFBS_end,Ctrl_24hr_KO_24hr_log2fc, Ctrl_72hr_KO_72hr_log2fc,KO_24hr_bound,KO_72hr_bound)
  target_grange <- makeGRangesFromDataFrame(target_red,
                                         keep.extra.columns=T)
  
  target_anno <- annotatePeak(target_grange, tssRegion=c(-3000, 3000),
                                     TxDb=txd, annoDb= db)
  target_anno <- data.frame(target_anno)
  
  ### select for the rows that are high confidence up or down
  Tobias_heUP <- target_anno[target_anno$SYMBOL%in%HE_up,]
  Tobias_heDown <- target_anno[target_anno$SYMBOL%in%HE_down,]
  
  ### summarize the bound sites
  Tobias_heUP_summed <- Tobias_heUP |> 
    group_by(across(SYMBOL)) |>
    summarise(Ctrl_24hr_KO_24hr_log2fc = sum(Ctrl_24hr_KO_24hr_log2fc)*-1
    )|>
    ungroup()  
  Tobias_heDown_summed <- Tobias_heDown |> 
    group_by(across(SYMBOL)) |>
    summarise(Ctrl_24hr_KO_24hr_log2fc = sum(Ctrl_24hr_KO_24hr_log2fc)*-1
    )|>
    ungroup()  
  
  ### Merge with the ATAC data and correlate
  x <- left_join(atac_heUP_summed,Tobias_heUP_summed, by = "SYMBOL")
  x[is.na(x[,3]),3] <- 0
  up_cor <- cor(x$ctrl_24_atac_l2fc,x$Ctrl_24hr_KO_24hr_log2fc)
  
  y <- left_join(atac_heDown_summed,Tobias_heDown_summed, by = "SYMBOL")
  y[is.na(y[,3]),3] <- 0
  down_cor <- cor(y$ctrl_24_atac_l2fc,y$Ctrl_24hr_KO_24hr_log2fc)
  correlation_frame <- correlation_frame %>% add_row(TF = TF,
                                                     up_correlation = up_cor,
                                                     down_correlation = down_cor)
}
plot(correlation_frame$up_correlation,correlation_frame$down_correlation)

hist(correlation_frame$up_correlation)
hist(correlation_frame$down_correlation)
#write.csv(correlation_frame, "correlation_TOBIAS_atac.csv")

correlation_frame <- read.csv("correlation_TOBIAS_atac.csv")
Gene_to_family <- read.csv("~/TOBIAS/JASPAR_FULL/TF_families.csv")

correlation_frame$TF <- toupper(gsub("_.*","",correlation_frame$TF))
Gene_to_family$TF <- toupper(Gene_to_family$name)

correlation_families <- inner_join(correlation_frame, Gene_to_family, by="TF")
correlation_families[correlation_families$Family %in% c("Fos-related factors", "Jun-related factors"),6] <- "AP-1 factors"
library(ggplot2)
library("dplyr")
library(ggrepel)
library(ggpmisc)

correlation_frame <- correlation_frame[order(abs(correlation_frame$down_correlation),decreasing = T),]
correlation_frame$TF <- with(correlation_frame,reorder(TF,down_correlation,decreasing = T))
ggplot(data = correlation_frame, aes(y=TF, x=down_correlation, colour = Family)) +
  geom_point() +
  theme_classic()

correlation_frame <- correlation_frame[order(abs(correlation_frame$up_correlation),decreasing = T),]
correlation_frame$TF <- with(correlation_frame,reorder(TF,up_correlation,decreasing = T))
ggplot(data = correlation_frame, aes(y=TF, x=up_correlation)) +
  geom_point() +
  theme_classic()

x <- table(correlation_families$Family) > 7
keep <- names(x[x])
correlation_families[!correlation_families$Family%in%keep,6] <- "Other"


correlation_families <- correlation_families[order(correlation_families$down_correlation,decreasing = T),]
correlation_families$TF <- with(correlation_families,reorder(TF,down_correlation,decreasing = T))
correlation_families$Rank <- as.numeric(1:nrow(correlation_families))
ggplot(data = correlation_families, aes(x=Rank, y=down_correlation, colour = Family)) +
  geom_text_repel(data=head(na.omit(correlation_families,40),20), aes(hjust=0.75, vjust=0,label=name),max.overlaps = 15) +
  geom_point() +
  theme_classic()

correlation_families <- correlation_families[order(correlation_families$up_correlation,decreasing = T),]
correlation_families$TF <- with(correlation_families,reorder(TF,up_correlation,decreasing = T))
correlation_families$Rank <- as.numeric(1:nrow(correlation_families))
ggplot(data = correlation_families, aes(x=Rank, y=up_correlation, colour = Family)) +
  geom_text_repel(data=head(na.omit(correlation_families,40),20), aes(hjust=0.75, vjust=0,label=name),max.overlaps = 15) +
  geom_point() +
  theme_classic()

rank_correlation_plot <- ggplot(merged_RNA_tobias, aes(x=ERGWT.Rank, y=ERGKO.Rank)) +
  geom_point(data=merged_RNA_tobias, aes(colour=-log2(ERGWT_ERGKO_pvalue),shape=WTvsKOsig,size=-log2(ERGWT_ERGKO_pvalue_tobias))) + 
  theme_minimal() +
  #ylim(9,12.25) +
  geom_text_repel(data=head(na.omit(merged_RNA_tobias[(merged_RNA_tobias[,15] == T),]),40), aes( hjust=0.75, vjust=0, label=name))+
  stat_poly_line() +
  stat_poly_eq() 
