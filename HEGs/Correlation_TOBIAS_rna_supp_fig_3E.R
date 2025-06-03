#### running a correlation between TOBIAS and RNA-seq expression.
Gene_to_family <- read.csv("~/TOBIAS/JASPAR_FULL/TF_families.csv")
TOBIAS_data <- data.frame(vroom::vroom("~/TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/BINDetect_output_all_motifs/bindetect_results.txt"))
TOBIAS_data_red <- TOBIAS_data[,c(2,10,12,14,16)]
TOBIAS_data_red$name <- toupper(TOBIAS_data_red$name)
es <- floor(read.csv("~/mnt/largeprojects/Kai/collab/ERG/RNA/ERG_siRNA/01RawData/Steve/processed_data/salmon_counts_matrix.csv", row.names = 1))

x <- read.delim("~/mnt/largeprojects/Kai/collab/ERG/RNA/ERG_siRNA/01RawData/Steve/siRNA_KD_metadata.txt", sep="\t")
set.seed(1234)
x$Group
colnames(es) <- x$Group
es <- es[,c(1,5,9,13,3,7,11,15,2,6,10,14,4,8,12,16)]
library(limma)
library(dplyr)
es <- normalizeBetweenArrays(log2(es), method="quantile")

es <- es[order(rowMeans(es), decreasing=TRUE), ]

#es <- es[head(order(rowMeans(es), decreasing=TRUE), 12000), ]
head(es)


es <- es[rowSums(is.na(es)) == 0,]
es <- es[rowSums(is.infinite(es)) == 0,]
es <- data.frame(es)
es_mod <- mutate(es, Ctrl.24h_mean = rowMeans(select(es,c(1:4)), na.rm = TRUE), 
                 ctrl.72h_mean = rowMeans(select(es,c(5:9)) ),
                 Erg.24h_mean = rowMeans(select(es,c(10:13)) ),
                 Erg.72h_mean = rowMeans(select(es,c(13:16)) )
                 )
es_red <- es_mod[,c(17:20)]                 
es_red$name <- row.names(es_red)
es_redder <- es_red[(es_red$name%in%TOBIAS_data_red$name),]
rownames(es_redder) <- es_redder$name
es_redder <- es_redder[ order(row.names(es_redder)),c(1:4) ]

colnames(TOBIAS_data_red)
TOBIAS_data_red <- TOBIAS_data_red |> 
  group_by(across(name)) |>
  summarise(Ctrl_24hr_mean_score = sum(Ctrl_24hr_mean_score),
            Ctrl_72hr_mean_score = sum(Ctrl_72hr_mean_score),
            KO_24hr_mean_score = sum(KO_24hr_mean_score),
            KO_72hr_mean_score = sum(KO_72hr_mean_score)) |>
  ungroup()
TOBIAS_data_red <- data.frame(TOBIAS_data_red)
rownames(TOBIAS_data_red) <- TOBIAS_data_red$name
TOBIAS_data_redder <- TOBIAS_data_red[(TOBIAS_data_red$name%in%es_red$name),]

A <- as.matrix(TOBIAS_data_redder[,c(2:5)])
B <- as.matrix(es_redder)
es_redder$cor_values <- sapply(seq.int(dim(A)[1]), function(i) cor(A[i,], B[i,]))
tfs_for_plotting <- inner_join(Gene_to_family, TOBIAS_data_red, by = "name")


tfs_for_plotting <- tfs_for_plotting[ order(tfs_for_plotting$name),]

es_redder$name <- row.names(es_redder)
tfs_for_plotting <- inner_join(es_redder,tfs_for_plotting, by = "name")
tfs_for_plotting[tfs_for_plotting$Family %in% c("Fos-related factors", "Jun-related factors"),7] <- "AP-1 factors"



x <- table(tfs_for_plotting$Family) > 7
keep <- names(x[x])
tfs_for_plotting <- tfs_for_plotting[tfs_for_plotting$Family%in%keep,]
tfs_for_plotting <- tfs_for_plotting[ order(tfs_for_plotting$Family),]

#### read in RNA24 hr
rna_de_24 <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/RNA/ERG_siRNA/01RawData/res_Ctrl24_Erg24.csv")
colnames(rna_de_24)[1] <- "name"
colnames(rna_de_24)[7] <- "RNA24hr_-log10(padj)"
rna_de_24[,7] <- log10(rna_de_24[,7])*-1  
rna_de_24 <- rna_de_24[,c(1,7)]

#### read in RNA72 hr
rna_de_72 <- vroom::vroom("~/mnt/largeprojects/Kai/collab/ERG/RNA/ERG_siRNA/01RawData/res_Ctrl72_Erg72.csv")
colnames(rna_de_72)[1] <- "name"
colnames(rna_de_72)[7] <- "RNA72hr_-log10(padj)"
rna_de_72[,7] <- log10(rna_de_72[,7])*-1  
rna_de_72 <- rna_de_72[,c(1,7)]

rna <- inner_join(rna_de_24, rna_de_72, by = "name")
rna$average_log10padj <- rowMeans(rna[,c(2,3)])
tfs_for_plotting_de <- data.frame(inner_join(rna, tfs_for_plotting, by = "name"))

library(ggplot2)
setwd(dir = "TOBIAS/TOBIAS_tmp_ERG_run/timecourse_2024/")
pdf("correlation_TOBIAS_RNA-seq_with_legend.pdf", width = 10)
ggplot(data = tfs_for_plotting_de, aes(x = Family, y = cor_values, size = abs(average_log10padj) )) +
  geom_point() +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
   ggrepel::geom_label_repel(data = subset(tfs_for_plotting, abs(cor_values) > 
                                           0.75), aes(label = name), size = 3, label.size = NA, 
                             max.overlaps = 30, box.padding = 0.4)

dev.off()

unique(tfs_for_plotting$Family)






#####
p <- ggplot2::ggplot(res, aes(x = new, y = Correlation)) + 
  ggplot2::geom_hline(yintercept = -0.5, linetype = "dashed", 
                      color = "grey") + ggplot2::geom_hline(yintercept = 0.5, 
                                                            linetype = "dashed", color = "grey") + 
  ggplot2::geom_jitter(aes(color = as.factor(new)), 
                                                                                                                        size = abs(res$Correlation) * 4, width = 0.4) + 
  ggplot2::scale_color_manual(values = color) + ggpubr::theme_pubr() + 
  ggplot2::labs(x = "", y = "Correlation of TF footprint score and TF expression") + 
  ggrepel::geom_label_repel(data = subset(res, is_annotate == 
                                            "yes"), aes(label = label), size = 3, label.size = NA, 
                            max.overlaps = 30, box.padding = 0.4) + ggplot2::geom_tile(data = dfcol, 
                                                                                       aes(x = x, y = y), height = 0.2, color = "black", 
                                                                                       fill = paletteer::paletteer_d("ggthemes::Classic_Cyclic"), 
                                                                                       alpha = 0.8, show.legend = F) + ggplot2::geom_text(data = dfcol, 
                                                                                                                                          aes(x = x, y = y, label = x), size = 4, color = "white") + 
  ggplot2::theme(panel.border = element_rect(colour = "black", 
                                             fill = NA)) + ggplot2::theme(legend.position = "none")
matrix()
