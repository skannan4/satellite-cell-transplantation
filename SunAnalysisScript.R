library(scran)
library(ggplot2)
library(ggrepel)
library(monocle)
library(pheatmap)
library(scater)
library(singleCellNet)
library(Matrix)

#Note to users - change this working directory to wherever you downloaded the relevant files
setwd("~/Documents/Research/SunnySun/FinalManuscript")
#Load in data_____________________________________________________________________________________________________________________
#Note that these files have already been converted to gene symbols using biomaRt

data_no_ercc = read.csv("data_no_ercc.csv", as.is = TRUE, row.names = 1)
data_w_ercc = read.csv("data_w_ercc.csv", as.is = TRUE, row.names = 1)
phenotype = read.csv("phenotype.csv", as.is = TRUE, row.names = 1)

#Pre-normalization quality control________________________________________________________________________________________________
#We used the following criteria:
#1) Remove cells with low reads (<200K reads per cell)
#2) Remove potential doublets (>100K UMIs)
#3) Remove potentially lysed cells (>40% ERCC reads)
#4) Remove any remaining low quality cells (<2000 genes)
#These filtering criteria were selected based on plotting parameters such as reads, UMIs, and ERCC reads and looking for reasonable cutoffs; please feel free to set different thresholds as appropriate

good_cells = rownames(phenotype)[phenotype$reads > 200000 & phenotype$umicounts < 100000 & phenotype$reads_ercc_percentage < 40 & phenotype$genes > 2000 & phenotype$genes < 12000 ]

#Scater normalization by ERCCs_____________________________________________________________________________________________________
#Our rationale here was that, unlike in some biological questions/experimental contexts, RNA content differences here could be relevant to our phenomenon of interest. We spiked in external control RNAs (ERCCs), and used scran to normalize our samples to these ERCCs, rather than to the library size (which assumes uniform RNA content). This result will produce different sets of differentially expressed genes! The interpretation of that must be relevant to the biological question; here, we felt that absolute RNA content differences were likely to be relevant.
#Please note that, as part of a calibration, the 1 month in vitro samples were spiked with exactly twice as much ERCCs as other samples; hence, for accurate quantitation, the size factors for these samples are divided by two. It is easy to validate that this produces acceptable results.

sce = SingleCellExperiment(list(counts=as.matrix(data_w_ercc[, good_cells])), colData = phenotype[good_cells, ])
isSpike(sce, "MySpike") = grep("^ERCC", rownames(sce))
sce <- computeSpikeFactors(sce)
spike_factors = sce@int_colData$size_factor_MySpike
spike_factors[phenotype[good_cells, ]$celltype == "in_vitro_1_month" ] = spike_factors[phenotype[good_cells, ]$celltype == "in_vitro_1_month" ]/2
sce@int_colData$size_factor_MySpike = spike_factors
sce@int_colData$size_factor = spike_factors

#Clustering and differential gene expression analysis_______________________________________________________________________________
#We now perform differential gene expression testing in Monocle.
#EDIT 07/11: We used Monocle 2 for this analysis. While Monocle 3 is now available, for this sample size we found Monocle 2 to be sufficient. Initial testing suggested similar results across both, however.

#Create Monocle object
HSMM = newCellDataSet(as.matrix(data_no_ercc[, good_cells]), phenoData = new("AnnotatedDataFrame", data = phenotype[good_cells, ]), featureData = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(data_no_ercc), col1 = rownames(data_no_ercc), col2 = rownames(data_no_ercc))), expressionFamily=negbinomial.size())
HSMM = estimateSizeFactors(HSMM)
HSMM = estimateDispersions(HSMM)

#Now that we have estimated dispersions, we can replace the Monocle-calculated size factors with our own scran-normalized size factors
pData(HSMM)$Size_Factor = spike_factors

#tSNE Plot
disp_table = dispersionTable(HSMM)
#For clustering, we wanted to use only the most dispersed genes; however, based on how dispersions are calculated in Monocle, the most dispersed genes are inevitably also the lowest expressed. So we used two steps - first, filter out only reasonably highly expressed genes, then selected the top 8000 most dispersed genes. The 8000 number was arbitrary, mostly for visual results - anywhere from 2000-10000 produced similar results, however.
unsup_clustering_genes = subset(disp_table, mean_expression >= 0.1)
unsup_clustering_genes = unsup_clustering_genes[order(unsup_clustering_genes$dispersion_empirical, decreasing = TRUE), ]
HSMM = setOrderingFilter(HSMM, head(unsup_clustering_genes$gene_id, 8000))
plot_ordering_genes(HSMM)
HSMM = reduceDimension(HSMM, max_components = 2, num_dim = 4,
                        reduction_method = 'tSNE', verbose = T)
HSMM = clusterCells(HSMM, num_clusters = 4)
plot_cell_clusters(HSMM, 1, 2, color = "celltype")

#Differential gene expression testing
#We found batch effects (here, lane effects) to be relatively minimal in our samples of interest, so we used a fairly simple model formula string
#We performed differential testing on only the satellite cell groups, so we first filtered out the iPSCs
HSMM_sc = HSMM[, HSMM$celltype != "ipsc"]
HSMM_sc = detectGenes(HSMM_sc, min_expr = 0.1)
diff_gene_test= differentialGeneTest(HSMM_sc, fullModelFormulaStr = "~group")
diff_gene_test = diff_gene_test[order(diff_gene_test$qval), ]
good_genes = rownames(fData(HSMM_sc))[fData(HSMM_sc)$num_cells_expressed > 30] #While we ran testing on everything, we filter out genes with relatively high expression for future steps
#We append on fold changes to further assist with calculation
normalized_data = sweep(exprs(HSMM), 2, pData(HSMM)$Size_Factor, "/")
diff_gene_test$mean_invivo = rowMeans(normalized_data[rownames(diff_gene_test), pData(HSMM)$group == "injected"])
diff_gene_test$mean_invitro = rowMeans(normalized_data[rownames(diff_gene_test), pData(HSMM)$group == "in_vitro"])
diff_gene_test$mean_satcells = rowMeans(normalized_data[rownames(diff_gene_test), pData(HSMM)$group != "esc"])
diff_gene_test$log2fc = log2((diff_gene_test$mean_invivo + 0.001)/(diff_gene_test$mean_invitro + 0.001))
diff_gene_test$good = rownames(diff_gene_test) %in% good_genes
diff_gene_test$diff = diff_gene_test$qval < 0.05
diff_genes = rownames(diff_gene_test)[diff_gene_test$diff == TRUE & diff_gene_test$good == TRUE]

#Plots______________________________________________________________________________________________________________________________

#tSNE Plot
tsne = as.data.frame(t(HSMM@reducedDimA))
tsne$celltype = pData(HSMM)$celltype
ggplot(tsne, aes(x = V1, y = V2, color = celltype)) + geom_point() + scale_color_manual(values=c("red", "orange", "blue", "green4"), labels = c(paste("1 mon in vitro PAX7::GFP+"), "in vitro PAX7::GFP+", "1 mon in vivo PAX7::GFP+", "hESC OCT4::GFP+")) + labs(x = "Component 1", y = "Component 2") + theme_classic() + theme(legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),  panel.border = element_rect(colour = "black", fill=NA), aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(),legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), axis.title = element_text(size = 12)) + ggtitle("t-SNE Plot")

#Normalized RNA Content
pData(HSMM)$normalized_rna = colSums(normalized_data)
ggplot(pData(HSMM), aes(x = celltype, y = normalized_rna, fill = celltype)) + geom_boxplot() + scale_fill_manual(values=c("red", "orange", "blue", "green4"), labels = c(paste("1 mon in vitro PAX7::GFP+"), "in vitro PAX7::GFP+", "1 mon in vivo PAX7::GFP+", "hESC OCT4::GFP+")) + labs(x = "Celltype", y = "Normalized ReNA Content") + theme_classic() + theme(legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),  panel.border = element_rect(colour = "black", fill=NA), aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(),legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), axis.title = element_text(size = 12))

#MA Plot
ggplot(diff_gene_test) + geom_point(data = diff_gene_test[diff_gene_test$qval > 0.05 | diff_gene_test$good == FALSE, ], aes(x = log10(mean_satcells), y = log2fc), color = "black") + geom_point(data = diff_gene_test[diff_gene_test$qval < 0.05 & diff_gene_test$log2fc < 0 & diff_gene_test$good == TRUE, ], aes(x = log10(mean_satcells), y = log2fc) , color = "red") + geom_point(data = diff_gene_test[diff_gene_test$qval < 0.05 & diff_gene_test$log2fc > 0 & diff_gene_test$good == TRUE, ], aes(x = log10(mean_satcells), y = log2fc) , color = "blue") + geom_hline(yintercept = 0, color = "red", linetype = 2, size = 1) + geom_label_repel(data = diff_gene_test[c("NDRG2", "NFIA", "SPRY1", "DUSP1", "NOTCH3", "HEYL"), ], ylim = c(0, NA), nudge_y = 1, aes(label = col1, x = log10(mean_satcells), y = log2fc), box.padding   = 0.35, point.padding = 1, segment.color = 'black') + geom_label_repel(data = diff_gene_test[c("MYOG", "MEF2C", "MYOD1"), ], ylim = c(NA, 0), nudge_y = -1, aes(label = col1, x = log10(mean_satcells), y = log2fc), box.padding   = 0.35, point.padding = 1, segment.color = 'black') +    theme(legend.position = c(0, 1),legend.justification = c(0, 1))+scale_color_manual(values = c("blue","red")) + labs(x = expression(log[10](mean~expression)), y = expression(log[2](fold~change))) + annotate("text", label = "260 genes", x = 1.5, y = 10, color = "blue", size = 8) + annotate("text", label = "11,368 genes", x = 1.5, y = -10, color = "red", size = 8) + ggtitle("MA Plot for in vivo vs in vitro PAX7::GFP+ MPCs") + theme_classic() + theme(legend.title = element_blank(),legend.spacing.y = unit(0, "mm"), panel.border = element_rect(colour = "black", fill=NA), aspect.ratio = 0.6, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5))

#Volcano Plot
ggplot(diff_gene_test) + geom_point(data = diff_gene_test[diff_gene_test$qval > 0.05 | diff_gene_test$good == FALSE, ], aes(x = log2fc, y = -log10(qval)), color = "black") + geom_point(data = diff_gene_test[diff_gene_test$qval < 0.05 & diff_gene_test$log2fc < 0 & diff_gene_test$good == TRUE, ], aes(x = log2fc, y = -log10(qval)) , color = "red") + geom_point(data = diff_gene_test[diff_gene_test$qval < 0.05 & diff_gene_test$log2fc > 0 & diff_gene_test$good == TRUE, ], aes(x = log2fc, y = -log10(qval)) , color = "blue") + geom_hline(yintercept = -log10(0.05), color = "red", linetype = 2, size = 1) + geom_label_repel(data = diff_gene_test[c("NDRG2", "NFIA", "SPRY1", "DUSP1", "NOTCH3", "HEYL"), ], xlim = c(0, NA), nudge_y = 30, aes(label = col1, x = log2fc, y = -log10(qval)), box.padding   = 0.35, point.padding = 1, segment.color = 'black') + geom_label_repel(data = diff_gene_test[c("MYOG", "MEF2C", "MYOD1"), ], xlim = c(NA, 0), nudge_y = -10, aes(label = col1, x = log2fc, y = -log10(qval)), box.padding   = 0.35, point.padding = 1, segment.color = 'black') + theme(legend.position = c(0, 1),legend.justification = c(0, 1))+scale_color_manual(values = c("blue","red")) + labs(x = expression(log[2](fold~change)), y = expression(-log[10](q-value))) + annotate("text", label = "260 genes", x = 7.5, y = 100, color = "blue", size = 8) + annotate("text", label = "11,368 genes", x = -7.5, y = 100, color = "red", size = 8) + ggtitle("Volcano Plot for in vivo vs in vitro PAX7::GFP+ MPCs") + theme_classic() + theme(legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),  panel.border = element_rect(colour = "black", fill=NA), aspect.ratio = 0.6, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(),legend.box.background = element_rect(colour = "black"), plot.title = element_text(hjust = 0.5))

#Heatmaps
#insert genes here
breaks_pval = read.table("breaks_pval.txt", as.is = TRUE)$V1
colors_pval = read.table("colors_pval.txt", as.is = TRUE)$V1 #ps the colors are accidentally backwards; just flip them below
breaks_fc = seq(-100, 0, length=101)
#colors_fc = c(rep("darkgreen", 51), colorRampPalette(c("darkgreen", "white"))(50))
colors_fc = c(rep("gray42", 51), colorRampPalette(c("gray42", "white"))(50))
pheatmap(t(data.frame(row.names = genes, fc = diff_gene_test[genes, ]$log2fc)), show_rownames = FALSE, show_colnames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE , cellheight = 20, cellwidth = 20, color = rev(colors_pval), breaks = breaks_pval)
pheatmap(t(data.frame(row.names = genes, fc = log10(diff_gene_test[genes, ]$qval))), show_rownames = FALSE, show_colnames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE , cellheight = 20, cellwidth = 20, color = colors_fc, breaks = breaks_fc)

#Full Heatmap for Sunny
sunny_genes = c("MYF5", "CD34", "VCAM1", "GLI2", "MEOX2", "ACTA1", "BMP7", "MYL3", "NUPR1", "SOX6", "SOX8", "TGFB2", "ACTC1", "ANKRD1", "CAV1", "CAV3", "GATA4", "MEF2C", "MYOD1", "SOX7", "TBX3", "TNNI1", "TNNI3", "VEGFA", "NOTCH3", "HEYL", "SFRP4", "FZD4", "FZD8", "WNT5A", "WNT5B", "DKK1", "FZD1", "FZD6", "AXIN1", "COL1A1", "COL1A2", "COL4A1", "LAMA2", "LAMA4", "LAMA1", "LAMA3", "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "CCNA1", "CCNB1", "CCND1", "CCNE2", "MCM2", "MCM3", "MCM4", "DUSP1", "DUSP10", "DUSP16", "FOX", "JUN", "JUND", "TP53", "CCL2", "CXCL2", "IL11RA", "CXCL13", "CXCL14", "CCR3", "IL21R")
pheatmap(normalized_data[sunny_genes[sunny_genes %in% diff_genes], sce$celltype != "esc"][, order(sce[, sce$celltype != "esc"]$celltype)], show_rownames = TRUE, show_colnames = FALSE, cluster_rows = TRUE, cluster_cols = FALSE, color = colors_pval, breaks = breaks_pval/5, scale = "row", annotation = phenotype[rownames(phenotype) %in% colnames(sce)[sce$celltype != "esc"], c(4, 5)] )

#Analysis of timepoint data_______________________________________________________________________________________________________

timepoint_no_ercc = read.csv("timepoint_data_no_ercc.csv", as.is = TRUE, row.names = 1)
timepoint_w_ercc = read.csv("timepoint_data_w_ercc.csv", as.is = TRUE, row.names = 1)
timepoint_phenotype = read.csv("timepoint_phenotype.csv", as.is = TRUE, row.names = 1)
colnames(timepoint_phenotype)[1:4] = c("i7barcode", "cellbarcode", "libraryid", "celltype")

#We again prefilter out bad cells, but using slightly different criteria as this experiment was done with lower depth
good_timepoint_cells = rownames(timepoint_phenotype)[timepoint_phenotype$umis > 1000 & timepoint_phenotype$umis < 100000 & timepoint_phenotype$genes > 400 & timepoint_phenotype$percentage_ERCC_reads < 30]

#To compare this dataset to the previous dataset, we will use mnnCorrect as part of scater/scran. However, initial exploration suggested a potential batch effect between the different experimental batches in the timepoint data. Thus, we first split up the timepoint data by batch and normalize each to their respective spike-ins. We then use mnnCorrect, treating the original data as reference data and adding on each new batch one at a time.

sce_tp_1 = SingleCellExperiment(list(counts=as.matrix(timepoint_w_ercc[, rownames(timepoint_phenotype)[timepoint_phenotype$libraryid == 1 & rownames(timepoint_phenotype) %in% good_timepoint_cells]])), colData = timepoint_phenotype[timepoint_phenotype$libraryid == 1 & rownames(timepoint_phenotype) %in% good_timepoint_cells, ])
isSpike(sce_tp_1, "MySpike") = grep("^ERCC", rownames(sce_tp_1))
sce_tp_1 <- computeSpikeFactors(sce_tp_1)
sce_tp_2 = SingleCellExperiment(list(counts=as.matrix(timepoint_w_ercc[, rownames(timepoint_phenotype)[timepoint_phenotype$libraryid == 2 & rownames(timepoint_phenotype) %in% good_timepoint_cells]])), colData = timepoint_phenotype[timepoint_phenotype$libraryid == 2 & rownames(timepoint_phenotype) %in% good_timepoint_cells, ])
isSpike(sce_tp_2, "MySpike") = grep("^ERCC", rownames(sce_tp_2))
sce_tp_2 <- computeSpikeFactors(sce_tp_2)
sce_tp_3 = SingleCellExperiment(list(counts=as.matrix(timepoint_w_ercc[, rownames(timepoint_phenotype)[timepoint_phenotype$libraryid == 3 & rownames(timepoint_phenotype) %in% good_timepoint_cells]])), colData = timepoint_phenotype[timepoint_phenotype$libraryid == 3 & rownames(timepoint_phenotype) %in% good_timepoint_cells, ])
isSpike(sce_tp_3, "MySpike") = grep("^ERCC", rownames(sce_tp_3))
sce_tp_3 <- computeSpikeFactors(sce_tp_3)
sce_tp_4 = SingleCellExperiment(list(counts=as.matrix(timepoint_w_ercc[, rownames(timepoint_phenotype)[timepoint_phenotype$libraryid == 4 & rownames(timepoint_phenotype) %in% good_timepoint_cells]])), colData = timepoint_phenotype[timepoint_phenotype$libraryid == 4 & rownames(timepoint_phenotype) %in% good_timepoint_cells, ])
isSpike(sce_tp_4, "MySpike") = grep("^ERCC", rownames(sce_tp_4))
sce_tp_4 <- computeSpikeFactors(sce_tp_4)
sce_tp_5 = SingleCellExperiment(list(counts=as.matrix(timepoint_w_ercc[, rownames(timepoint_phenotype)[timepoint_phenotype$libraryid == 5 & rownames(timepoint_phenotype) %in% good_timepoint_cells]])), colData = timepoint_phenotype[timepoint_phenotype$libraryid == 5 & rownames(timepoint_phenotype) %in% good_timepoint_cells, ])
isSpike(sce_tp_5, "MySpike") = grep("^ERCC", rownames(sce_tp_5))
sce_tp_5 <- computeSpikeFactors(sce_tp_5)

#normalize all datasets
logcounts(sce) = log(sweep(counts(sce), 2, spike_factors, "/") + 1)
sce_tp_1 = normalize(sce_tp_1)
sce_tp_2 = normalize(sce_tp_2)
sce_tp_3 = normalize(sce_tp_3)
sce_tp_4 = normalize(sce_tp_4)
sce_tp_5 = normalize(sce_tp_5)

#Typically, you cannot use all of the genes for mnnCorrect. The authors generally suggest using genes with highest variance. We tried several options, keeping in mind that the studies have a significant difference in depth, and thus including too many genes would likely result in false clustering due to worse gene capture (exacerbated by the fact that the second study did not have any in vitro cells). We tried two approaches: either selecting the most differentially expressed genes, or the genes with highest expression. In both cases, as long as the first ~2000 or so were selected, we observed similar results. Both options are presented below.

#option 1
unsup = rownames(diff_gene_test)[rownames(diff_gene_test) %in% rownames(sce_tp_1)][1:2000]
#option 2
unsup = unique(c(rownames(diff_gene_test[order(diff_gene_test$mean_invitro, decreasing = TRUE), ])[1:1000], rownames(diff_gene_test[order(diff_gene_test$mean_invivo, decreasing = TRUE), ])[1:1000]))
unsup = unsup[unsup %in% rownames(sce_tp_1)]

#mnnCorrect
original <- list(logcounts(sce)[unsup,],
                 logcounts(sce_tp_1)[unsup,],
                 logcounts(sce_tp_2)[unsup,],
                 logcounts(sce_tp_3)[unsup,],
                 logcounts(sce_tp_4)[unsup,],
                 logcounts(sce_tp_5)[unsup,])
corrected <- do.call(mnnCorrect, c(original, list(k=20, sigma=0.1, cos.norm.out = FALSE)))

#visualize mnnCorrect
corrected_phenotype = as.data.frame(rbind(sce@colData[, 1:4], sce_tp_1@colData[, 1:4], sce_tp_2@colData[, 1:4], sce_tp_3@colData[, 1:4], sce_tp_4@colData[, 1:4], sce_tp_5@colData[, 1:4]))
omat <- do.call(cbind, original)
mat <- do.call(cbind, corrected$corrected)
colnames(omat) = rownames(corrected_phenotype)
colnames(mat) = rownames(corrected_phenotype)
sce_corrected <- SingleCellExperiment(list(original=omat, corrected=mat), colData = corrected_phenotype)
osce <- runTSNE(sce_corrected, exprs_values="original", rand_seed=100)
plotTSNE(osce, colour_by="celltype") + ggtitle("Original")
csce <- runTSNE(sce_corrected, exprs_values="corrected", rand_seed=100)
plotTSNE(csce, colour_by="celltype") + ggtitle("Corrected")

#just a couple of placeholders for plotting
expression = as.data.frame(t(mat))
expression$celltype = corrected_phenotype$celltype
corrected_phenotype_order = rbind(corrected_phenotype[corrected_phenotype$celltype == "1_wk", ], corrected_phenotype[corrected_phenotype$celltype == "2_wk", ], corrected_phenotype[corrected_phenotype$celltype == "3_wk", ], corrected_phenotype[corrected_phenotype$celltype == "4_wk", ], corrected_phenotype[corrected_phenotype$celltype == "esc", ], corrected_phenotype[corrected_phenotype$celltype == "injected", ], corrected_phenotype[corrected_phenotype$celltype == "in_vitro_ctrl", ], corrected_phenotype[corrected_phenotype$celltype == "in_vitro_1_month", ])
corrected_phenotype_order$study = "Study 1"
corrected_phenotype_order[1:332, ]$study = "Study 2"
corrected_phenotype_order$group = "1 wk in vivo PAX7::GFP+"
corrected_phenotype_order[corrected_phenotype_order$celltype == "2_wk", ]$group = "2 wk in vivo PAX7::GFP+"
corrected_phenotype_order[corrected_phenotype_order$celltype == "3_wk", ]$group = "3 wk in vivo PAX7::GFP+"
corrected_phenotype_order[corrected_phenotype_order$celltype %in% c("4_wk", "injected"), ]$group = "4 wk in vivo PAX7::GFP+"
corrected_phenotype_order[corrected_phenotype_order$celltype == "in_vitro_ctrl", ]$group = "in vitro PAX7::GFP+"
corrected_phenotype_order[corrected_phenotype_order$celltype == "in_vitro_1_month", ]$group = "4 wk in vitro PAX7::GFP+"
corrected_phenotype_order[corrected_phenotype_order$celltype == "esc", ]$group = "hESC OCT4::GFP+"
corrected_phenotype_order$group = factor(corrected_phenotype_order$group, levels = c("1 wk in vivo PAX7::GFP+", "2 wk in vivo PAX7::GFP+", "3 wk in vivo PAX7::GFP+", "4 wk in vivo PAX7::GFP+",  "in vitro PAX7::GFP+", "4 wk in vitro PAX7::GFP+", "hESC OCT4::GFP+"))
mat_order= mat[, match(rownames(corrected_phenotype_order), colnames(mat))]

#Plotting heatmap of top 30 upregulated and downregulated genes
up = rownames(diff_gene_test)[diff_gene_test$qval < 0.05 & diff_gene_test$log2fc > 0 & diff_gene_test$good == TRUE]
up = up[up %in% rownames(mat)]
down = rownames(diff_gene_test)[diff_gene_test$qval < 0.05 & diff_gene_test$log2fc < 0 & diff_gene_test$good == TRUE]
down = down[down %in% rownames(mat)]
breaks_comb = unique(c(seq(-5, -2, 1), seq(-2, 2, 0.1), seq(2, 5, 1)))
colors_comb = c(rep("darkblue", 3), colorRampPalette(c("darkblue", 'white'))(21), colorRampPalette(c("white", 'darkred'))(20), rep("darkred", 3))
pheatmap(mat_order[c(up[1:15], down[1:15]), corrected_phenotype_order$celltype != "esc"], cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, scale = "row", annotation = corrected_phenotype_order[corrected_phenotype_order$celltype != "esc", c(6,5)], breaks = breaks_comb, color = colors_comb)

#Calculation of differentially expressed genes between the various timepoints
spike_factors2 = c(sce_tp_1@int_colData$size_factor_MySpike, sce_tp_2@int_colData$size_factor_MySpike, sce_tp_3@int_colData$size_factor_MySpike, sce_tp_4@int_colData$size_factor_MySpike, sce_tp_5@int_colData$size_factor_MySpike)
HSMM_timepoint = newCellDataSet(as.matrix(timepoint_no_ercc[, good_timepoint_cells]), phenoData = new("AnnotatedDataFrame", data = timepoint_phenotype[good_timepoint_cells, ]), featureData = new("AnnotatedDataFrame", data = data.frame(row.names = rownames(timepoint_no_ercc), col1 = rownames(timepoint_no_ercc), col2 = rownames(timepoint_no_ercc))), expressionFamily=negbinomial.size())
HSMM_timepoint = estimateSizeFactors(HSMM_timepoint)
HSMM_timepoint = estimateDispersions(HSMM_timepoint)
pData(HSMM_timepoint)$Size_Factor = spike_factors2
diff_gene_test_tp= differentialGeneTest(HSMM_timepoint, fullModelFormulaStr = "~celltype + libraryid", reducedModelFormulaStr = "~libraryid")
diff_gene_test_tp = diff_gene_test_tp[order(diff_gene_test_tp$qval), ]
HSMM_timepoint = detectGenes(HSMM_timepoint, min_expr = 0.1)
good_timepoint_genes = rownames(fData(HSMM_timepoint))[fData(HSMM_timepoint)$num_cells_expressed > 60] 
diff_gene_test_tp$good = rownames(diff_gene_test_tp) %in% good_timepoint_genes
diff_timepoint_genes = rownames(diff_gene_test_tp)[diff_gene_test_tp$qval < 0.05 & diff_gene_test_tp$good == TRUE]

#SingleCellNet Classification_____________________________________________________________________________________________________

#Here, we wanted to compare the gene expression signatures of our iPSC-derived satellite cells to in vivo counterparts, particularly focusing on quiescent vs. activated profiles. We used the data from Dell'Orso et al. (Development 2019) as our comparison. For simplicity, we have prepared the data and phenotype tables and simply load them in here.
#For comparison of gene expression signatures, we use SingleCellNet (Tan and Cahan, Cell Systems 2019). Please see their vignette for more details. As an aside, because of the low number of phenotypes tested here, and the way the "rands" work in SingleCellNet, the results are somewhat stochastic in terms of absolute classification values. I have tried to capture an "average" run, though perhaps more rigorous aggregation could be used. However, in all iterations, the trends that we observe in the manuscript are maintained, suggesting a common biology.

#load in files
dellorso_data = readMM("dellorso_data.mtx")
rownames(dellorso_data) = read.table("dellorso_gene.txt",as.is = TRUE)$V1
colnames(dellorso_data) = read.table("dellorso_cells.txt",as.is = TRUE)$V1
dellorso_pheno = read.table("dellorso_pheno.csv",as.is = TRUE)

#Define the category "cell" for later SingleCellNet functions
dellorso_pheno$cell = rownames(dellorso_pheno)
phenotype$cell = rownames(phenotype)

#This portion is taken from the SingleCellNet vignette. To define our query data, we use the raw input data (minus ERCCs) without including the ESCs. To define the training group, we used two separate approaches. First, we tested using only the homeostatic MuSC groups (eA and cQ from the original manuscript). This is done below.
stQuery = phenotype[phenotype$celltype != "esc" & rownames(phenotype) %in% good_cells, ]
stQuery$celltype = factor(stQuery$celltype, levels = c("in_vitro_ctrl", "in_vitro_1_month", "injected"))
expQuery = data_no_ercc[, phenotype$celltype != "esc" & rownames(phenotype) %in% good_cells]
stTM = dellorso_pheno[dellorso_pheno$pheno != "injured", ]
stTM = droplevels(stTM)
expTMraw = dellorso_data[, dellorso_pheno$pheno != "injured"]
oTab = utils_loadObject("human_mouse_genes_Jul_24_2018.rda")
aa = csRenameOrth(expQuery, expTMraw, oTab)
expQueryOrth = aa[['expQuery']]
expTrainOrth = aa[['expTrain']]
expTMraw2 = expTrainOrth[,rownames(stTM)]
stList = splitCommon(stTM, ncells=100, dLevel="clean_pheno")
stTrain = stList[[1]]
expTrain = expTMraw2[,rownames(stTrain)]
system.time(class_info2<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 5, nTrees = 1000, nTopGenePairs = 25, dLevel = "clean_pheno", colName_samp = "cell")) #train classifier
stTestList = splitCommon(stList[[2]], ncells=100, dLevel="clean_pheno") 
stTest = stTestList[[1]]
expTest = expTMraw[,rownames(stTest)]
system.time(classRes_val_all2 <- scn_predict(class_info2[['cnProc']], expTest, nrand = 50)) #test on held-out data
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all2, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "clean_pheno", classQuery = "clean_pheno", nRand = 50)
plot_PRs(tm_heldoutassessment) #plot PR curve; note that SingleCellNet also offers a heatmap approach to visualize, but we skip this here
nqRand = 50
system.time(crHS <- scn_predict(class_info2[['cnProc']], expQueryOrth, nrand=nqRand)) #apply to test data
sgrp = as.vector(stQuery$celltype)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS, sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap
sc_violinClass(sampTab = stQuery, classRes = crHS, cellIDCol = "cell", dLevel = "celltype") #plotviolin classification
#I prefer to plot without the rands, and to clean up the group labels; done below
stQuery$nice_label = "1 mon in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "in_vitro_ctrl", ]$nice_label = "in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "injected", ]$nice_label = "1 mon in vivo PAX7::GFP+"
stQuery$nice_label = factor(stQuery$nice_label, levels = c("in vitro PAX7::GFP+", "1 mon in vitro PAX7::GFP+", "1 mon in vivo PAX7::GFP+"))
sc_violinClass(sampTab = stQuery, classRes = crHS[!rownames(crHS) == "rand", ], cellIDCol = "cell", dLevel = "nice_label") +scale_fill_manual(values=c("red", "orange", "royalblue")) + ylab("Classification Score")
sgrp = as.factor(stQuery$nice_label)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS[!rownames(crHS) == "rand", ], sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap

#We now repeat the above, except this type including all of the cells (including the cells from the post injury sample) in the training dataset
stTM = dellorso_pheno
stTM = droplevels(stTM)
expTMraw = dellorso_data
oTab = utils_loadObject("human_mouse_genes_Jul_24_2018.rda")
aa = csRenameOrth(expQuery, expTMraw, oTab)
expQueryOrth = aa[['expQuery']]
expTrainOrth = aa[['expTrain']]
expTMraw2 = expTrainOrth[,rownames(stTM)]
stList = splitCommon(stTM, ncells=100, dLevel="clean_pheno")
stTrain = stList[[1]]
expTrain = expTMraw2[,rownames(stTrain)]
system.time(class_info2<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 10, nTrees = 1000, nTopGenePairs = 25, dLevel = "clean_pheno", colName_samp = "cell")) #train classifier
stTestList = splitCommon(stList[[2]], ncells=100, dLevel="clean_pheno") 
stTest = stTestList[[1]]
expTest = expTMraw[,rownames(stTest)]
system.time(classRes_val_all2 <- scn_predict(class_info2[['cnProc']], expTest, nrand = 50)) #test on held-out data
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all2, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "clean_pheno", classQuery = "clean_pheno", nRand = 50)
plot_PRs(tm_heldoutassessment) #plot PR curve; note that SingleCellNet also offers a heatmap approach to visualize, but we skip this here
nqRand = 50
system.time(crHS <- scn_predict(class_info2[['cnProc']], expQueryOrth, nrand=nqRand)) #apply to test data
sgrp = as.vector(stQuery$celltype)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS, sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap
sc_violinClass(sampTab = stQuery, classRes = crHS, cellIDCol = "cell", dLevel = "celltype") #plotviolin classification
#Again, with cleaner plotting
stQuery$nice_label = "1 mon in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "in_vitro_ctrl", ]$nice_label = "in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "injected", ]$nice_label = "1 mon in vivo PAX7::GFP+"
stQuery$nice_label = factor(stQuery$nice_label, levels = c("in vitro PAX7::GFP+", "1 mon in vitro PAX7::GFP+", "1 mon in vivo PAX7::GFP+"))
sc_violinClass(sampTab = stQuery, classRes = crHS[!rownames(crHS) == "rand", ], cellIDCol = "cell", dLevel = "nice_label") +scale_fill_manual(values=c("red", "orange", "royalblue")) + ylab("Classification Score")
sgrp = as.factor(stQuery$nice_label)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS[!rownames(crHS) == "rand", ], sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap

#Validation of single cell species by mapping_____________________________________________________________________________________
#To validate that our cells (particularly the injected cells) are human in origin, we mapped the samples to both human and mouse genomes. We quantified the number of reads mapped to each species to validate that cells were of human origin.

mouse_mapping = read.table("mouse_mapping.txt", as.is = TRUE)
human_mapping = read.table("human_mapping.txt", as.is = TRUE)
mouse_mapping = mouse_mapping[good_cells, ]
human_mapping = human_mapping[good_cells, ]
mouse_mapping$cellname = rownames(mouse_mapping)
human_mapping$cellname = rownames(human_mapping)
mouse_mapping = melt(mouse_mapping)
human_mapping = melt(human_mapping)
mouse_mapping$species = "Mouse"
human_mapping$species = "Human"
mapping = rbind(human_mapping, mouse_mapping)
mapping$celltype = phenotype[mapping$cellname, ]$celltype
mapping$nice_label = "1 mon in vitro PAX7::GFP+"
mapping[mapping$celltype == "in_vitro_ctrl", ]$nice_label = "in vitro PAX7::GFP+"
mapping[mapping$celltype == "injected", ]$nice_label = "1 mon in vivo PAX7::GFP+"
mapping[mapping$celltype == "esc", ]$nice_label = "hESC OCT4::GFP+"
stQuery$nice_label = factor(stQuery$nice_label, levels = c("in vitro PAX7::GFP+", "1 mon in vitro PAX7::GFP+", "1 mon in vivo PAX7::GFP+", "hESC OCT4::GFP+"))

ggplot(mapping, aes(x = variable, y = value * 100, fill = species)) + geom_boxplot()  + labs(x = "Read Classification", y = "% Reads per Cell") + theme_classic() + theme(legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),  panel.border = element_rect(colour = "black", fill=NA), aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(),legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), axis.title = element_text(size = 14))
ggplot(mapping, aes(x = variable, y = value * 100, fill = species)) + geom_boxplot()  + labs(x = "Read Classification", y = "% Reads per Cell") + theme_classic() + theme(legend.title = element_blank(),legend.spacing.y = unit(0, "mm"),  panel.border = element_rect(colour = "black", fill=NA), axis.text = element_text(colour = 1, size = 12), legend.background = element_blank(),legend.box.background = element_blank(), plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 12), axis.title = element_text(size = 14)) + facet_wrap(~nice_label)

#Comparison to Charville data
#_________________________________________________________________________________________________________________________________
#We cross-compare our analysis to the data gathered in Charville et al. (2015). Note that this data is also used to colour the genes in the heatmap; we selected log2fc threshold of +/- 1.5.

charville_data = read.table("~/Documents/Research/SunnySun/FinalManuscript/charville_data.txt", as.is = TRUE, row.names = 1, header = TRUE)
charville_data$log2fc = log2((charville_data$Average.QSC + 1)/(charville_data$Average.ASC + 1))
#for simplicity, filter diff_gene_test to only "good" genes using our previous threshold, and select only the ones that are also in the Charville data
diff_comparison = diff_gene_test[diff_gene_test$good == TRUE, ]
diff_comparison = diff_comparison[rownames(diff_comparison) %in% rownames(charville_data), ]
diff_comparison$log2fc_charville = charville_data[rownames(diff_comparison), ]$log2fc
diff_comparison$diff_charville = abs(diff_comparison$log2fc_charville) >= 1.5
diff_comparison$diff_direction = "None"
diff_comparison[diff_comparison$diff == TRUE & diff_comparison$log2fc > 0, ]$diff_direction = "Quiescent"
diff_comparison[diff_comparison$diff == TRUE & diff_comparison$log2fc < 0, ]$diff_direction = "Activated"
diff_comparison$diff_charville_direction = "None"
diff_comparison[diff_comparison$log2fc_charville >= 1.5, ]$diff_charville_direction = "Quiescent"
diff_comparison[diff_comparison$log2fc_charville <= -1.5, ]$diff_charville_direction = "Activated"

#SingleCellNet Classification against staged SMPC/SCs______________________________________________________________________________

#Here, we wanted to compare the gene expression signatures of our iPSC-derived satellite cells to in vivo counterparts, particularly focusing on various developmental stages. We used the data from Xi et al. (Cell Stem Cell 2020) as our comparison - this data contains, among other datasets, skeletal myogenic progenitors and stem cells from embryonic, fetal, and postnatal human muscle. The authors classify these, based on diffusion mapping output, into five stages, from least to most mature. For simplicity, we have prepared the data and phenotype tables and simply load them in here.
#For comparison of gene expression signatures, we use SingleCellNet (Tan and Cahan, Cell Systems 2019). Please see their vignette for more details, and see the caveats in the previous SingleCellNet comparison section above.

#load in files
xi_data = read.table("xi_data.tsv.gz", as.is = TRUE, row.names = 1, header = TRUE)
xi_pheno = read.table("xi_pheno.tsv",as.is = TRUE, row.names = 1, header = TRUE)

#Define the category "cell" for later SingleCellNet functions
xi_pheno$cell = rownames(xi_pheno)
phenotype$cell = rownames(phenotype)

#This portion is taken from the SingleCellNet vignette. To define our query data, we use the raw input data (minus ERCCs) without including the ESCs. We defined the training group from Fig. 4A of Xi et al.
stQuery = phenotype[phenotype$celltype != "esc" & rownames(phenotype) %in% good_cells, ]
stQuery$celltype = factor(stQuery$celltype, levels = c("in_vitro_ctrl", "in_vitro_1_month", "injected"))
expQuery = data_no_ercc[, phenotype$celltype != "esc" & rownames(phenotype) %in% good_cells]
stTM = xi_pheno
stTM = droplevels(stTM)
expTMraw = xi_data
mutual_genes = intersect(rownames(expQuery), rownames(expTMraw))
expQuery = as.matrix(expQuery[mutual_genes, ])
expTMraw = as.matrix(expTMraw[mutual_genes, ])
stList = splitCommon(stTM, ncells=100, dLevel="Dev.Stage")
stTrain = stList[[1]]
expTrain = expTMraw[,rownames(stTrain)]
system.time(class_info2<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 25, nRand = 5, nTrees = 1000, nTopGenePairs = 25, dLevel = "Dev.Stage", colName_samp = "cell")) #train classifier
stTestList = splitCommon(stList[[2]], ncells=100, dLevel="Dev.Stage") 
stTest = stTestList[[1]]
expTest = expTMraw[,rownames(stTest)]
system.time(classRes_val_all2 <- scn_predict(class_info2[['cnProc']], expTest, nrand = 50)) #test on held-out data
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all2, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "Dev.Stage", classQuery = "Dev.Stage", nRand = 50)
plot_PRs(tm_heldoutassessment) #plot PR curve; note that SingleCellNet also offers a heatmap approach to visualize, but we skip this here
nqRand = 50
system.time(crHS <- scn_predict(class_info2[['cnProc']], expQuery, nrand=nqRand)) #apply to test data
sgrp = as.vector(stQuery$celltype)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS, sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap
sc_violinClass(sampTab = stQuery, classRes = crHS, cellIDCol = "cell", dLevel = "celltype") #plotviolin classification
#I prefer to plot without the rands, and to clean up the group labels; done below
stQuery$nice_label = "1 mon in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "in_vitro_ctrl", ]$nice_label = "in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "injected", ]$nice_label = "1 mon in vivo PAX7::GFP+"
stQuery$nice_label = factor(stQuery$nice_label, levels = c("in vitro PAX7::GFP+", "1 mon in vitro PAX7::GFP+", "1 mon in vivo PAX7::GFP+"))
sc_violinClass(sampTab = stQuery, classRes = crHS[!rownames(crHS) == "rand", ], cellIDCol = "cell", dLevel = "nice_label") +scale_fill_manual(values=c("red", "orange", "royalblue")) + ylab("Classification Score")
sgrp = as.factor(stQuery$nice_label)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS[!rownames(crHS) == "rand", ], sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap
#I like an attribution plot for this as well
colnames(stQuery)[12] = "group2" #had to do this to let the SCN code work, oops
plot_attr(crHS[ , 1:376], stQuery, nrand=0, sid="cell", dLevel="nice_label") + xlab("Group") + ylab("Count")


stQuery = rbind(phenotype[phenotype$celltype != "esc" & rownames(phenotype) %in% good_cells, 1:4], timepoint_phenotype_good[, 1:4])
stQuery$celltype = factor(stQuery$celltype, levels = c("in_vitro_ctrl", "in_vitro_1_month", "1_wk", "2_wk", "3_wk", "4_wk", "injected"))
common_exp_genes = intersect(rownames(exprs(HSMM_sc)), rownames(timepoint_no_ercc))
expQuery = cbind(exprs(HSMM_sc)[common_exp_genes, ], timepoint_no_ercc[common_exp_genes, good_timepoint_cells])
stTM = xi_pheno
stTM = droplevels(stTM)
expTMraw = xi_data
mutual_genes = intersect(rownames(expQuery), rownames(expTMraw))
expQuery = as.matrix(expQuery[mutual_genes, ])
expTMraw = as.matrix(expTMraw[mutual_genes, ])
stList = splitCommon(stTM, ncells=100, dLevel="Dev.Stage")
stTrain = stList[[1]]
expTrain = expTMraw[,rownames(stTrain)]
system.time(class_info2<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 25, nRand = 5, nTrees = 1000, nTopGenePairs = 25, dLevel = "Dev.Stage", colName_samp = "cell")) #train classifier
stTestList = splitCommon(stList[[2]], ncells=100, dLevel="Dev.Stage") 
stTest = stTestList[[1]]
expTest = expTMraw[,rownames(stTest)]
system.time(classRes_val_all2 <- scn_predict(class_info2[['cnProc']], expTest, nrand = 50)) #test on held-out data
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all2, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "Dev.Stage", classQuery = "Dev.Stage", nRand = 50)
plot_PRs(tm_heldoutassessment) #plot PR curve; note that SingleCellNet also offers a heatmap approach to visualize, but we skip this here
nqRand = 50
system.time(crHS <- scn_predict(class_info2[['cnProc']], expQuery, nrand=0)) #apply to test data
colnames(crHS) = rownames(stQuery)
stQuery$cell = rownames(stQuery)
sgrp = as.vector(stQuery$celltype)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS, sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap
sc_violinClass(sampTab = stQuery, classRes = crHS, cellIDCol = "cell", dLevel = "celltype") #plotviolin classification
#I prefer to plot without the rands, and to clean up the group labels; done below
stQuery$nice_label = "1 mon in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "in_vitro_ctrl", ]$nice_label = "in vitro PAX7::GFP+"
stQuery[stQuery$celltype == "injected", ]$nice_label = "1 mon in vivo PAX7::GFP+"
stQuery[stQuery$celltype == "1_wk", ]$nice_label = "1 wk in vivo PAX7::GFP+"
stQuery[stQuery$celltype == "2_wk", ]$nice_label = "2 wk in vivo PAX7::GFP+"
stQuery[stQuery$celltype == "3_wk", ]$nice_label = "3 wk in vivo PAX7::GFP+"
stQuery[stQuery$celltype == "4_wk", ]$nice_label = "1 mon in vivo PAX7::GFP+"
stQuery$nice_label = factor(stQuery$nice_label, levels = c("in vitro PAX7::GFP+", "1 mon in vitro PAX7::GFP+", "1 wk in vivo PAX7::GFP+", "2 wk in vivo PAX7::GFP+", "3 wk in vivo PAX7::GFP+", "1 mon in vivo PAX7::GFP+"))
sc_violinClass(sampTab = stQuery, classRes = crHS[!rownames(crHS) == "rand", ], cellIDCol = "cell", dLevel = "nice_label") + scale_fill_manual(values=c("red", "orange", "skyblue", "deepskyblue", "dodgerblue", "royalblue")) + ylab("Classification Score")
sgrp = as.factor(stQuery$nice_label)
names(sgrp) = as.vector(stQuery$cell)
sc_hmClass(crHS[!rownames(crHS) == "rand", ], sgrp, max=5000, isBig=TRUE, cCol=F, font=8) #plot heatmap
#I like an attribution plot for this as well
#colnames(stQuery)[12] = "group2" #had to do this to let the SCN code work, oops
plot_attr(crHS[ , 1:708], stQuery, nrand=0, sid="cell", dLevel="nice_label") + xlab("Group") + ylab("Count")

#For answering some reviewer questions, I made a quick heatmap of the SCN classifier genes for the reference and the query data. I did it as averaged pseudo-bulk just to make the heatmap easier to see.
classifier_genes = c(class_info2$cgenes_list$Stage_1, class_info2$cgenes_list$Stage_2, class_info2$cgenes_list$Stage_3, class_info2$cgenes_list$Stage_4, class_info2$cgenes_list$Stage_5)
xi_data_agg = data.frame(row.names = rownames(xi_data), Stage_1 = rowMeans(xi_data[, xi_pheno$Dev.Stage == "Stage_1"]), Stage_2 = rowMeans(xi_data[, xi_pheno$Dev.Stage == "Stage_2"]), Stage_3 = rowMeans(xi_data[, xi_pheno$Dev.Stage == "Stage_3"]), Stage_4 = rowMeans(xi_data[, xi_pheno$Dev.Stage == "Stage_4"]), Stage_5 = rowMeans(xi_data[, xi_pheno$Dev.Stage == "Stage_5"]))
classifier_heatmap = pheatmap(log(xi_data_agg[class_info2$cnProc$cgenes, ]+1), cluster_rows = TRUE, cluster_cols = FALSE, scale = "row", show_colnames = FALSE)
classifier_heatmap_genes = cutree(classifier_heatmap$tree_row, k = 2)
pheatmap(log(xi_data_agg[class_info2$cnProc$cgenes[c(classifier_heatmap$tree_row$order[44:114], classifier_heatmap$tree_row$order[1:43])], ]+1), cluster_rows = FALSE, cluster_cols = FALSE, scale = "row", show_colnames = TRUE)
data_agg = data.frame(row.names = rownames(expQuery), in_vitro_ctrl = rowMeans(expQuery[, stQuery$celltype == "in_vitro_ctrl"]), in_vitro_1_month = rowMeans(expQuery[, stQuery$celltype == "in_vitro_1_month"]), injected = rowMeans(expQuery[, stQuery$celltype == "injected"]))
pheatmap(log(data_agg[class_info2$cnProc$cgenes[c(classifier_heatmap$tree_row$order[44:114], classifier_heatmap$tree_row$order[1:43])], ]+1), cluster_rows = FALSE, cluster_cols = FALSE, scale = "row", show_colnames = TRUE)
