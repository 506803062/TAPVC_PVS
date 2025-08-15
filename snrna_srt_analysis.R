######################### TAPVC PVS #########################

############## SnRNA-seq ##############

library(dplyr)
library(Seurat)
library(patchwork)
set.seed(123)
library(harmony)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(SCENIC)
library(monocle)
library(CytoTRACE2)
library(AUCell)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(CellChat)
library(iTALK)

#Dimensionality reduction and annotation

## Import of original data

counts_1 <- Read10X(data.dir = "path/PVS/rawdata/H1/filter_matrix", gene.column = 1)
counts_2 <- Read10X(data.dir = "path/PVS/rawdata/H2/filter_matrix", gene.column = 1)
counts_3 <- Read10X(data.dir = "path/PVS/rawdata/H3/filter_matrix", gene.column = 1)
counts_4 <- Read10X(data.dir = "path/PVS/rawdata/PVS1/filter_matrix", gene.column = 1)
counts_5 <- Read10X(data.dir = "path/PVS/rawdata/PVS2/filter_matrix", gene.column = 1)
counts_6 <- Read10X(data.dir = "path/PVS/rawdata/TAPVC1/filter_matrix", gene.column = 1)
counts_7 <- Read10X(data.dir = "path/PVS/rawdata/TAPVC2/filter_matrix", gene.column = 1)

scRNA_1 <- CreateSeuratObject(counts_1, project = "H1", min.cells = 1)
scRNA_2 <- CreateSeuratObject(counts_2, project = "H2", min.cells = 1)
scRNA_3 <- CreateSeuratObject(counts_3, project = "H3", min.cells = 1)
scRNA_4 <- CreateSeuratObject(counts_4, project = "PVS1", min.cells = 1)
scRNA_5 <- CreateSeuratObject(counts_5, project = "PVS2", min.cells = 1)
scRNA_6 <- CreateSeuratObject(counts_6, project = "TAPVC1", min.cells = 1)
scRNA_7 <- CreateSeuratObject(counts_7, project = "TAPVC2", min.cells = 1)

scRNA <- merge(scRNA_1, y=c(scRNA_2, scRNA_3,scRNA_4, scRNA_5, 
                                   scRNA_6, scRNA_7))
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
sample.ids <- c("H1", "H2", "H3", "PVS1", "PVS2", "TAPVC1", "TAPVC2")
group.ids <- c(rep("H", 3), rep("PVS", 2), rep("TAPVC", 2))
scRNA@meta.data$group = plyr::mapvalues(x = scRNA@meta.data[,"orig.ident"], from = sample.ids, to = group.ids)

mycolour <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
              '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
              '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
              '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
              '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
              '#968175')
mycolour2 <- c('#af2337', '#ecc342', '#2967a0', '#2f3c28', '#96b437',
               '#da93ab','#e58932', '#80598f', '#7e331f', '#3b855a',
               '#c0b286', '#a9c9ed', '#ec977f', '#848482', '#604628',
               '#d26034', '#a64c6b', '#dbd245', '#eba83b', '#5d5092',
               '#222222', '#f2f3f4')

##QC

scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)

pdf(file = "path/PVS/image/QC.pdf", width = 6, height = 7)
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, pt.size = 0, cols = mycolour)
dev.off()

##Data standardization
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)

##PCA
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
VizDimLoadings(scRNA, dims = 1:2, reduction = "pca")
scRNA <- JackStraw(scRNA, num.replicate = 100, dims = 50)
scRNA <- ScoreJackStraw(scRNA, dims = 1:50)
JackStrawPlot(scRNA, dims = 1:50)
ElbowPlot(scRNA, ndims = 50)
pdf(file = "path/PVS/image/All_PCA.pdf", width = 5.5, height = 4)
DimPlot(scRNA, reduction = "pca", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

##Eliminate batch effect

scRNA <- scRNA %>% RunHarmony("orig.ident", plot_convergence = TRUE)
pdf(file = "path/PVS/image_3/All_harmony.pdf", width = 5.5, height = 4)
DimPlot(scRNA, reduction = "harmony") + theme_bw() + theme(panel.grid=element_blank())
dev.off()

scRNA <- JackStraw(scRNA, num.replicate = 100, dims = 50)
scRNA <- ScoreJackStraw(scRNA, dims = 1:50)
JackStrawPlot(scRNA, dims = 1:50)
ElbowPlot(scRNA, ndims = 50)

## UMAP

scRNA <- scRNA %>% RunUMAP(reduction = "harmony",dims = 1:50) %>%
  FindNeighbors(reduction = "harmony",dims = 1:50) %>%
  FindClusters(resolution = 0.9) %>% identity()
  
## Annotation of cell population

current.cluster.ids <- c(0:34)
new.cluster.ids <- c("VSMC", "FB", "VSMC", "EC", "VSMC", "FB", "FB", "Macrophage", "MyoFB", "FB", "EC", "EPI",  "EC", 
                     "EC", "VSMC", "EC", "EC", "VSMC", "FB", "VSMC", "MyoFB", "MyoFB", "EC", "OLG", "FB", "T/NK","EC", 
                     "EC", "VSMC", "EC", "VSMC", "VSMC", "Mast", "Mix", "EPI")
scRNA@meta.data$celltype = plyr::mapvalues(x = scRNA@meta.data[,"RNA_snn_res.0.9"], from = current.cluster.ids, to = new.cluster.ids)

scRNA$celltype <- factor(scRNA$celltype,level = c ("VSMC", "FB", "Macrophage", "T/NK", "EC", "Mast", "MyoFB","OLG", "EPI","Mix"))
scRNA$group <- factor(scRNA$group,level = c ("H", "TAPVC", "PVS"))
scRNA$orig.ident <- factor(scRNA$orig.ident,level = c ("H1", "H2", "H3", "TAPVC1", "TAPVC2", "PVS1", "PVS2"))

scRNA <- subset(scRNA, celltype %in% c("VSMC", "FB", "Macrophage", "T/NK", "EC", "Mast", "MyoFB", "OLG", "EPI"))
scRNA$celltype <- factor(scRNA$celltype,level = c ("VSMC", "FB", "Macrophage", "T/NK", "EC", "Mast", "MyoFB", "OLG", "EPI"))

Idents(scRNA) = "celltype"

features <- c("ACTA2", "MYH11", "LUM", "DCN", "CD68", "AIF1", "CD3D", "CD3E", "KLRD1", "NKG7", "PECAM1", "VWF", "TPSAB1", "TPSB2", "PLP1", "NRXN1", "KRT19", "MSLN")

pdf(file = "path/PVS/image/Plot_markergene_0.9_celltype.pdf", width = 6.3, height = 2.6)
DotPlot(scRNA, feature = features, group.by = "celltype", cols = c("white", '#990101')) + theme_bw() + theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(file = "path/PVS/image/All_harmony_UMAP_cluster_0.9_celltype_last.pdf", width = 5.5, height = 4)
DimPlot(scRNA, reduction = "umap", group.by = "celltype", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/PVS/image/All_harmony_TSNE_cluster_0.9_celltype_last.pdf", width = 5.5, height = 4)
DimPlot(scRNA, reduction = "tsne", group.by = "celltype", cols = mycolour2) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

##Analyze the differentially expressed genes in cell populations
Idents(scRNA) = "celltype"
scRNA <- JoinLayers(scRNA)
diff.wilcox = FindAllMarkers(scRNA, min.pct = 0.25, logfc.threshold = 0.25)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(diff.wilcox, "path/PVS/result/diff_genes_wilcox_celltype.csv")
write.csv(top20, "path/PVS/result/top20_diff_genes_wilcox_celltype.csv")

###The proportion of cell groups
df <- table(scRNA$celltype,scRNA$group) %>% melt()
colnames(df) <- c("Cluster","Group","Number")
df$Cluster <- factor(df$Cluster,level = c ("VSMC", "FB", "Macrophage", "T/NK", "EC", "Mast", "MyoFB", "OLG", "EPI"))
df$Group <- factor(df$Group,level = c ("H", "TAPVC", "PVS"))
pdf(file = "path/PVS/image/celltype_group_percent.pdf", width = 4.5, height = 3)
ggplot(data = df, aes(x =Number, y = Group, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8, position= "fill")+
  scale_fill_manual(values=mycolour2) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
dev.off()


########### MyoFB cell subcluster ###########

MyoFBcell <- subset(scRNA, celltype %in% c("MyoFB"))

#sudo-bulk level

##Analyze the differential genes of MyoFB in different groups 
diff_MyoFB_cell_HvsTAPVC <- FindMarkers(MyoFBcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "group",ident.1 = "H", ident.2 = "TAPVC")
diff_MyoFB_cell_HvsPVS <- FindMarkers(MyoFBcell, min.pct = 0.25, logfc.threshold = 0.25, group.by = "group",ident.1 = "H", ident.2 = "PVS")
write.csv(diff_MyoFB_cell_HvsTAPVC, "path/PVS/result/MyoFBcell_diff_genes_wilcox_H_vs_TAPVC.csv", row.names = T)
write.csv(diff_MyoFB_cell_HvsPVS, "path/PVS/result/MyoFBcell_diff_genes_wilcox_H_vs_PVS.csv", row.names = T)

## H vs TAPVC
library(DOSE)
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)

##Draw the H_vs_TAPVC volcano map
diff_MyoFB_cell_HvsTAPVC[which(diff_MyoFB_cell_HvsTAPVC$avg_log2FC  >= 0.25 & diff_MyoFB_cell_HvsTAPVC$p_val < 0.05),'sig'] <- 'up'
diff_MyoFB_cell_HvsTAPVC[which(diff_MyoFB_cell_HvsTAPVC$avg_log2FC  <= -0.25 & diff_MyoFB_cell_HvsTAPVC$p_val < 0.05),'sig'] <- 'down'
diff_MyoFB_cell_HvsTAPVC$symbol <- rownames(diff_MyoFB_cell_HvsTAPVC)

##Screen the top3 genes with significant upregulation
up_data <- filter(diff_MyoFB_cell_HvsTAPVC, sig == 'up') %>% 
  distinct(symbol, .keep_all = TRUE) %>%            
  top_n(3, avg_log2FC)                        

##Screen the genes with the top3 significance in the down-regulation
down_data <- filter(diff_MyoFB_cell_HvsTAPVC, sig == 'down') %>%  
  distinct(symbol, .keep_all = TRUE) %>%                  
  top_n(-3, avg_log2FC) 

library(ggrepel)
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_H_vs_TAPVC_volcano.pdf", width = 4, height = 4)
ggplot(diff_MyoFB_cell_HvsTAPVC, aes(avg_log2FC , -log10(p_val), col = sig)) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_color_manual(values = c("#00314f", "#93292d")) +
  labs(x="log2(FoldChange)",y="-log10 (p_val)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.25, 0.25), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  geom_text_repel(data = up_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) +
  geom_text_repel(data = down_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) 
dev.off()

go_res <- enrichGO(
  rownames(diff_MyoFB_cell_HvsTAPVC),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)
goresult <- data.frame(go_res@result)

go_select <- goresult %>% 
  filter(ONTOLOGY == "BP" & pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  head(10)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

pdf(file = "path/PVS/image/MyoFB/diff_MyoFB_HvsTAPVC_GO.pdf", width = 4, height = 4)
ggplot(go_select, aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", fill = "#fee4b7", width = 0.8) +
  geom_text(aes(x = 0.2, label = Description), hjust = 0, vjust = 1, size = 4.5) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.75),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "-log10(pvalue)", y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = 'off')
dev.off()

## H vs PVS

diff_MyoFB_cell_HvsPVS[which(diff_MyoFB_cell_HvsPVS$avg_log2FC  >= 0.25 & diff_MyoFB_cell_HvsPVS$p_val < 0.05),'sig'] <- 'up'
diff_MyoFB_cell_HvsPVS[which(diff_MyoFB_cell_HvsPVS$avg_log2FC  <= -0.25 & diff_MyoFB_cell_HvsPVS$p_val < 0.05),'sig'] <- 'down'
diff_MyoFB_cell_HvsPVS$symbol <- rownames(diff_MyoFB_cell_HvsPVS)

up_data <- filter(diff_MyoFB_cell_HvsPVS, sig == 'up') %>% 
  distinct(symbol, .keep_all = TRUE) %>%               
  top_n(3, avg_log2FC)                          

down_data <- filter(diff_MyoFB_cell_HvsPVS, sig == 'down') %>% 
  distinct(symbol, .keep_all = TRUE) %>%                  
  top_n(-3, avg_log2FC) 

library(ggrepel)
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_H_vs_PVS_volcano.pdf", width = 4, height = 4)
ggplot(diff_MyoFB_cell_HvsPVS, aes(avg_log2FC , -log10(p_val), col = sig)) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_color_manual(values = c("#00314f", "#93292d")) +
  labs(x="log2(FoldChange)",y="-log10 (p_val)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-0.25, 0.25), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) +
  geom_text_repel(data = up_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) + 
  geom_text_repel(data = down_data, aes(x = avg_log2FC, y = -log10(p_val), label = symbol)) 
dev.off()


go_res <- enrichGO(
  rownames(diff_MyoFB_cell_HvsPVS),
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)
goresult <- data.frame(go_res@result)

go_select <- goresult %>% 
  filter(ONTOLOGY == "BP" & pvalue < 0.05) %>% 
  arrange(pvalue) %>% 
  head(10)  %>% 
  dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))

pdf(file = "path/PVS/image/MyoFB/diff_MyoFB_HvsPVS_GO.pdf", width = 4, height = 4)
ggplot(go_select, aes(x = -log10(pvalue), y = Description)) +
  geom_bar(stat = "identity", fill = "#aedcef", width = 0.8) +
  geom_text(aes(x = 0.2, label = Description), hjust = 0, vjust = 1, size = 4.5) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    axis.line = element_line(linewidth = 0.75, color = "black"),
    axis.ticks.x = element_line(linewidth = 0.75),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, color = "black"),
    panel.grid = element_blank()
  ) +
  labs(x = "-log10(pvalue)", y = "") +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(clip = 'off')
dev.off()

#Single cell level

MyoFBcell <- FindVariableFeatures(MyoFBcell, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MyoFBcell)
MyoFBcell <- ScaleData(MyoFBcell, features = all.genes)

##PCA
MyoFBcell <- RunPCA(MyoFBcell, features = VariableFeatures(object = MyoFBcell))
VizDimLoadings(MyoFBcell, dims = 1:2, reduction = "pca")
MyoFBcell <- JackStraw(MyoFBcell, num.replicate = 100, dims = 50)
MyoFBcell <- ScoreJackStraw(MyoFBcell, dims = 1:50)
p1 <- JackStrawPlot(MyoFBcell, dims = 1:50)
p2 <- ElbowPlot(MyoFBcell, ndims = 50)
pdf(file = "path/PVS/image/MyoFB/ElbowPlot.pdf", width = 13, height = 3)
p1+p2
dev.off()

### UMAP
MyoFBcell$group <- factor(MyoFBcell$group,level = c ("H", "TAPVC", "PVS"))
MyoFBcell <- MyoFBcell %>% RunUMAP(reduction = "harmony",dims = 1:30) %>%
  FindNeighbors(reduction = "harmony",dims = 1:30) %>%
  FindClusters(resolution = 0.05) %>% identity()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_UMAP_cluster_r0.1.pdf", width = 5, height = 4)
DimPlot(MyoFBcell, reduction = "umap", cols = mycolour, label = TRUE) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_UMAP_sample.pdf", width = 5.5, height = 4)
DimPlot(MyoFBcell, reduction = "umap", group.by = "orig.ident", cols = mycolour) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_UMAP_group.pdf", width = 5.5, height = 4)
DimPlot(MyoFBcell, reduction = "umap", group.by = "group", cols = c('#2967a0', '#ecc342', '#af2337')) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

current.cluster.ids <- c(0:2)
new.cluster.ids <- c("MyoFB1", "MyoFB2", "MyoFB3")
MyoFBcell@meta.data$subcelltype = plyr::mapvalues(x = MyoFBcell@meta.data[,"RNA_snn_res.0.05"], from = current.cluster.ids, to = new.cluster.ids)

Idents(MyoFBcell) = "subcelltype"

##Display cell group 
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_harmony_UMAP_subcelltype.pdf", width = 5, height = 4)
DimPlot(MyoFBcell, reduction = "umap", group.by = "subcelltype", cols = c("#D51F26", "#00A08A", "#F2AD00")) + theme_bw() + theme(panel.grid=element_blank())
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_harmony_TSNE_subcelltype.pdf", width = 5, height = 4)
DimPlot(MyoFBcell, reduction = "tsne", group.by = "subcelltype", cols = c("#D51F26", "#00A08A", "#F2AD00")) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

##Analyze the differentially expressed genes in cell populations
Idents(MyoFBcell) = "subcelltype"
diff.wilcox = FindAllMarkers(MyoFBcell, min.pct = 0.25, logfc.threshold = 0.25)
all.markers = diff.wilcox %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
write.csv(all.markers, "path/PVS/result/MyoFBcell_diff_genes_wilcox_subcelltype.csv", row.names = F)

markergene <- c("RPS23","RPL26","ATP5F1E","RPL27","COX7C",
                "CHRM2", "RBM20","MLIP","DAPK2","PPARGC1A",
                "GPC6","NTRK3","ADGRL3","SGCZ","KCNQ5")

pdf(file = "path/PVS/image/MyoFB/Plot_markergene_MyoFB_celltype_2.pdf", width = 7, height = 2)
DotPlot(MyoFBcell, feature = markergene, cols = c("white", '#990101')) + RotatedAxis()+ theme_bw() + theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

###The proportion of subpopulations in the cell population
df <- table(MyoFBcell$subcelltype,MyoFBcell$group) %>% melt()
colnames(df) <- c("Cluster","Group","Number")
df$Cluster <- factor(df$Cluster,level = c ("MyoFB1", "MyoFB2", "MyoFB3"))
df$Group <- factor(df$Group,level = c ("H", "TAPVC", "PVS"))
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_celltype_group_percent.pdf", width = 5, height = 3)
ggplot(data = df, aes(x =Number, y = Group, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c("#D51F26", "#00A08A", "#F2AD00")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )
dev.off()

# GO analysis

top410 = all.markers %>% group_by(cluster) %>% top_n(n = 410, wt = avg_log2FC)
group <- data.frame(gene=top410$gene,
                    group=top410$cluster)

Gene_ID <- bitr(top410$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

data_GO_sim <- clusterProfiler::simplify(data_GO, 
                                         cutoff=0.7, 
                                         by="p.adjust", 
                                         select_fun=min)

dotplot(data_GO_sim, color = "p.adjust", size = "Count", showCategory=3, font.size = 5)+ 
  scale_color_gradientn(colours = c('red','grey'))

go_results <- setReadable(data_GO_sim, 
                          OrgDb = "org.Hs.eg.db", 
                          keyType = "ENTREZID")

go_table <- go_results@compareClusterResult
table(go_table$Cluster)
selected_cell_types <- c("MyoFB1", "MyoFB2", "MyoFB3")
go_table <- go_table[go_table$Cluster %in% selected_cell_types, ]

go_table <- go_table %>% 
  group_by(Cluster) %>% 
  top_n(n = 3, wt = -qvalue) 
table(go_table$Cluster)

go_table$Cluster <- factor(go_table$Cluster, levels = selected_cell_types)

go_table <- go_table[order(go_table$Cluster), ]

go_table$Description <- factor(go_table$Description, levels = go_table$Description)

go_table$geneID <- sapply(strsplit(go_table$geneID, "/"), function(x) paste(x[1:5], collapse = "/"))

p1 = ggplot(go_table, aes(x = -log10(qvalue), y = rev(Description), fill = Cluster)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(x = 0.1, y = rev(Description), label = Description), size = 3.5, hjust = 0) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = 'black', linewidth = 0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12),
        legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#D51F26", "#00A08A", "#F2AD00", "#F98400")) +
  geom_text(data = go_table,
            aes(x = 0.1, y = rev(Description), label = geneID, color = Cluster),
            size = 4,
            fontface = 'italic', 
            hjust = 0,
            vjust = 2.3) +
  scale_color_manual(values = c("#D51F26", "#00A08A", "#F2AD00", "#F98400")) +
  scale_y_discrete(expand = c(0.1, 0)) +
  ggtitle("Enrichment of celltype marker genes")

pathway_count <- length(levels(go_table$Description))

cluster_counts <- go_table %>%
  group_by(Cluster) %>%
  summarize(count = n())

df_anno <- data.frame(Cluster = unique(go_table$Cluster))
df_anno$Cluster <- factor(df_anno$Cluster, levels = levels(go_table$Cluster))

df_anno <- merge(df_anno, cluster_counts, by = "Cluster")

df_anno <- df_anno %>%
  arrange(desc(match(Cluster, levels(go_table$Cluster)))) %>%
  mutate(
    position = cumsum(c(0, head(count, -1))) + count/2,
    ymin = cumsum(c(0, head(count, -1))),
    ymax = cumsum(count)
  )
p2 = ggplot() +
  geom_rect(data = df_anno, 
            aes(xmin = 0.5, xmax = 1.5, 
                ymin = ymin + 0.1, 
                ymax = ymax - 0.1, 
                fill = Cluster), 
            alpha = 0.3) +
  geom_text(data = df_anno,
            aes(x = 1, y = position, label = Cluster),
            size = 4,
            #fontface = 'italic', 
            angle = 90) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#D51F26", "#00A08A", "#F2AD00")) +
  scale_y_continuous(limits = c(0, pathway_count), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))
combined_plot <- p2 + p1 + plot_layout(widths = c(0.05, 1))

pdf(file = "path/PVS/image/MyoFB/MyoFB_subtype_GO_2.pdf", width = 7, height = 5.5)
print(combined_plot)
dev.off()

#SCENIC

##Prepare the cell meta information

cellInfo <- data.frame(MyoFBcell@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")] <- "sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")] <- "cluster"
colnames(cellInfo)[which(colnames(cellInfo)=="subcelltype")] <-"subcelltype"
cellInfo <- cellInfo[,c("sample","cluster","subcelltype")]
saveRDS(cellInfo, file="path/PVS/result/SCENIC/int/cellInfo.Rds")

##Prepare the expression matrix
exprMat <- as.matrix(MyoFBcell@assays$RNA$counts)
##Set up the analysis environment
mydbDIR <- "path/PVS/result/SCENIC/cisTarget"
mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather", "hg19-tss-centered-10kb-7species.mc9nr.feather")
names(mydbs) <- c("500bp","10kb")
scenicOptions <- initializeScenic(org="hgnc",nCores=100,dbDir=mydbDIR,datasetTitle ="MyoFBcell")
saveRDS(scenicOptions, "path/PVS/result/SCENIC/int/scenicOptions.rds")

##==Inference of transcriptional regulatory networks==##
##Gene filtering
genesKept <- geneFiltering(exprMat, scenicOptions,minCountsPerGene=3 * 0.01 * ncol(exprMat),minSamples=ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
##Calculate the correlation matrix
runCorrelation(exprMat_filtered, scenicOptions)

##Correlation regression analysis of TF-Targets
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)

##Infer co-expression modules
runSCENIC_1_coexNetwork2modules(scenicOptions)

data(list="motifAnnotations_hgnc_v9", package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget"))

##==Activity scoring and visualization of regulon==##
##regulons calculate the AUC value and conduct downstream analysis
exprMat_all<- as.matrix(MyoFBcell@assays$RNA$counts)
exprMat_all <- log2(exprMat_all+1)
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat_all)

runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat_all)

setwd("path/PVS/result/SCENIC")

##import regulonAUC matrix
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(MyoFBcell, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(MyoFBcell, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

library(pheatmap)
cellInfo <- readRDS("path/PVS/result/SCENIC/int/cellInfo.Rds")
celltype = subset(cellInfo,select='subcelltype')
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)

celltypecolor <- c("#D51F26", "#00A08A", "#F2AD00")
names(celltypecolor) <- c("MyoFB1", "MyoFB2", "MyoFB3")
ann_colors <- list(subcelltype = celltypecolor)
pdf(file = "path/PVS/image/MyoFB/SCENIC/MyoFBcell_regulon_AUC.pdf", width = 5, height = 4)
pheatmap(AUCmatrix, show_colnames=F, annotation_col=celltype, annotation_colors = ann_colors)
dev.off()

pdf(file = "path/PVS/image/MyoFB/SCENIC/MyoFBcell_regulon_Bin.pdf", width = 5, height = 4)
pheatmap(BINmatrix, show_colnames=F, annotation_col=celltype, annotation_colors = ann_colors,color = colorRampPalette(colors = c("white","black"))(100))
dev.off()

pdf(file = "path/PVS/image/MyoFB/SCENIC/MyoFBcell_regulon_vlnplot.pdf", width = 7, height = 4)
VlnPlot(MyoFBcell,features=c("MEIS2", "FOXO3", "GLI2", "EGR1","JUN", "FOS", "FOSB", "JUNB", "ELF2","ERG","FLI1","PPARG",
                             "MECOM", "ETS1","FOXP2","WT1","ZMAT4","PAX5","POU2F2","RARB","THRB","NR1H4","MTA3","NFE2L1",
                             "EZH2","MAZ","CREB3L2","FOXP1","SOX5","TCF12"),
        group.by="subcelltype",
        stack=T,cols=mycolour)+NoLegend()
dev.off()

#Pseudo-time analysis

## Cytotrace2

cytotrace2_result_sce <- cytotrace2(MyoFBcell, 
                                is_seurat = TRUE, 
                                slot_type = "counts", 
                                species = 'human',
                                seed = 1234)

data <- as(as.matrix(MyoFBcell@assays$RNA$counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = MyoFBcell@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

diff.genes <- read.csv('path/PVS/result/MyoFBcell_diff_genes_wilcox_subcelltype.csv')
diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
plot_ordering_genes(mycds)

#Dimensionality reduction
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')

#Sorting
mycds <- orderCells(mycds)
mycds$group <- factor(mycds$group,level = c ("H", "TAPVC", "PVS"))
mycds$orig.ident <- factor(mycds$orig.ident,level = c ("H1", "H2", "H3", "TAPVC1", "TAPVC2", "PVS1", "PVS2"))

pdf(file = "path/PVS/image/MyoFB/MyoFBcell_subcelltype_trajectory.pdf", width = 5, height = 4)
plot_cell_trajectory(mycds, color_by = "subcelltype", cell_size = 0.1) + scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/PVS/image/MyoFB/MyoFBcell_group_trajectory.pdf", width = 5, height = 4)
plot_cell_trajectory(mycds, color_by = "group", cell_size = 0.1)+ scale_colour_manual(values = c('#2967a0', '#ecc342', '#af2337'))+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/PVS/image/MyoFB/MyoFBcell_sample_trajectory.pdf", width = 5, height = 4)
plot_cell_trajectory(mycds, color_by = "orig.ident", cell_size = 0.1) + scale_colour_manual(values = c("#015c92", "#3CA4E5", "#53A7D8", "#ffd842", "#f4b528", "#ed1c24", "#b91e45")) + theme_bw() + theme(panel.grid=element_blank())
dev.off()

pdf(file = "path/PVS/image/MyoFB/MyoFBcell_state_trajectory.pdf", width = 5, height = 4)
plot_cell_trajectory(mycds, color_by = "State",cell_size = 0.1)+ scale_colour_manual(values = mycolour2)+ theme_bw() + theme(panel.grid=element_blank())
dev.off()

cols<-brewer.pal(10,"Spectral")
pal<-colorRampPalette(cols)
addcol<-pal(50)
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_cytotrace2_trajectory.pdf", width = 5, height = 3)
plot_cell_trajectory(mycds, color_by = "CytoTRACE2_Score",cell_size = 0.1)+ theme_bw() + theme(panel.grid=element_blank())+
  scale_colour_gradientn(colours = addcol)
dev.off()

pdf(file = "path/PVS/image/MyoFB/MyoFBcell_Pseudotime_ACTA2.pdf", width = 4.5, height = 3)
plot_genes_in_pseudotime(mycds["ACTA2",], color_by = "subcelltype")+ scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_Pseudotime_TGFB1.pdf", width = 4.5, height = 3)
plot_genes_in_pseudotime(mycds["TGFB1",], color_by = "subcelltype")+ scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_Pseudotime_CCN2.pdf", width = 4.5, height = 3)
plot_genes_in_pseudotime(mycds["CCN2",], color_by = "subcelltype")+ scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_Pseudotime_COL1A1.pdf", width = 4.5, height = 3)
plot_genes_in_pseudotime(mycds["COL1A1",], color_by = "subcelltype")+ scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_Pseudotime_COL3A1.pdf", width = 4.5, height = 3)
plot_genes_in_pseudotime(mycds["COL3A1",], color_by = "subcelltype")+ scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))
dev.off()
pdf(file = "path/PVS/image/MyoFB/MyoFBcell_Pseudotime_COL4A1.pdf", width = 4.5, height = 3)
plot_genes_in_pseudotime(mycds["COL4A1",], color_by = "subcelltype")+ scale_colour_manual(values = c("#D51F26", "#00A08A", "#F2AD00"))
dev.off()

EMT_genes <- read.csv("path/PVS/result/EMT_genes.txt",header = F)
gene <- as.list(EMT_genes)
MyoFBcell <- AddModuleScore(MyoFBcell,
                            features = gene,
                            ctrl = 100,
                            name = "EMT_score")
head(MyoFBcell@meta.data)
library(ggpubr)
pdf(file = "path/PVS/image/MyoFB_VlnPlot_EMT_score.pdf", width = 6, height = 4)
VlnPlot(MyoFBcell, features = "EMT_score1", ncol = 1, pt.size = 0, cols = c("#D51F26", "#00A08A", "#F2AD00"), group.by = "subcelltype")+
  geom_boxplot(width = 0.2, col = "black", fill = "white")+ theme_bw() + theme(panel.grid=element_blank())+
  stat_compare_means(comparisons = my_comparisons,method = "t.test", label = "p.signif")+ 
  theme_bw() + theme(panel.grid=element_blank()) + ylim(0,1)
dev.off()

resident_genes <- read.csv("path/PVS/result/resident_genes.txt",header = F)
gene <- as.list(resident_genes)
MyoFBcell <- AddModuleScore(MyoFBcell,
                            features = gene,
                            ctrl = 100,
                            name = "resident_score")
head(MyoFBcell@meta.data)
my_comparisons <- list( c("MyoFB1", "MyoFB2"), c("MyoFB1", "MyoFB3"), c("MyoFB2", "MyoFB3"))

pdf(file = "path/PVS/image/MyoFB_VlnPlot_resident_score.pdf", width = 6, height = 4)
VlnPlot(MyoFBcell, features = "resident_score1", ncol = 1, pt.size = 0, cols = c("#D51F26", "#00A08A", "#F2AD00"), group.by = "subcelltype")+
  geom_boxplot(width = 0.2, col = "black", fill = "white")+
  stat_compare_means(comparisons = my_comparisons,method = "t.test", label = "p.signif")+ 
  theme_bw() + theme(panel.grid=element_blank()) + ylim(0,2)
dev.off()

##Correlation analysis

Idents(MyoFBcell)<- MyoFBcell$subcelltype
expr <- GetAssayData(object = MyoFBcell, assay = "RNA", layer = "data")

av.exp<- AverageExpression(MyoFBcell, group.by = "subcelltype")
# av.exp<- av.exp[which(row.names(av.exp)%in% features),]
df <- as.data.frame(av.exp[["RNA"]])
tdc <- cor (df, method="pearson")

library("RColorBrewer")
cols<-brewer.pal(9,"YlGnBu")
pal<-colorRampPalette(cols)
addcol<-pal(50)

pdf(file = "path/PVS/image/MyoFB_subtype_corr.pdf", width = 4, height = 4)
corrplot(tdc, method = "ellipse", type = "upper", is.corr = FALSE,
         tl.col = "black", tl.cex = 0.8,tl.pos = "lt",col = addcol, col.lim = c(0.5, 1))
corrplot(tdc, method = "number", type = "lower", is.corr = FALSE,
         tl.col = "n", tl.cex = 0.8, tl.pos = "n",
         add = T,col = addcol, col.lim = c(0.3, 1))
dev.off()

#Co-expression network analysis

MyoFBcell <- SetupForWGCNA(
  MyoFBcell,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
MyoFBcell <- MetacellsByGroups(
  seurat_obj = MyoFBcell,
  group.by = c("subcelltype", "orig.ident"), # specify the columns in MyoFBcell@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'subcelltype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
MyoFBcell <- NormalizeMetacells(MyoFBcell)

###### MyoFB
MyoFBcell <- SetDatExpr(
  MyoFBcell,
  group_name = c("MyoFB3"), # the name of the group of interest in the group.by column
  group.by='subcelltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Test different soft powers:
MyoFBcell <- TestSoftPowers(
  MyoFBcell,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(MyoFBcell)

# assemble with patchwork
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_wrap_plots.pdf", width = 7, height = 5)
wrap_plots(plot_list, ncol=2)
dev.off()
power_table <- GetPowerTable(MyoFBcell)
head(power_table)

# construct co-expression network:
MyoFBcell <- ConstructNetwork(
  MyoFBcell, soft_power=5,
  setDatExpr=FALSE,
  tom_name = 'MyoFB3', overwrite_tom = TRUE # name of the topoligical overlap matrix written to disk
)

pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_Dendrogram.pdf", width = 5, height = 3)
PlotDendrogram(MyoFBcell, main='MyoFB3 hdWGCNA Dendrogram')
dev.off()

# need to run ScaleData first or else harmony throws an error:
MyoFBcell <- ScaleData(MyoFBcell, features=VariableFeatures(MyoFBcell))

# compute all MEs in the full single-cell dataset
MyoFBcell <- ModuleEigengenes(
  MyoFBcell,
  group.by.vars="orig.ident"
) 

hMEs <- GetMEs(MyoFBcell)
head(hMEs)

# module eigengenes:
MEs <- GetMEs(MyoFBcell, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
MyoFBcell <- ModuleConnectivity(
  MyoFBcell,
  group.by = 'subcelltype', group_name = 'MyoFB3'
)

# plot genes ranked by kME for each module
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_PlotKMEs.pdf", width = 8, height = 3.5)
PlotKMEs(MyoFBcell, ncol=4)
dev.off()

# get the module assignment table:
modules <- GetModules(MyoFBcell)

hub_df <- GetHubGenes(MyoFBcell, n_hubs = 10)

saveRDS(MyoFBcell, file='path/PVS/result/hdWGCNA/hdWGCNA_MyoFBcell.rds')

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  MyoFBcell,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_hMEs.pdf", width = 7, height = 3)
wrap_plots(plot_list, ncol=4)
dev.off()
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_ModuleCorrelogram.pdf", width = 4, height = 4)
ModuleCorrelogram(MyoFBcell)
dev.off()


# get hMEs from seurat object
MEs <- GetMEs(MyoFBcell, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
MyoFBcell@meta.data <- cbind(MyoFBcell@meta.data, MEs)

# plot with Seurat's DotPlot function
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_DotPlot_2.pdf", width = 4.2, height = 2)
DotPlot(MyoFBcell, features=c("turquoise", "blue", "yellow", "brown", "black", "green", "red"), group.by = 'subcelltype', cols = c('white','#990101')) + theme_bw() + theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# hubgene network
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_HubGeneNetworkPlot.pdf", width = 7, height = 7)
HubGeneNetworkPlot(
  MyoFBcell,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()

options(future.globals.maxSize = 80000 * 1024^2)
MyoFBcell <- RunModuleUMAP(
  MyoFBcell,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(MyoFBcell)

# plot with ggplot
pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_umap_df.pdf", width = 5, height = 5)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()

pdf(file = "path/PVS/image/MyoFB/hdWGCNA/MyoFB_hdWGCNA_ModuleUMAPPlot.pdf", width = 5, height = 5)
ModuleUMAPPlot(
  MyoFBcell,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE)
dev.off()

#### Renaming hdWGCNA modules
MyoFBcell <- ResetModuleNames(
  MyoFBcell,
  new_name = "M" # the base name for the new modules
)

# print out the new module names
modules <- GetModules(MyoFBcell)
print(levels(modules$module))

##### Change the colors of all modules

# get the module table
modules <- GetModules(MyoFBcell)
mods <- unique(modules$module)

# make a table of the module-color pairings
mod_colors_df <- dplyr::select(modules, c(module, color)) %>%
  distinct %>% arrange(module)
rownames(mod_colors_df) <- mod_colors_df$module

# print the dataframe
mod_colors_df

library(MetBrewer)

# get a table of just the module and it's unique color
mod_color_df <- GetModules(MyoFBcell) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)

# the number of unique modules (subtract 1 because the grey module stays grey):
n_mods <- nrow(mod_color_df) - 1

# using the "Signac" palette from metbrewer, selecting for the number of modules
new_colors <- paste0(met.brewer("Signac", n=n_mods))

# reset the module colors
MyoFBcell <- ResetModuleColors(MyoFBcell, new_colors)

############ GO analysis of Module

module_colors <- setdiff(unique(modules$module), "grey")

modules_sub <- subset(modules, module %in% c("M1", "M2", "M3", "M4", "M5", "M6", "M7"))

group <- data.frame(gene=modules_sub$gene_name,
                    group=modules_sub$module)

Gene_ID <- bitr(modules_sub$gene_name, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  #qvalueCutoff = 0.05
)

data_GO_sim <- clusterProfiler::simplify(data_GO, 
                                         cutoff=0.7, 
                                         by="p.adjust", 
                                         select_fun=min)

pdf(file = "path/PVS/image/MyoFB/hdWGCNA/hdWGCNA_module_GO_2.pdf", width = 4.5, height = 8)
dotplot(data_GO_sim, color = "p.adjust", size = "Count", showCategory=5, font.size = 5)+
  geom_point(shape=21, colour='black')+ 
  scale_color_gradientn(colours = c('#990101','white')) + RotatedAxis() + theme_bw() + theme(panel.grid=element_blank())
dev.off()

#Cellchat

MyoFB3 <- subset(MyoFBcell, subcelltype %in% c("MyoFB3"))
scRNA@meta.data$cellchat_celltype <- scRNA@meta.data$celltype
MyoFB3_id <- rownames(MyoFB3@meta.data)

test <- scRNA@meta.data
index <- rownames(test) %in% as.character(MyoFB3_id)
df[index, "cellchat_celltype"] <- "MyoFB3"

index <- rownames(test) %in% as.character(MyoFB3_id)
if ("cellchat_celltype" %in% colnames(test)) {
  if (is.factor(test$cellchat_celltype)) {
    levels(test$cellchat_celltype) <- c(levels(test$cellchat_celltype), "MyoFB3")
  }
  test[index, "cellchat_celltype"] <- "MyoFB3"
} else {
  warning("列 'cellchat_celltype' 不存在")
}

scRNA@meta.data <- test

cellchatcells <- subset(scRNA, cellchat_celltype %in% c("VSMC", "FB", "Macrophage", "T/NK", "EC", "Mast", "MyoFB3", "OLG", "EPI"))
cellchatcells$cellchat_celltype <- factor(cellchatcells$cellchat_celltype,level = c ("VSMC", "FB", "Macrophage", "T/NK", "EC", "Mast", "MyoFB3", "OLG", "EPI"))

cellchat <- createCellChat(object=cellchatcells, meta = cellchatcells@meta.data, group.by = "cellchat_celltype")
CellChatDB <- CellChatDB.human
str(CellChatDB) 

colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)

showDatabaseCategory(CellChatDB)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB

cellchat <- smoothData(cellchat, adj =PPI.human) 


#### Infer the ligand-receptor level cellular communication network

cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "path/PVS/result/cellchat/net_lr.csv")

#### Infer the cellular communication network at the signal pathway level
cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "path/PVS/result/cellchat/net_pathway.csv")


cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
cell_colors <- c("VSMC" = '#af2337', 
                 "FB" = '#ecc342', 
                 "Macrophage" = '#2967a0', 
                 "T/NK" = '#2f3c28', 
                 "EC" = '#96b437', 
                 "Mast" = '#da93ab', 
                 "MyoFB3" = '#e58932', 
                 "OLG" = '#80598f', 
                 "EPI" = '#7e331f')
pdf(file = "path/PVS/image/cellchat/PVS_all_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions", color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_all_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors)
dev.off()

##### Check the signals emitted by each type of cell

mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
pdf(file = "path/PVS/image/cellchat/PVS_single_Number_of_interactions.pdf")
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
pdf(file = "path/PVS/image/cellchat/PVS_single_Interaction_weights_strength.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()

#### Visualization of multiple ligand-receptor-mediated cellular interactions
levels(cellchat@idents)
# show all the significant interactions (L-R pairs)
####### Extract the top 20 important ligand-receptor interactions

interaction_weights <- cellchat@net$prob
interaction_strength <- apply(interaction_weights, 3, sum)
top20_interactions <- names(sort(interaction_strength, decreasing = TRUE))[1:20]

pdf(file = "path/PVS/image/cellchat/PVS_all_plot_L-R_pairs_top20.pdf",width = 4, height = 5)
pairLR.use <- data.frame(
  interaction_name = top20_interactions)
netVisual_bubble(cellchat, sources.use = c(7,5), 
                 targets.use = c(1,2,5,7), remove.isolate = FALSE, pairLR.use = pairLR.use)
dev.off()

####  Analysis of Communication Network systems
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP",thresh = 0.01)
pdf(file = "path/PVS/image/cellchat/PVS_all_outgoing.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_all_incoming.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_all_incoming.pdf",width = 4, height = 3)
netAnalysis_signalingRole_scatter(cellchat, color.use = cell_colors)
dev.off()

###### Group Level

cellchatcells_H <- subset(cellchatcells, group %in% "H")
cellchatcells_TAPVC <- subset(cellchatcells, group %in% "TAPVC")
cellchatcells_PVS <- subset(cellchatcells, group %in% "PVS")

cellchat_H <- createCellChat(object=cellchatcells_H, meta = cellchatcells_H@meta.data, group.by = "cellchat_celltype")
cellchat_TAPVC <- createCellChat(object=cellchatcells_TAPVC, meta = cellchatcells_TAPVC@meta.data, group.by = "cellchat_celltype")
cellchat_PVS <- createCellChat(object=cellchatcells_PVS, meta = cellchatcells_PVS@meta.data, group.by = "cellchat_celltype")

cellchat_H@DB <- CellChatDB
cellchat_TAPVC@DB <- CellChatDB
cellchat_PVS@DB <- CellChatDB

cellchat_H <- subsetData(cellchat_H)
cellchat_TAPVC <- subsetData(cellchat_TAPVC)
cellchat_PVS <- subsetData(cellchat_PVS)

cellchat_H <- identifyOverExpressedGenes(cellchat_H)
cellchat_H <- identifyOverExpressedInteractions(cellchat_H)
cellchat_H <- smoothData(cellchat_H, adj =PPI.human) 

cellchat_TAPVC <- identifyOverExpressedGenes(cellchat_TAPVC)
cellchat_TAPVC <- identifyOverExpressedInteractions(cellchat_TAPVC)
cellchat_TAPVC <- smoothData(cellchat_TAPVC, adj =PPI.human) 

cellchat_PVS <- identifyOverExpressedGenes(cellchat_PVS)
cellchat_PVS <- identifyOverExpressedInteractions(cellchat_PVS)
cellchat_PVS <- smoothData(cellchat_PVS, adj =PPI.human) 

cellchat_H <- computeCommunProb(cellchat_H, raw.use = FALSE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_H <- filterCommunication(cellchat_H, min.cells = 10)
df.net <- subsetCommunication(cellchat_H)
write.csv(df.net, "path/PVS/result/cellchat/net_lr_H.csv")
cellchat_H <- computeCommunProbPathway(cellchat_H)
df.netp <- subsetCommunication(cellchat_H, slot.name = "netP")
write.csv(df.netp, "path/PVS/result/cellchat/net_pathway_H.csv")
cellchat_H <- aggregateNet(cellchat_H)
groupSize <- as.numeric(table(cellchat_H@idents))
par(mfrow = c(1,2), xpd=TRUE)

cellchat_TAPVC <- computeCommunProb(cellchat_TAPVC, raw.use = FALSE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_TAPVC <- filterCommunication(cellchat_TAPVC, min.cells = 10)
df.net <- subsetCommunication(cellchat_TAPVC)
write.csv(df.net, "path/PVS/result/cellchat/net_lr_H.csv")
cellchat_TAPVC <- computeCommunProbPathway(cellchat_TAPVC)
df.netp <- subsetCommunication(cellchat_TAPVC, slot.name = "netP")
write.csv(df.netp, "path/PVS/result/cellchat/net_pathway_H.csv")
cellchat_TAPVC <- aggregateNet(cellchat_TAPVC)
groupSize <- as.numeric(table(cellchat_TAPVC@idents))
par(mfrow = c(1,2), xpd=TRUE)

cellchat_PVS <- computeCommunProb(cellchat_PVS, raw.use = FALSE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_PVS <- filterCommunication(cellchat_PVS, min.cells = 10)
df.net <- subsetCommunication(cellchat_PVS)
write.csv(df.net, "path/PVS/result/cellchat/net_lr_H.csv")
cellchat_PVS <- computeCommunProbPathway(cellchat_PVS)
df.netp <- subsetCommunication(cellchat_PVS, slot.name = "netP")
write.csv(df.netp, "path/PVS/result/cellchat/net_pathway_H.csv")
cellchat_PVS <- aggregateNet(cellchat_PVS)
groupSize <- as.numeric(table(cellchat_PVS@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf(file = "path/PVS/image/cellchat/PVS_H_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_circle(cellchat_H@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions", color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_H_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_circle(cellchat_H@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_circle(cellchat_TAPVC@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions", color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_circle(cellchat_TAPVC@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_PVS_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_circle(cellchat_PVS@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions", color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_PVS_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_circle(cellchat_PVS@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength", color.use = cell_colors)
dev.off()

##### Check the signals emitted by each type of cell

mat <- cellchat_H@net$count
par(mfrow = c(3,3), xpd=TRUE)
pdf(file = "path/PVS/image/cellchat/PVS_H_single_Number_of_interactions.pdf")
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()
mat <- cellchat_H@net$weight
par(mfrow = c(3,3), xpd=T)
pdf(file = "path/PVS/image/cellchat/PVS_H_single_Interaction_weights_strength.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()

mat <- cellchat_TAPVC@net$count
par(mfrow = c(3,3), xpd=TRUE)
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_single_Number_of_interactions.pdf")
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()
mat <- cellchat_TAPVC@net$weight
par(mfrow = c(3,3), xpd=T)
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_single_Interaction_weights_strength.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()

mat <- cellchat_PVS@net$count
par(mfrow = c(3,3), xpd=TRUE)
pdf(file = "path/PVS/image/cellchat/PVS_PVS_single_Number_of_interactions.pdf")
for (i in 1:nrow(mat)) {
  # i = 1
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()
mat <- cellchat_PVS@net$weight
par(mfrow = c(3,3), xpd=T)
pdf(file = "path/PVS/image/cellchat/PVS_PVS_single_Interaction_weights_strength.pdf")
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use = cell_colors)
}
dev.off()

#### Visualization of multiple ligand-receptor-mediated cellular interactions
levels(cellchat_H@idents)
# show all the significant interactions (L-R pairs)
####### Extract the top 20 important ligand-receptor interactions

interaction_weights <- cellchat_H@net$prob
interaction_strength <- apply(interaction_weights, 3, sum)
top20_interactions <- names(sort(interaction_strength, decreasing = TRUE))[1:20]
pairLR.use <- data.frame(interaction_name=top20_interactions)

#Receptor cells and ligand cells need to be specified
pdf(file = "path/PVS/image/cellchat/PVS_H_plot_L-R_pairs_top20.pdf",width = 4, height = 5)
pairLR.use <- data.frame(
  interaction_name = top20_interactions)
netVisual_bubble(cellchat, sources.use = c(7,5), 
                 targets.use = c(1,2,5,7), remove.isolate = FALSE, pairLR.use = pairLR.use)
dev.off()

#### Visualization of multiple ligand-receptor-mediated cellular interactions
levels(cellchat_TAPVC@idents)
# show all the significant interactions (L-R pairs)
interaction_weights <- cellchat_TAPVC@net$prob
interaction_strength <- apply(interaction_weights, 3, sum)
top20_interactions <- names(sort(interaction_strength, decreasing = TRUE))[1:20]
pairLR.use <- data.frame(interaction_name=top20_interactions)
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_plot_L-R_pairs_top20.pdf",width = 4, height = 5)
pairLR.use <- data.frame(
  interaction_name = top20_interactions)
netVisual_bubble(cellchat, sources.use = c(7,5), 
                 targets.use = c(1,2,5,7), remove.isolate = FALSE, pairLR.use = pairLR.use)
dev.off()

levels(cellchat_PVS@idents)
# show all the significant interactions (L-R pairs)
interaction_weights <- cellchat_PVS@net$prob
interaction_strength <- apply(interaction_weights, 3, sum)
top20_interactions <- names(sort(interaction_strength, decreasing = TRUE))[1:20]
pairLR.use <- data.frame(interaction_name=top20_interactions)

pdf(file = "path/PVS/image/cellchat/PVS_PVS_plot_L-R_pairs_top20.pdf",width = 4, height = 5)
pairLR.use <- data.frame(
  interaction_name = top20_interactions)
netVisual_bubble(cellchat, sources.use = c(7,5), 
                 targets.use = c(1,2,5,7), remove.isolate = FALSE, pairLR.use = pairLR.use)
dev.off()

####Analysis of Communication Network systems
cellchat_H <- netAnalysis_computeCentrality(cellchat_H, slot.name = "netP",thresh = 0.01)
pdf(file = "path/PVS/image/cellchat/PVS_H_outgoing.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_H, pattern = "outgoing", font.size = 5, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_H_incoming.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_H, pattern = "incoming", font.size = 5, color.use = cell_colors)
dev.off()
cellchat_TAPVC <- netAnalysis_computeCentrality(cellchat_TAPVC, slot.name = "netP",thresh = 0.01)
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_outgoing.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_TAPVC, pattern = "outgoing", font.size = 5, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_incoming.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_TAPVC, pattern = "incoming", font.size = 5, color.use = cell_colors)
dev.off()
cellchat_PVS <- netAnalysis_computeCentrality(cellchat_PVS, slot.name = "netP",thresh = 0.01)
pdf(file = "path/PVS/image/cellchat/PVS_PVS_outgoing.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_PVS, pattern = "outgoing", font.size = 5, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/PVS_PVS_incoming.pdf",width = 6, height = 7)
netAnalysis_signalingRole_heatmap(cellchat_PVS, pattern = "incoming", font.size = 5, color.use = cell_colors)
dev.off()

cellchat_H <- netAnalysis_computeCentrality(cellchat_H, slot.name = "netP")
pdf(file = "path/PVS/image/cellchat/PVS_H_sources_targets.pdf",width = 4, height = 3)
netAnalysis_signalingRole_scatter(cellchat_H)
dev.off()
cellchat_TAPVC <- netAnalysis_computeCentrality(cellchat_TAPVC, slot.name = "netP")
pdf(file = "path/PVS/image/cellchat/PVS_TAPVC_sources_targets.pdf",width = 4, height = 3)
netAnalysis_signalingRole_scatter(cellchat_TAPVC)
dev.off()
cellchat_PVS <- netAnalysis_computeCentrality(cellchat_PVS, slot.name = "netP")
pdf(file = "path/PVS/image/cellchat/PVS_PVS_sources_targets.pdf",width = 4, height = 3)
netAnalysis_signalingRole_scatter(cellchat_PVS)
dev.off()


###### H vs TAPVC

cellchat_H_TAPVC.list <- list(H=cellchat_H, TAPVC=cellchat_TAPVC)
cellchat_H_TAPVC <- mergeCellChat(cellchat_H_TAPVC.list, add.names = names(cellchat_H_TAPVC.list), cell.prefix = TRUE)

pdf(file = "path/PVS/image/cellchat/HvsTAPVC_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_diffInteraction(cellchat_H_TAPVC, weight.scale = T, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/HvsTAPVC_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_diffInteraction(cellchat_H_TAPVC, weight.scale = T, measure = "weight", color.use = cell_colors)
dev.off()

pathway.union <- union(cellchat_H_TAPVC.list[[1]]@netP$pathways, cellchat_H_TAPVC.list[[2]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_H_TAPVC.list[[1]]))
ht2 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC.list[[2]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_H_TAPVC.list[[2]]))
pdf(file = "path/PVS/image/cellchat/HvsTAPVC_netAnalysis_signalingRole_heatmap.pdf", width = 12, height = 5)
ht1+ht2
dev.off()

###### H vs PVS
cellchat_H_PVS.list <- list(H=cellchat_H, PVS=cellchat_PVS)
cellchat_H_PVS <- mergeCellChat(cellchat_H_PVS.list, add.names = names(cellchat_H_PVS.list), cell.prefix = TRUE)

pdf(file = "path/PVS/image/cellchat/HvsPVS_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_diffInteraction(cellchat_H_PVS, weight.scale = T, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/HvsPVS_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_diffInteraction(cellchat_H_PVS, weight.scale = T, measure = "weight", color.use = cell_colors)
dev.off()

pathway.union <- union(cellchat_H_PVS.list[[1]]@netP$pathways, cellchat_H_PVS.list[[2]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(cellchat_H_PVS.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_H_PVS.list[[1]]))
ht2 = netAnalysis_signalingRole_heatmap(cellchat_H_PVS.list[[2]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_H_PVS.list[[2]]))
pdf(file = "path/PVS/image/cellchat/HvsPVS_netAnalysis_signalingRole_heatmap.pdf", width = 12, height = 5)
ht1+ht2
dev.off()

###### TAPVC vs PVS
cellchat_TAPVC_PVS.list <- list(H=cellchat_TAPVC, PVS=cellchat_PVS)
cellchat_TAPVC_PVS <- mergeCellChat(cellchat_TAPVC_PVS.list, add.names = names(cellchat_TAPVC_PVS.list), cell.prefix = TRUE)

pdf(file = "path/PVS/image/cellchat/TAPVCvsPVS_Number_of_interactions.pdf", width = 5, height = 5)
netVisual_diffInteraction(cellchat_TAPVC_PVS, weight.scale = T, color.use = cell_colors)
dev.off()
pdf(file = "path/PVS/image/cellchat/TAPVCvsPVS_Interaction_weights_strength.pdf", width = 5, height = 5)
netVisual_diffInteraction(cellchat_TAPVC_PVS, weight.scale = T, measure = "weight", color.use = cell_colors)
dev.off()

pathway.union <- union(cellchat_TAPVC_PVS.list[[1]]@netP$pathways, cellchat_TAPVC_PVS.list[[2]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(cellchat_TAPVC_PVS.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_TAPVC_PVS.list[[1]]))
ht2 = netAnalysis_signalingRole_heatmap(cellchat_TAPVC_PVS.list[[2]], pattern = "all", signaling = pathway.union, 
                                        title = names(cellchat_TAPVC_PVS.list[[2]]))
pdf(file = "path/PVS/image/cellchat/TAPVCvsPVS_netAnalysis_signalingRole_heatmap.pdf", width = 12, height = 5)
ht1+ht2
dev.off()

cellchat_H_TAPVC_PVS.list <- list(H=cellchat_H, TAPVC=cellchat_TAPVC, PVS=cellchat_PVS)
cellchat_H_TAPVC_PVS <- mergeCellChat(cellchat_H_TAPVC_PVS.list, add.names = names(cellchat_H_TAPVC_PVS.list), cell.prefix = TRUE)

pathway.all <- intersect(c(cellchat_H_TAPVC_PVS.list[[1]]@netP$pathways, cellchat_H_TAPVC_PVS.list[[2]]@netP$pathways),cellchat_H_TAPVC_PVS.list[[3]]@netP$pathways) 

ht1 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC_PVS.list[[1]], pattern = "all", signaling = pathway.all, 
                                        title = names(cellchat_H_TAPVC_PVS.list[[1]]), color.use = cell_colors)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC_PVS.list[[2]], pattern = "all", signaling = pathway.all, 
                                        title = names(cellchat_H_TAPVC_PVS.list[[2]]), color.use = cell_colors)
ht3 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC_PVS.list[[3]], pattern = "all", signaling = pathway.all, 
                                        title = names(cellchat_H_TAPVC_PVS.list[[3]]), color.use = cell_colors)
pdf(file = "path/PVS/image/cellchat/HvsTAPVCvsPVS_netAnalysis_signalingRole_heatmap_all.pdf", width = 20, height = 5)
ht1+ht2+ht3
dev.off()

pathway.union <- union(c(cellchat_H_TAPVC_PVS.list[[1]]@netP$pathways, cellchat_H_TAPVC_PVS.list[[2]]@netP$pathways),cellchat_H_TAPVC_PVS.list[[3]]@netP$pathways)
pathway.diff <- pathway.union[!pathway.union %in% pathway.all] 

pathway.diff <- factor(pathway.diff,level = c ("VWF","CD34","SELPLG","NPR2","CLDN","ANGPTL","SEMA5","CD48","Glutamate","GP1BA",
                                               "FLT3","MSTN","5alphaP","PRL","IFN-I","CEACAM","SLITRK","TAC","NPNT","HH","SELL","27HC"))
pathway.diff_2 <- c("VWF","CD34","SELPLG","NPR2","CLDN","ANGPTL","SEMA5","CD48","Glutamate","GP1BA",
                    "FLT3","MSTN","5alphaP","PRL","IFN-I","CEACAM","SLITRK","TAC","NPNT","HH","SELL","27HC")

ht1 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC_PVS.list[[1]], pattern = "all", signaling = pathway.diff_2, 
                                        title = names(cellchat_H_TAPVC_PVS.list[[1]]), color.use = cell_colors)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC_PVS.list[[2]], pattern = "all", signaling = pathway.diff_2, 
                                        title = names(cellchat_H_TAPVC_PVS.list[[2]]), color.use = cell_colors)
ht3 = netAnalysis_signalingRole_heatmap(cellchat_H_TAPVC_PVS.list[[3]], pattern = "all", signaling = pathway.diff_2, 
                                        title = names(cellchat_H_TAPVC_PVS.list[[3]]), color.use = cell_colors)
pdf(file = "path/PVS/image/cellchat/HvsTAPVCvsPVS_netAnalysis_signalingRole_heatmap_unique.pdf", width = 20, height = 5)
ht1+ht2+ht3
dev.off()

#iTALK

######## All group #######
cellchatcells <- readRDS("path/PVS/result/cellchat/cellchatcells.rds")
iTalk_data <- as.data.frame(t(cellchatcells@assays$RNA$counts+1))
iTalk_data$cell_type <- cellchatcells@meta.data$cellchat_celltype
iTalk_data$compare_group <- cellchatcells@meta.data$group		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/PVS/result/iTALK/highly_exprs_genes.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/PVS/result/iTALK/ALL_LRpairs_Overview.csv")

pdf(file = "path/PVS/image/iTALK/All_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1.2,edge.label.cex = 0.9,vertex.size=30,arrow.width=3,edge.max.width=10,margin=0.2)
dev.off()
pdf(file = "path/PVS/image/iTALK/All_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

######## H vs TAPVC ########
data2 <- subset(iTalk_data, subset=iTalk_data$compare_group=="H"|iTalk_data$compare_group=="TAPVC")
deg_VSMC<-DEG(data2 %>% filter(cell_type=="VSMC"),method='DESeq2',contrast=c("H","TAPVC"))
deg_FB<-DEG(data2 %>% filter(cell_type=='FB'),method='DESeq2',contrast=c("H","TAPVC"))
deg_Macrophage<-DEG(data2 %>% filter(cell_type=='Macrophage'),method='DESeq2',contrast=c("H","TAPVC"))
deg_T_NK<-DEG(data2 %>% filter(cell_type=="T/NK"),method='DESeq2',contrast=c("H","TAPVC"))
deg_EC<-DEG(data2 %>% filter(cell_type=="EC"),method='DESeq2',contrast=c("H","TAPVC"))
deg_Mast<-DEG(data2 %>% filter(cell_type=="Mast"),method='DESeq2',contrast=c("H","TAPVC"))
deg_MyoFB<-DEG(data2 %>% filter(cell_type=="MyoFB3"),method='DESeq2',contrast=c("H","TAPVC"))
deg_OLG<-DEG(data2 %>% filter(cell_type=="OLG"),method='DESeq2',contrast=c("H","TAPVC"))
deg_EPI<-DEG(data2 %>% filter(cell_type=="EPI"),method='DESeq2',contrast=c("H","TAPVC"))
deg_all<-rbind(deg_VSMC,deg_FB,deg_Macrophage,deg_T_NK,deg_EC,deg_Mast,deg_MyoFB,deg_OLG,deg_EPI)
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(deg_all, datatype='DEG',comm_type=comm_type)
  res<-rbind(res,res_cat)
}
write.csv(res, "path/PVS/result/iTALK/H_vs_TAPVC_res.csv")
res_2<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:50,]
write.csv(res_2, "path/PVS/result/iTALK/H_vs_TAPVC_res_top50.csv", quote = FALSE)
dat <- read.csv("path/PVS/result/iTALK/H_vs_TAPVC_res_top50.csv")
pdf(file = "path/PVS/image/iTALK/H_vs_TAPVC_iTALK_LRPlot.pdf")
LRPlot(dat,datatype='DEG',link.arr.lwd=dat$cell_from_logFC,link.arr.width=dat$cell_to_logFC)
dev.off()

######## H vs PVS ########
data3 <- subset(iTalk_data, subset=iTalk_data$compare_group=="H"|iTalk_data$compare_group=="PVS")
deg_VSMC<-DEG(data3 %>% filter(cell_type=="VSMC"),method='DESeq2',contrast=c("H","PVS"))
deg_FB<-DEG(data3 %>% filter(cell_type=='FB'),method='DESeq2',contrast=c("H","PVS"))
deg_Macrophage<-DEG(data3 %>% filter(cell_type=='Macrophage'),method='DESeq2',contrast=c("H","PVS"))
deg_T_NK<-DEG(data3 %>% filter(cell_type=="T/NK"),method='DESeq2',contrast=c("H","PVS"))
deg_EC<-DEG(data3 %>% filter(cell_type=="EC"),method='DESeq2',contrast=c("H","PVS"))
deg_Mast<-DEG(data3 %>% filter(cell_type=="Mast"),method='DESeq2',contrast=c("H","PVS"))
deg_MyoFB<-DEG(data3 %>% filter(cell_type=="MyoFB3"),method='DESeq2',contrast=c("H","PVS"))
deg_OLG<-DEG(data3 %>% filter(cell_type=="OLG"),method='DESeq2',contrast=c("H","PVS"))
deg_EPI<-DEG(data3 %>% filter(cell_type=="EPI"),method='DESeq2',contrast=c("H","PVS"))
deg_all<-rbind(deg_VSMC,deg_FB,deg_Macrophage,deg_T_NK,deg_EC,deg_Mast,deg_MyoFB,deg_OLG,deg_EPI)
res<-NULL
for(comm_type in comm_list){
  res_cat<-FindLR(deg_all, datatype='DEG',comm_type=comm_type)
  res<-rbind(res,res_cat)
}
res_2<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),][1:50,]
write.csv(res_2, "path/PVS/result/iTALK/H_vs_PVS_res_top50.csv", quote = FALSE)
dat <- read.csv("path/PVS/result/iTALK/H_vs_PVS_res_top50.csv")
pdf(file = "path/PVS/image/iTALK/H_vs_PVS_iTALK_LRPlot.pdf")
LRPlot(dat,datatype='DEG',link.arr.lwd=dat$cell_from_logFC,link.arr.width=dat$cell_to_logFC)
dev.off()

##################### H #####################

H <- subset(cellchatcells, group %in% c("H"))
iTalk_data <- as.data.frame(t(H@assays$RNA$counts+1))
iTalk_data$cell_type <- H@meta.data$cellchat_celltype
iTalk_data$compare_group <- H@meta.data$group		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/PVS/result/iTALK/highly_exprs_genes_H.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/PVS/result/iTALK/ALL_LRpairs_Overview_H.csv")

pdf(file = "path/PVS/image/iTALK/H_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/PVS/image/iTALK/H_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

##################### TAPVC #####################

TAPVC <- subset(cellchatcells, group %in% c("TAPVC"))
iTalk_data <- as.data.frame(t(TAPVC@assays$RNA$counts+1))
iTalk_data$cell_type <- TAPVC@meta.data$cellchat_celltype
iTalk_data$compare_group <- TAPVC@meta.data$group		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/PVS/result/iTALK/highly_exprs_genes_TAPVC.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)
iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/PVS/result/iTALK/ALL_LRpairs_Overview_TAPVC.csv")

pdf(file = "path/PVS/image/iTALK/TAPVC_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/PVS/image/iTALK/TAPVC_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

##################### PVS #####################

PVS <- subset(cellchatcells, group %in% c("PVS"))
iTalk_data <- as.data.frame(t(PVS@assays$RNA$counts+1))
iTalk_data$cell_type <- PVS@meta.data$cellchat_celltype
iTalk_data$compare_group <- PVS@meta.data$group		 
highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
saveRDS(highly_exprs_genes, 'path/PVS/result/iTALK/highly_exprs_genes_PVS.rds')

# Communication type
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(mycolour2[1:length(cell_types)], names=cell_types)

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res, "path/PVS/result/iTALK/ALL_LRpairs_Overview_PVS.csv")

pdf(file = "path/PVS/image/iTALK/PVS_iTALK.pdf")
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,edge.label.cex = 0.001,vertex.size=30,arrow.width=0.5,edge.max.width=5,margin=0.2)
dev.off()
pdf(file = "path/PVS/image/iTALK/PVS_iTALK_LRPlot.pdf")
LRPlot(iTalk_res[1:50,],datatype='mean count',link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:50],link.arr.width=iTalk_res$cell_to_mean_exprs[1:50])
dev.off()

############## SRT-seq ##############

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(AUCell)
library(ggpubr)
library(cowplot)
library(reshape2)

#### Merge
ctl <- readRDS("path/PVS/space/result/ctl_original.rds")
Assays(ctl)
DefaultAssay(ctl) <- "Spatial.008um"
ctl_qc <- subset(ctl, subset = nFeature_Spatial.008um > 2)
ctl_qc <- NormalizeData(ctl_qc)
ctl_qc <- FindVariableFeatures(ctl_qc)
ctl_qc <- ScaleData(ctl_qc)

PVS <- readRDS("path/PVS/space/result/PVS_original.rds")
Assays(PVS)
DefaultAssay(PVS) <- "Spatial.008um"
PVS_qc <- subset(PVS, subset = nFeature_Spatial.008um > 2)
PVS_qc <- NormalizeData(PVS_qc)
PVS_qc <- FindVariableFeatures(PVS_qc)
PVS_qc <- ScaleData(PVS_qc)

ctl_qc@meta.data$group <- rep("ctl", 210791)
PVS_qc@meta.data$group <- rep("TAPVC_PVS", 191986)
data.merge <- merge(ctl_qc, PVS_qc, add.cell.ids = c("ctl", "PVS"))

DefaultAssay(data.merge) <- "Spatial.008um"
VariableFeatures(data.merge) <- c(VariableFeatures(ctl), VariableFeatures(PVS))

# normalize both 8um and 16um bins
DefaultAssay(data.merge) <- "Spatial.008um"
data.merge <- NormalizeData(data.merge)

data.merge <- FindVariableFeatures(data.merge)
data.merge <- ScaleData(data.merge)

# we select 50,0000 cells and create a new 'sketch' assay
data.merge <- SketchData(
  object = data.merge,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(data.merge) <- "sketch"
saveRDS(data.merge, file = "path/PVS/space/result/data.merge_original.rds")

# perform clustering workflow
data.merge <- readRDS("path/PVS/space/result/data.merge_original.rds")

data.merge <- FindVariableFeatures(data.merge)
data.merge <- ScaleData(data.merge)
data.merge <- RunPCA(data.merge, assay = "sketch", reduction.name = "pca.sketch")
data.merge <- FindNeighbors(data.merge, assay = "sketch", reduction = "pca.sketch", dims = 1:30)
data.merge <- FindClusters(data.merge, cluster.name = "seurat_cluster.sketched", resolution = 0.1)
data.merge <- RunUMAP(data.merge, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:30)

data.merge <- ProjectData(
  object = data.merge,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:30,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(data.merge) <- "sketch"
Idents(data.merge) <- "seurat_cluster.sketched"
p1 <- DimPlot(data.merge, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(data.merge) <- "Spatial.008um"
Idents(data.merge) <- "seurat_cluster.projected"
p2 <- DimPlot(data.merge, reduction = "full.umap.sketch", cols = mycolour, label = TRUE) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2


pdf(file = "path/PVS/space/image/merge_SpatialDimPlot_r0.1.pdf", width = 10, height = 4)
SpatialDimPlot(data.merge, group.by = "seurat_cluster.projected", label = T, repel = T, label.size = 4)
dev.off()


features <- c("ACTA2", "MYH11", "LUM", "DCN", "CD68", "AIF1", "CD3D", "CD3E", "PECAM1", "VWF", "KLRD1", "NKG7", "TPSAB1", "TPSB2", "PLP1", "NRXN1", "CD79A", "MS4A1", "MZB1", "JCHAIN", "KRT19", "MSLN", "ITGBL1",
              "CD44", "NCAM1","THY1", "ENG","CD34", "PTPRC")
pdf(file = "path/PVS/space/image/merge_Plot_markergene_0.1.pdf", width = 7, height = 5)
DotPlot(data.merge, feature = features, group.by = "seurat_cluster.projected", cols = c("lightgrey", "red")) + RotatedAxis()
dev.off()

DefaultAssay(data.merge) <- "Spatial.008um"
DefaultAssay(ctl_qc) <- "Spatial.008um"
DefaultAssay(PVS_qc) <- "Spatial.008um"
pdf(file = "path/PVS/space/image/ctl_ACTA2_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(ctl_qc, features = "ACTA2") + ggtitle("ACTA2 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(file = "path/PVS/space/image/PVS_ACTA2_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(PVS_qc, features = "ACTA2") + ggtitle("ACTA2 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(file = "path/PVS/space/image/ctl_MYH11_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(ctl_qc, features = "MYH11") + ggtitle("MYH11 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()
pdf(file = "path/PVS/space/image/PVS_MYH11_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(PVS_qc, features = "MYH11") + ggtitle("MYH11 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(file = "path/PVS/space/image/ctl_LUM_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(ctl_qc, features = "LUM") + ggtitle("MYH11 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()
pdf(file = "path/PVS/space/image/PVS_LUM_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(PVS_qc, features = "LUM") + ggtitle("LUM expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(file = "path/PVS/space/image/ctl_DCN_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(ctl_qc, features = "DCN") + ggtitle("DCN expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()
pdf(file = "path/PVS/space/image/PVS_DCN_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(PVS_qc, features = "DCN") + ggtitle("DCN expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(file = "path/PVS/space/image/ctl_PECAM1_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(ctl_qc, features = "PECAM1") + ggtitle("PECAM1 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()
pdf(file = "path/PVS/space/image/PVS_PECAM1_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(PVS_qc, features = "PECAM1") + ggtitle("PECAM1 expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

pdf(file = "path/PVS/space/image/ctl_VWF_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(ctl_qc, features = "VWF") + ggtitle("VWF expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()
pdf(file = "path/PVS/space/image/PVS_VWF_expression.pdf", width = 3, height = 4)
SpatialFeaturePlot(PVS_qc, features = "VWF") + ggtitle("VWF expression (8um)") + scale_fill_gradient(low = "lightgrey", high = "red")
dev.off()

### annotation
current.cluster.ids <- c(0:11)
new.cluster.ids <- c("MyoFB", "VSMC", "EC", "MyoFB", "FB", "FB", "VSMC", "EC", "others", "others",  "others", "others")
data.merge@meta.data$celltype = plyr::mapvalues(x = data.merge@meta.data[,"seurat_cluster.projected"], from = current.cluster.ids, to = new.cluster.ids)
data.merge$celltype <- factor(data.merge$celltype,level = c ("VSMC", "FB", "EC", "MyoFB", "others"))

my_cols <- c("VSMC" = "#af2337", "FB" = "#ecc342", "EC" = "#96b437", "MyoFB" = '#e58932', "others" = 'grey')

pdf(file = "path/PVS/space/image/merge_celltype_r0.1.pdf", width = 7, height = 4)
SpatialDimPlot(data.merge, group.by = "celltype", cols = my_cols)
dev.off()

features <- c("ACTA2", "MYH11", "LUM", "DCN", "CD68", "AIF1", "CD3D", "CD3E", "KLRD1", "NKG7", "PECAM1", "VWF", "TPSAB1", "TPSB2", "PLP1", "NRXN1", "KRT19", "MSLN")
pdf(file = "path/PVS/space/image/merge_celltype_DotPlot.pdf", width = 5.5, height = 2)
DotPlot(data.merge, feature = features, group.by = "celltype", cols = c("white", '#990101')) + theme_bw() + theme(panel.grid=element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ECM AUCell

ECM_genes <- read.csv("path/PVS/space/result/ECM_genelist.txt",header = F)
gene <- as.list(ECM_genes)

data.merge <- JoinLayers(data.merge)
data.merge <- AddModuleScore(data.merge,
                            features = gene,
                            ctrl = 100,
                            name = "ECM_score")

pdf(file = "path/PVS/space/image/merge_ECM_score.pdf", width = 7, height = 4)
SpatialFeaturePlot(data.merge, features = "ECM_score1") + ggtitle("ECM_score")
dev.off()

df <- table(data.merge$celltype,data.merge$lastgroup) %>% melt()
colnames(df) <- c("Cluster","Group","Number")
df$Cluster <- factor(df$Cluster,level = c ("VSMC", "FB", "EC", "MyoFB", "others"))
df$Group<- factor(df$Group,level = c ("ctl", "TAPVC", "PVS"))

pdf(file = "path/PVS/space/image/space_ctl_TAPVC_PVS_celltype_proration.pdf", width = 5, height = 3)
ggplot(data = df, aes(x =Number, y = Group, fill = Cluster)) +
  geom_bar(stat = "identity", width=0.8, position="fill")+
  scale_fill_manual(values=my_cols) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90)
  )