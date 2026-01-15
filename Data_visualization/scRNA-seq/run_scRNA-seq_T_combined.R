### Data analysis for scRNA-seq data in combined T samples
library(Seurat)
library(ggplot2)

# define function for read data
ReadData10X <- function(path, sample){
  count <- Read10X(data.dir = path)
  data_obj <- CreateSeuratObject(counts = count, project = "STOmics", 
                                 min.cells = 3, min.features = 200)
  data_obj[["sample"]] <- sample
  return(data_obj)
}

### load T data
data_A2T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/A2T_run/outs/filtered_feature_bc_matrix/","A2T")
data_A3T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/A3T_run/outs/filtered_feature_bc_matrix/","A3T")
data_B1T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/B1T_run/outs/filtered_feature_bc_matrix/","B1T")
data_B2T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/B2T_run/outs/filtered_feature_bc_matrix/","B2T")
data_C2T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/C2T_run/outs/filtered_feature_bc_matrix/","C2T")
data_C3T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/C3T_run/outs/filtered_feature_bc_matrix/","C3T")
data_D1T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/D1T_run/outs/filtered_feature_bc_matrix/","D1T")
data_D2T <- ReadData10X("/mnt/hpc/data/2results/01_scRNA/D2T2_run/outs/filtered_feature_bc_matrix/","D2T")
data_A2T
#An object of class Seurat
# 26387 features across 12033 samples within 1 assay 
# Active assay: RNA (26387 features, 0 variable features)

data_A2T@assays$RNA@counts[1:5,1:5]
#5 x 5 sparse Matrix of class "dgCMatrix"
#             AAACCCAAGAGGTTTA-1 AAACCCAAGATCCAAA-1 AAACCCACAATACCTG-1
#RPS6KA2                       .                  5                 16
#LOC101866160                  .                  .                  .
#rna76575                      .                  .                  .

### combine all samples
data_combined <- merge(x=data_A2T,y=c(data_A3T,data_B1T,data_B2T,data_C2T,data_C3T,data_D1T,data_D2T),
                        add.cell.ids=c("A2T","A3T","B1T","B2T","C2T","C3T","D1T","D2T"))
dataName="T_combined"

# plot QC metrics
feats <- c("nFeature_RNA", "nCount_RNA")
sample_colors <- as.vector(ArchR::ArchRPalettes$paired)
VlnPlot(data_combined, group.by = "sample", features = feats, ncol = 3, pt.size = 0)
FeatureScatter(data_combined, "nCount_RNA", "nFeature_RNA", group.by = "sample", pt.size = 0.5) + scale_color_manual(values=sample_colors)
# Filter cells
data_combined <- subset(data_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)
feats <- c("nFeature_RNA", "nCount_RNA")
VlnPlot(data_combined, group.by = "sample", features = feats, ncol = 3, pt.size = 0)
FeatureScatter(data_combined, "nCount_RNA", "nFeature_RNA", group.by = "sample", pt.size = 0.5) + scale_color_manual(values=sample_colors)

### integration with harmony ###
library(harmony)

data_combined <- NormalizeData(data_combined)
data_combined = FindVariableFeatures(data_combined, verbose = F)
data_combined = ScaleData(data_combined, features = rownames(data_combined), verbose = F)

data_combined <- RunPCA(data_combined, npcs = 50, verbose = FALSE)
ElbowPlot(data_combined, ndims=50)
DimPlot(data_combined, reduction = "pca", group.by="sample") + scale_color_manual(values=sample_colors) 

data_combined <- RunHarmony(data_combined, group.by.vars = "sample")

# UMAP and tSNE visualization
data_combined <- RunUMAP(object = data_combined, reduction="harmony", dims = 1:30, verbose = FALSE)
DimPlot(data_combined, reduction = "umap", group.by = "sample",raster = F) + scale_color_manual(values=sample_colors)
DimPlot(data_combined, reduction = "umap", group.by = "sample", split.by="sample",ncol = 4,raster = F) + scale_color_manual(values=sample_colors)
data_combined <- RunTSNE(object = data_combined, reduction = "harmony", dims = 1:30)
DimPlot(data_combined, reduction = "tsne", group.by = "sample",raster = F) + scale_color_manual(values=sample_colors)
DimPlot(data_combined, reduction = "tsne", group.by = "sample", split.by="sample",ncol = 4,raster = F) + scale_color_manual(values=sample_colors)

# cell clustering
data_combined <- FindNeighbors(data_combined, reduction = "harmony", dims = 1:30)
data_combined <- FindClusters(data_combined, resolution = 0.4)
cluster_colors <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                    RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                    RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                    RColorBrewer::brewer.pal(n = 8, name = "Dark2"))[-c(6,11)]
DimPlot(data_combined, reduction = "umap", label = T, label.size = 5,raster = F) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "umap", label = T, split.by="sample",label.size = 5,ncol=4,raster = F) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "tsne", label = T, label.size = 5,raster = F) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "tsne", label = T, split.by="sample",label.size = 5,ncol = 4,raster = F) + scale_color_manual(values=cluster_colors)

### Find cell cluster DE markers ###
DefaultAssay(data_combined) <- "RNA"

# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_combined, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
FindDE_all_cluster$cluster <- paste0("C", FindDE_all_cluster$cluster)
FindDE_all_cluster$link <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", rownames(FindDE_all_cluster))

# filtering p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]
# filter RPLXX (Ribosomal protein) and MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

# top10 DE gene visualization
DefaultAssay(data_combined) <- "RNA"
library(dplyr)
library(viridis)
library(RColorBrewer)

top10 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(data_combined, features = as.vector(top10$gene), group.colors=cluster_colors)
# top5 gene dot plot
top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(object = data_combined, features=unique(top5$gene), assay="RNA",col.min = 0,col.max = 1) + scale_color_viridis(direction=1) + theme(axis.text.x = element_text(angle = 45,hjust=1))

### rename cell idents ###
data_combined <- RenameIdents(object = data_combined,
            "19"="Spermatogonia", "7"="Spermatogonia", "18"="Spermatocytes", 
            "21"="RSs", "3"="RSs", "12"="RSs", "4"="ESs/Sperm", "5"="ESs/Sperm","9"="ESs/Sperm", "16"="ESs/Sperm", 
            "0"="Sertoli","6"="Sertoli", "11"="Sertoli", "23"="Sertoli",
            "1"="Leydig","22"="Leydig","8"="Myoid", "10"="Myoid", "24"="Myoid", 
            "20"="Endothelial","2"="Endothelial","15"="Endothelial","14"="Endothelial","17"="Endothelial",
            "13"="Immune")
data_combined$CellType = Idents(data_combined)
celltype_colors <- c("#A6CDE2","#1E78B4","#F59899","#E11E26","#F47E1F","#6A3E98","#FCBF6E","#CAB2D6","#B15928")

DimPlot(data_combined, reduction = "umap", label = F, label.size = 5, raster = F,repel = T) + scale_color_manual(values=celltype_colors)
DimPlot(data_combined, reduction = "umap", label = F, split.by="sample",label.size = 5,ncol=4, raster = F) + scale_color_manual(values=celltype_colors)
DimPlot(data_combined, reduction = "tsne", label = F, label.size = 5, raster = F,repel = T) + scale_color_manual(values=celltype_colors)
DimPlot(data_combined, reduction = "tsne", label = F, split.by="sample",label.size = 5,ncol=4, raster = F) + scale_color_manual(values=celltype_colors)

# marker gene violin plot
myfeatures <- c("STK31","BUD23","SPDYA","LYAR","ACRV1","TNP2","HMGB4","TPPP2","AMH","CLU","DCN","TSHZ2","ACTA2","TAGLN","VWF","ADGRL4","PTPRC","LYZ")
VlnPlot(object = data_combined, assay= "RNA", features=myfeatures, stack=T, 
        flip=F, cols=celltype_colors, fill.by="ident")

### Find cell type DE markers ###
DefaultAssay(data_combined) <- "RNA"
# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_combined, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
# filtering p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]
# filter RPLXX (Ribosomal protein) and MT-XXX (mitchrondrial genes)
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

# top10 DE gene visualization
DefaultAssay(data_combined) <- "RNA"
library(dplyr)
library(viridis)
library(RColorBrewer)

top10 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(data_combined, features = as.vector(top10$gene), group.colors=cluster_colors, draw.lines = F) 

# top5 gene dot plot
top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(object = data_combined, features=unique(top5$gene), assay="RNA", cols = "RdBu")

#### DEs GO and KEGG enrichment ###
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# ID
gene_df = bitr(FindDE_all_cluster$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
head(gene_df)
#   SYMBOL ENTREZID
# 1   AQP3      360
# 2  KRT6C   286887
# 3  LYPD3    27076
data <- merge(FindDE_all_cluster,gene_df,by.x="gene",by.y="SYMBOL")

### GO BP enrichment
formula_res <- compareCluster(ENTREZID~cluster,
                              data=data,
                              OrgDb = "org.Hs.eg.db",
                              ont = "BP", 
                              readable = T,
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              fun="enrichGO")
head(as.data.frame(formula_res))
eggosim <- simplify(formula_res, cutoff = 0.8, by = "p.adjust",
                    select_fun = min, measure = "Wang")
dotplot(eggosim, showCategory=6) + theme_bw(base_size=15) + coord_flip() + theme(axis.text.x = element_text(angle = 60,hjust=1))

### sample-cluster stats ###
sample_stats <- table(data_combined$sample, data_combined$CellType)
library(reshape2)
sample_stats <- melt(sample_stats)
colnames(sample_stats) <- c("Sample","Celltype","Value")

ggplot(sample_stats,aes(Celltype,Value,fill=Sample)) + 
  geom_bar(stat="identity",position="fill") + 
  geom_hline(yintercept=0.5,color="black",linetype=2) + 
  theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_fill_manual(values=sample_colors)
ggplot(sample_stats,aes(Sample,Value,fill=Celltype)) + 
  geom_bar(stat="identity",position="fill") + 
  geom_hline(yintercept=0.5,color="black",linetype=2) + 
  theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_fill_manual(values=celltype_colors)
