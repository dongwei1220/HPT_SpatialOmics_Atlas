### Data analysis for scRNA-seq data in combined H samples
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

### load H data
data_A2H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/A2H_run/outs/filtered_feature_bc_matrix/","A2H")
data_A3H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/A3H_run/outs/filtered_feature_bc_matrix/","A3H")
data_B1H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/B1H_run/outs/filtered_feature_bc_matrix/","B1H")
data_B2H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/B2H_run/outs/filtered_feature_bc_matrix/","B2H")
data_C2H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/C2H_run/outs/filtered_feature_bc_matrix/","C2H")
data_C3H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/C3H_run/outs/filtered_feature_bc_matrix/","C3H")
data_D1H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/D1H_run/outs/filtered_feature_bc_matrix/","D1H")
data_D2H <- ReadData10X("/mnt/hpc/users/data/2results/2025/01_scRNA/D2H_run/outs/filtered_feature_bc_matrix/","D2H")
data_A2H
# An object of class Seurat 
# 26057 features across 15315 samples within 1 assay 
# Active assay: RNA (26057 features, 0 variable features)

data_A2H@assays$RNA@counts[1:5,1:5]
#5 x 5 sparse Matrix of class "dgCMatrix"
#             AAACCCAAGAGGTTTA-1 AAACCCAAGATCCAAA-1 AAACCCACAATACCTG-1
#RPS6KA2                       .                  5                 16
#LOC101866160                  .                  .                  .
#rna76575                      .                  .                  .

### combine all samples
data_combined <- merge(x=data_A2H,y=c(data_A3H,data_B1H,data_B2H,data_C2H,data_C3H,data_D1H,data_D2H),
                        add.cell.ids=c("A2H","A3H","B1H","B2H","C2H","C3H","D1H","D2H"))
dataName="H_combined"

feats <- c("nFeature_RNA", "nCount_RNA")
sample_colors <- c("#1f78b4","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#cab2d6","#b15928","#ff7f00","#6a3d9a","#3D3D3D", "#B4B883")
VlnPlot(data_combined, group.by = "sample", features = feats, ncol = 3, pt.size = 0)
FeatureScatter(data_combined, "nCount_RNA", "nFeature_RNA", group.by = "sample", pt.size = 0.5) + scale_color_manual(values=sample_colors)
# Filter cells
data_combined <- subset(data_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 7000)

### integration with harmony ###
library(harmony)

data_combined <- NormalizeData(data_combined)
data_combined = FindVariableFeatures(data_combined, verbose = F)
data_combined = ScaleData(data_combined, features = rownames(data_combined), verbose = F)

data_combined <- RunPCA(data_combined, npcs = 50, verbose = FALSE)
ElbowPlot(data_combined, ndims=50)
DimPlot(data_combined, reduction = "pca", group.by="sample") + scale_color_manual(values=sample_colors) 

data_combined <- RunHarmony(data_combined, group.by.vars = "sample")
data_combined
# An object of class Seurat 
# 28136 features across 106669 samples within 1 assay 
# Active assay: RNA (28136 features, 2000 variable features)
#  2 dimensional reductions calculated: pca, harmony

# UMAP and tSNE visualization
data_combined <- RunUMAP(object = data_combined, reduction="harmony", dims = 1:30, verbose = FALSE)
DimPlot(data_combined, reduction = "umap", group.by = "sample",raster = F) + scale_color_manual(values=sample_colors)
DimPlot(data_combined, reduction = "umap", group.by = "sample", split.by="sample",ncol = 4, raster = F) + scale_color_manual(values=sample_colors)
data_combined <- RunTSNE(object = data_combined, reduction = "harmony", dims = 1:30)
DimPlot(data_combined, reduction = "tsne", group.by = "sample",raster = F) + scale_color_manual(values=sample_colors)
DimPlot(data_combined, reduction = "tsne", group.by = "sample", split.by="sample",ncol = 4, raster = F) + scale_color_manual(values=sample_colors)

# cell clustering
data_combined <- FindNeighbors(data_combined, reduction = "harmony", dims = 1:30)
data_combined <- FindClusters(data_combined, resolution = 0.4)

cluster_colors <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                    RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                    RColorBrewer::brewer.pal(n = 8, name = "Set2"),
                    RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
                    RColorBrewer::brewer.pal(n = 8, name = "Accent"))[-c(6,11)]
DimPlot(data_combined, reduction = "umap", label = T, label.size = 5,raster = F) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "umap", label = T, split.by="sample",label.size = 5,ncol=4,raster = F) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "tsne", label = T, label.size = 5,raster = F) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "tsne", label = T, split.by="sample",label.size = 5,ncol = 4,raster = F) + scale_color_manual(values=cluster_colors)

### Find cell cluster DE markers ###
DefaultAssay(data_combined) <- "RNA"

# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_combined, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
FindDE_all_cluster$cluster <- paste0("C", FindDE_all_cluster$cluster)

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
DoHeatmap(data_combined, features = as.vector(top10$gene), group.colors=cluster_colors,draw.lines = F) + scale_fill_viridis()
# top5 gene dot plot
top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(object = data_combined, features=unique(top5$gene), assay="RNA") + scale_color_viridis(direction=1)
DotPlot(object = data_combined, features = c("GNRH1","GNRH2","GNRHR","KISS1","KISS1R","TAC1"), assay="RNA") + scale_color_viridis(direction=1) + theme(axis.text.x = element_text(angle = 45,hjust=1))

### rename cell idents ###
data_combined <- RenameIdents(object = data_combined,
            '2'='NEU','4'='NEU','6'='NEU','7'='NEU','8'='NEU','9'='NEU','10'='NEU','11'='NEU',
            '13'='NEU','17'='NEU','18'='NEU','25'='NEU','26'='NEU','27'='NEU',
            '5'='AST','23'='AST','19'='AST',
            '3'='MIC','30'='MIC','22'='MIC','28'='MIC',
            '0'='OLI','31'='OLI','1'='OPC','21'='OPC','29'='OPC','32'='OPC','24'='Immune',
            "12"="VAS",'15'='TANY','20'='EPEN',
            '14'='Unknown','16'='Unknown')
data_combined$CellType = Idents(data_combined)

celltype_colors <- c("#D52126","#88CCEE","#117733","#2F8AC4","#99C945","#89288F","#CC61B0","#E68316","#F97B72","#DDCC77")

DimPlot(data_combined, reduction = "umap", label = F, label.size = 5, raster = F) + scale_color_manual(values=celltype_colors)
DimPlot(data_combined, reduction = "umap", label = F, split.by="sample",label.size = 5,ncol=4, raster = F) + scale_color_manual(values=celltype_colors)
DimPlot(data_combined, reduction = "tsne", label = F, label.size = 5, raster = F) + scale_color_manual(values=celltype_colors)
DimPlot(data_combined, reduction = "tsne", label = F, split.by="sample",label.size = 5,ncol=4, raster = F) + scale_color_manual(values=celltype_colors)

# marker gene violin plot
myfeatures <- c("STMN2","SYT1","SNAP25","SLC17A6","CELF4","SLC1A2","SLC4A4","NHSL1","C1QA","ARHGAP15","PLXDC2","PLP1","MBP","SLC24A2","PDGFRA","VCAN","SOX6",
"PTPRC","ITGA4","CD44","CEMIP","DCN","FN1","COL25A1","PDE7B","SLCO1C1","DTHD1","DNAH12","DYNLRB2")
VlnPlot(object = data_combined, assay= "RNA", features=myfeatures, stack=T, flip=F, cols=celltype_colors, fill.by="ident")

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
DotPlot(object = data_combined, features=unique(top5$gene), assay="RNA", cols = "RdBu",dot.min=0.1) 

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