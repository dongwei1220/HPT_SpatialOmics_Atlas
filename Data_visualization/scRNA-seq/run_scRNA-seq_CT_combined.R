### Data analysis for scRNA-seq data in combined CT samples
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

### load CT data
data_A2CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/A2CT_run/outs/filtered_feature_bc_matrix/","A2CT")
data_A3CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/A3CT_run/outs/filtered_feature_bc_matrix/","A3CT")
data_B1CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/B1CT_run/outs/filtered_feature_bc_matrix/","B1CT")
data_B2CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/B2CT_run/outs/filtered_feature_bc_matrix/","B2CT")
data_C2CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/C2CT_run/outs/filtered_feature_bc_matrix/","C2CT")
data_C3CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/C3CT_run/outs/filtered_feature_bc_matrix/","C3CT")
data_D1CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/D1CT_run/outs/filtered_feature_bc_matrix/","D1CT")
data_D2CT <- ReadData10X("/mnt/hpc/users/data/2results/scRNA/D2CT_run/outs/filtered_feature_bc_matrix/","D2CT")
data_A2CT
# An object of class Seurat 
# 24701 features across 8222 samples within 1 assay 
# Active assay: RNA (24701 features, 0 variable features)
data_A2CT@assays$RNA@counts[1:5,1:5]
#5 x 5 sparse Matrix of class "dgCMatrix"
#             AAACCCAAGAGGTTTA-1 AAACCCAAGATCCAAA-1 AAACCCACAATACCTG-1
#RPS6KA2                       .                  5                 16
#LOC101866160                  .                  .                  .
#rna76575                      .                  .                  .

### combine all samples
data_combined <- merge(x=data_A2CT,y=c(data_A3CT,data_B1CT,data_B2CT,data_C2CT,data_C3CT,data_D1CT,data_D2CT),
                        add.cell.ids=c("A2CT","A3CT","B1CT","B2CT","C2CT","C3CT","D1CT","D2CT"))
dataName="CT_combined"

# plot QC metrics
feats <- c("nFeature_RNA", "nCount_RNA")
sample_colors <- as.vector(ArchR::ArchRPalettes$stallion)[-c(2,6)]
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
data_combined <- RunUMAP(object = data_combined, reduction="harmony", dims = 1:30, verbose = FALSE)
DimPlot(data_combined, reduction = "umap", group.by = "sample",raster=F) + scale_color_manual(values=sample_colors)

# cell clustering
data_combined <- FindNeighbors(data_combined, reduction = "harmony", dims = 1:30)
data_combined <- FindClusters(data_combined, resolution = 0.5)

cluster_colors <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                    RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                    RColorBrewer::brewer.pal(n = 8, name = "Set2"))[-c(6,11)]
DimPlot(data_combined, reduction = "umap", label = T, label.size = 5,raster=FALSE) + scale_color_manual(values=cluster_colors)
DimPlot(data_combined, reduction = "umap", label = T, split.by="sample",label.size = 5,ncol=4,raster=FALSE) + scale_color_manual(values=cluster_colors)

### Find cell cluster DE markers ###
DefaultAssay(data_combined) <- "RNA"

# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_combined, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
# filtering p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

# top10 DE gene visualization
DefaultAssay(data_combined) <- "RNA"
top10 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(data_combined, features = as.vector(top10$gene), group.colors=cluster_colors)
# top5 gene dot plot
top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(object = data_combined, features=unique(top5$gene), assay="RNA", group.by = "seurat_clusters")

### rename cell idents ###
data_combined <- RenameIdents(object = data_combined,
            "0"="Somatotrope", "3"="Somatotrope","1"="Lactotrope","11"="Lactotrope","19"="Lactotrope", 
            "8"="Corticotrope", "20"="Corticotrope", "6"="Gonadotrope", "16"="Gonadotrope", "10"="Thyrotrope",
            "7"="Melanotrope", "2"="Folliculostellate", "18"="Folliculostellate", "23"="Folliculostellate", 
            "4"="Microglia","15"="Microglia","17"="Microglia","5"="Pituicytes","14"="Pituicytes","22"="Pituicytes",
            "13"="Endo","21"="Endo","9"="Stromal","12"="Immune") 
data_combined$CellType = Idents(data_combined)

celltype_colors <- c("#89288F","#2F8AC4","#208A42","#D51F26","#F47D2B","#8A9FD1","#3BBCA8","#C06CAB","#F37B7D","#90D5E4","#89C75F","#D8A767")
DimPlot(data_combined, reduction = "umap", group.by="CellType", label = F, label.size = 5, raster = F) + scale_color_manual(values=celltype_colors)

### Find cell type DE markers ###
DefaultAssay(data_combined) <- "RNA"

# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_combined, only.pos = T, logfc.threshold = 0.5, min.pct = 0.25)
# filtering
# p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]
tmp <- grep("^RP[SL]",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}
tmp <- grep("^MT-",rownames(FindDE_all_cluster), perl=T)
if(length(tmp) > 0){FindDE_all_cluster <- FindDE_all_cluster[-tmp, ]}

DefaultAssay(data_combined) <- "RNA"
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

sample_colors <- as.vector(ArchR::ArchRPalettes$circus)
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
