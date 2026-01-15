### CTR_combined sample data analysis ###
library(AnnotationHub)
library(ensembldb)
library(Seurat)
library(ggplot2)

options(future.globals.maxSize = 3e+09)

### get barcode information
barcode_Numb <- "192"

## barcode file
BCA_file <- paste0("/spatial/barcode_96_192/barcode_", barcode_Numb, "_A.txt")
BCB_file <- paste0("/spatial/barcode_96_192/barcode_", barcode_Numb, "_B.txt")

BCA_Data <- read.table(BCA_file, sep="\t", header=FALSE)
BCB_Data <- read.table(BCB_file, sep="\t", header=FALSE)

BCA_Data["index_num"] <- rownames(BCA_Data)
names(BCA_Data) <- c("Barcode", "index_num")

BCB_Data["index_num"] <- rownames(BCB_Data)
names(BCB_Data) <- c("Barcode", "index_num")

head(BCA_Data)
#   Barcode index_num
# 1 AACGTGAT         1
# 2 CAGATCTG         2
# 3 AACAACCA         3
head(BCB_Data)
#    Barcode index_num
# 1 AACGTGAT         1
# 2 CAGATCTG         2
# 3 AACAACCA         3

# load expresion matrix
path_root <- "/mnt/hpc/users/data/2results/count/"
# selected sample list
sample_list <- c("A1CTR2","A2CTR2","A3CTR2","B1CTR4","B2CTR4","C1CTR4","C2CTR4","C3CTR2","D1CTR2","D2CTR2","D3CTR4")

for(dataName in sample_list){

counts_RNA_path <- file.path(path_root, paste0(dataName, "_stdata.under.tissue.expression-matrix.tsv"))
counts_RNA <- read.table(counts_RNA_path, sep="\t", header=TRUE)

### barcode convert ###
for(line_num in c(1:dim(counts_RNA)[1])){

    BC_INFO <- strsplit(counts_RNA[line_num,"X"], "x")[[1]]
    BC_A <- BC_INFO[1]
    BC_B <- BC_INFO[2]

    BC_A <- BCA_Data[BCA_Data["index_num"] == BC_A, "Barcode"]
    BC_B <- BCB_Data[BCB_Data["index_num"] == BC_B, "Barcode"]
    # combine barcode (BC_B-BC_A-1)
    BC_Merge <- paste0(BC_B, BC_A, "-1")

    counts_RNA[line_num, "X"] <- BC_Merge

}

rownames(counts_RNA) <- counts_RNA[, 1]
counts_RNA <- counts_RNA[, -1]
counts_RNA <- t(counts_RNA)

# create a Seurat object containing the RNA adata
data_RNA <- CreateSeuratObject(counts = counts_RNA, assay = "RNA", project = dataName)

### SCTransform normaization
data_RNA <- SCTransform(data_RNA,vars.to.regress = "nFeature_RNA", verbose = FALSE)
data_RNA

saveRDS(data_RNA, paste0("data_RNA_",dataName,".rds"))

}

### Load seurat object data
sample_list <- c("A1CTR2","A2CTR2","A3CTR2","B1CTR4","B2CTR4","C1CTR4","C2CTR4","C3CTR2","D1CTR2","D2CTR2","D3CTR4")
data_A1CT <- readRDS("./data_RNA_A1CTR2.rds")
data_A2CT <- readRDS("./data_RNA_A2CTR2.rds")
data_A3CT <- readRDS("./data_RNA_A3CTR2.rds")
data_B1CT <- readRDS("./data_RNA_B1CTR4.rds")
data_B2CT <- readRDS("./data_RNA_B2CTR4.rds")
data_C1CT <- readRDS("./data_RNA_C1CTR4.rds")
data_C2CT <- readRDS("./data_RNA_C2CTR4.rds")
data_C3CT <- readRDS("./data_RNA_C3CTR2.rds")
data_D1CT <- readRDS("./data_RNA_D1CTR2.rds")
data_D2CT <- readRDS("./data_RNA_D2CTR2.rds")
data_D3CT <- readRDS("./data_RNA_D3CTR4.rds")
DefaultAssay(data_A1CT) <- "SCT"
DefaultAssay(data_A2CT) <- "SCT"
DefaultAssay(data_A3CT) <- "SCT"
DefaultAssay(data_B1CT) <- "SCT"
DefaultAssay(data_B2CT) <- "SCT"
DefaultAssay(data_C1CT) <- "SCT"
DefaultAssay(data_C2CT) <- "SCT"
DefaultAssay(data_C3CT) <- "SCT"
DefaultAssay(data_D1CT) <- "SCT"
DefaultAssay(data_D2CT) <- "SCT"
DefaultAssay(data_D3CT) <- "SCT"

data_list <- list(data_A1CT,data_A2CT,data_A3CT,data_B1CT,data_B2CT,data_C1CT,data_C2CT,data_C3CT,data_D1CT,data_D2CT,data_D3CT)
var.features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)

data_RNA <- merge(x=data_A1CT,y=c(data_A2CT,data_A3CT,data_B1CT,data_B2CT,data_C1CT,data_C2CT,data_C3CT,data_D1CT,data_D2CT,data_D3CT),
                  add.cell.ids=c("A1CT","A2CT","A3CT","B1CT","B2CT","C1CT","C2CT","C3CT","D1CT","D2CT","D3CT"))

dataName <- "CTR_combined"

# Visualize QC metrics as a violin plot
sample_colors <- as.vector(ArchR::ArchRPalettes$stallion)[-c(2,6)]
VlnPlot(data_RNA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3,group.by = "orig.ident",pt.size = 0)
FeatureScatter(data_RNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "orig.ident") + scale_color_manual(values=sample_colors)

# cell filtering
data_RNA <- subset(data_RNA, subset = nFeature_RNA < 10000 & nFeature_RNA > 20)
sample_colors <- as.vector(ArchR::ArchRPalettes$stallion)[-c(2,6)]
VlnPlot(data_RNA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3,group.by = "orig.ident",pt.size = 0) 

VariableFeatures(data_RNA) <- var.features
data_RNA <- RunPCA(data_RNA, verbose = FALSE)
ElbowPlot(data_RNA, ndims = 50)

library(harmony)
data_RNA <- RunHarmony(data_RNA, assay.use="SCT", group.by.vars = "orig.ident")
data_RNA <- RunUMAP(data_RNA, reduction = "harmony", dims = 1:20)
data_RNA <- FindNeighbors(data_RNA, reduction = "harmony", dims = 1:20)
DefaultAssay(data_RNA) <- "SCT"
data_RNA <- FindClusters(data_RNA, resolution = 0.7)

cluster_colors <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                    RColorBrewer::brewer.pal(n = 12, name = "Set3"),
                    RColorBrewer::brewer.pal(n = 8, name = "Set2"))[-c(6,11)]
DimPlot(data_RNA, reduction = "umap", label=TRUE, label.size=7) + scale_color_manual(values=cluster_colors)
DimPlot(data_RNA, reduction = "umap", group.by = "orig.ident",label=F, label.size=7) + scale_color_manual(values=sample_colors)
DimPlot(data_RNA, reduction = "umap", split.by = "orig.ident",label=TRUE, label.size=7,ncol=4) + scale_color_manual(values=cluster_colors)

### MAGIC data imputation ###
library(Rmagic)

# Run MAGIC and reset the active assay
data_RNA <- magic(data_RNA,assay="RNA")
# [1] "Added MAGIC output to MAGIC_RNA. To use it, pass assay='MAGIC_RNA' to downstream methods or set seurat_object@active.assay <- 'MAGIC_RNA'."

### Find DE markers ###
DefaultAssay(data_RNA) <- "RNA"
data_RNA <- NormalizeData(data_RNA)

# DE analysis for every cluster
FindDE_all_cluster <- FindAllMarkers(data_RNA, only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
FindDE_all_cluster$cluster <- paste0("C", FindDE_all_cluster$cluster)
# filtering p val adjust < 0.05
FindDE_all_cluster <- FindDE_all_cluster[FindDE_all_cluster$p_val_adj < 0.05, ]

# top10 DE gene visualization
DefaultAssay(data_RNA) <- "SCT"
library(dplyr)
library(viridis)
library(RColorBrewer)

top10 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(data_RNA, features = as.vector(top10$gene), group.colors=cluster_colors,raster = F) + scale_fill_viridis()

# marker gene feature plot
FeaturePlot(data_RNA,features = top10$gene[1:6],order = T,min.cutoff = "q5",max.cutoff = "q95",cols=c("gray","red"),ncol = 3)
# marker gene violin plot
VlnPlot(object = data_RNA, assay= "RNA", features=top10$gene[1:6], stack=T, 
        flip=T, cols=cluster_colors[as.numeric(levels(data_RNA))+1], fill.by="ident") 
# top5 gene dot plot
top5 <- FindDE_all_cluster %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(object = data_RNA, features=unique(top5$gene), assay="RNA", col.min=0,col.max=1) + scale_color_viridis(direction=1) + theme(axis.text.x = element_text(angle = 45,hjust=1))

### rename cell idents ###
data_RNA <- RenameIdents(object = data_RNA,
            "0"="Somatotrope", "2"="Somatotrope", "6"="Lactotrope", "10"="Corticocrope",
            "4"="Gonado/Folliculo", "17"="Thyrotrope","3"="Melanotrope", "12"="Microglia",
            "1"="Pituicytes","5"="Pituicytes", "7"="Pituicytes", 
            "14"="Pituicytes","8"="Pituicytes","11"="Endo","15"="Endo",
            "9"="Stromal","16"="Stromal","13"="Unknown")
DotPlot(object = data_RNA, features=features, assay="RNA", cluster.idents=F, col.min=0.0, col.max=1.0, dot.min=0) + scale_color_viridis(direction=1) + RotatedAxis()
data_RNA$CellType = Idents(data_RNA)

celltype_colors <- c("#89288F","#2F8AC4","#208A42","#D51F26", "#F47D2B","#8A9FD1","#C06CAB","#F37B7D","#90D5E4","#89C75F","#D8A767")
DimPlot(data_RNA, reduction = "umap", label = F, label.size = 5, raster = F,repel = T) + scale_color_manual(values=celltype_colors)
DimPlot(data_RNA, reduction = "umap", label = F, split.by="orig.ident",label.size = 5,ncol=4, raster = F) + scale_color_manual(values=celltype_colors)

for(i in unique(data_RNA$orig.ident)){
    data_sample <- subset(data_RNA, subset= orig.ident == i )
    data_sample <- RenameCells(data_sample,new.names = gsub("^.*_","",Cells(data_sample)))
    pdf(paste0(dataName,"_tissue/Spatial_DimPlot_umap_",i,"_by_celltype.pdf"),height=7,width=9.7)
    p <- spatial_DimPlot(data_sample, reduction = "umap", 192) + scale_color_manual(values=celltype_colors) + 
                        labs(title = i) + guides(color = guide_legend(override.aes = list(size=4)))
    print(p)
    dev.off()
}

### sample-cluster stats ###
sample_stats <- table(data_RNA$orig.ident, data_RNA$CellType)
library(reshape2)
sample_stats <- melt(sample_stats)
colnames(sample_stats) <- c("Sample","Celltype","Value")

sample_colors <- as.vector(ArchR::ArchRPalettes$circus)[-3]
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

data_RNA@meta.data$Stage <- "NA"
data_RNA@meta.data[data_RNA@meta.data$orig.ident %in% c("A1CTR2","A2CTR2","A3CTR2"),"Stage"] <- "Juvenile_Pituitary"
data_RNA@meta.data[data_RNA@meta.data$orig.ident %in% c("B1CTR4","B2CTR4"),"Stage"] <- "Pubertal_Pituitary"
data_RNA@meta.data[data_RNA@meta.data$orig.ident %in% c("C1CTR4","C2CTR4","C3CTR2"),"Stage"] <- "Adult_Pituitary"
data_RNA@meta.data[data_RNA@meta.data$orig.ident %in% c("D1CTR2","D2CTR2","D3CTR4"),"Stage"] <- "Old_Pituitary"
data_RNA$Stage <- factor(data_RNA$Stage,levels = c("Juvenile_Pituitary","Pubertal_Pituitary","Adult_Pituitary","Old_Pituitary"))

sample_stats <- table(data_RNA$Stage, data_RNA$CellType)
library(reshape2)
sample_stats <- melt(sample_stats)
colnames(sample_stats) <- c("Stage","Celltype","Value")

ggplot(sample_stats,aes(Celltype,Value,fill=Stage)) + 
  geom_bar(stat="identity",position="fill") + 
  geom_hline(yintercept=0.5,color="black",linetype=2) + 
  theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_fill_manual(values=sample_colors)
ggplot(sample_stats,aes(Stage,Value,fill=Celltype)) + 
  geom_bar(stat="identity",position="fill") + 
  geom_hline(yintercept=0.5,color="black",linetype=2) + 
  theme_bw(base_size = 16) + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  scale_fill_manual(values=celltype_colors)
