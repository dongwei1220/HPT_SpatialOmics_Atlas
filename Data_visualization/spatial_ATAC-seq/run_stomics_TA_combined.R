### TA combined sample data analysis ###
library(ArchR)
llibrary(Seurat)
library(ggplot2)

addArchRThreads(threads = 16) 

### load reference genome ###
genomeAnnotation <- readRDS("./refGenome/Mfascicularis.T2TMFA8.genomeAnnotation.rds")
geneAnnotation <- readRDS("./refGenome/Mfascicularis.T2TMFA8.geneAnnotation.rds")
genomeAnnotation
# List of length 3
# names(3): genome chromSizes blacklist
geneAnnotation
# List of length 3
# names(3): genes exons TSS

# load scATAC-seq data
path_root <- "/mnt/hpc/users/data/2results/03_ATAC/"

### step1. 读取fragment输入文件并构建ArchRProject对象
# 读取fragment文件
inputFiles <- grep("outs/fragments.tsv.gz",list.files("/mnt/hpc/users/data/2results/03_ATAC/",recursive=T,pattern="fragments.tsv.gz$",full.names=T),value=T)
names(inputFiles) <- gsub("/outs/fragments.tsv.gz", "", 
                          grep("outs/fragments.tsv.gz",list.files("/mnt/hpc/users/data/2results/03_ATAC/",recursive=T,pattern="fragments.tsv.gz$"),value=T))
names(inputFiles) <- gsub("^.*/","",names(inputFiles))

### selected T samples 
T_samples <- c("A1TA4","A2TA4","A3TA4","B1TA4","B2TA4","C1TA4","C2TA4","C3TA2","D2TA5","D3TA5")
inputFiles_selected <- inputFiles[T_samples]

# addArchRChrPrefix(chrPrefix = FALSE)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles_selected,
  sampleNames = names(inputFiles_selected),
  minTSS = 0, 
  minFrags = 0,
  maxFrags = 1e+12,
  minFragSize = 0,
  maxFragSize = 3000,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  QCDir = "QualityControl",
  force = F,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

# 创建ArchRProject
projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  outputDirectory = paste0("/spatial/run_stomics/ArchR/T2T/TA_combined"),
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme1

df <- getCellColData(projHeme1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")
p

p1 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  baseSize = 7,
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 7,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  baseSize = 7,
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = projHeme1,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 7,
  addBoxPlot = TRUE
)

plotPDF(p1,p2,p3,p4, 
        name = "QC-Sample-Statistics.pdf", 
        ArchRProj = projHeme1, 
        addDOC = FALSE, 
        width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = projHeme1)
p2 <- plotTSSEnrichment(ArchRProj = projHeme1)

plotPDF(p1,p2, 
        name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = projHeme1, 
        addDOC = FALSE, 
        width = 5, height = 5)

library(Seurat)

data_T_combined <- readRDS("data_RNA_TR_combined.rds")
data_T_combined$cellNames <- Cells(data_T_combined)

cellsPass <- intersect(projHeme1$cellNames,data_T_combined$cellNames)
projHeme2 <- projHeme1[cellsPass, ]

projHeme2 <- addCellColData(projHeme2,data=as.vector(data_T_combined$CellType[gsub("A[0-9]#","_",cellsPass)]),
                            name="CellType",cells=cellsPass,force=TRUE)
projHeme2$CellType <- factor(projHeme2$CellType,levels = levels(data_T_combined))

df <- getCellColData(projHeme2, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "QC-Sample-TSS-vs-Frags_post.pdf", 
        ArchRProj = projHeme2, addDOC = FALSE)

p1 <- plotGroups(
  ArchRProj = projHeme2,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  baseSize = 7,
  plotAs = "ridges"
)
p2 <- plotGroups(
  ArchRProj = projHeme2,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 7,
  addBoxPlot = TRUE
)
p3 <- plotGroups(
  ArchRProj = projHeme2,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  baseSize = 7,
  plotAs = "ridges"
)
p4 <- plotGroups(
  ArchRProj = projHeme2,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  baseSize = 7,
  addBoxPlot = TRUE
)

plotPDF(p1,p2,p3,p4, 
        name = "QC-Sample-Statistics_post.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, 
        width = 4, height = 4)

p1 <- plotFragmentSizes(ArchRProj = projHeme2)
p2 <- plotTSSEnrichment(ArchRProj = projHeme2)

plotPDF(p1,p2, 
        name = "QC-Sample-FragSizes-TSSProfile_post.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, 
        width = 5, height = 5)

set.seed(2021)
library(BSgenome.Mfascicularis.NCBI.T2TMFA8)
                                                
projHeme3 <- addGroupCoverages(ArchRProj = projHeme2, groupBy = "CellType")

getChromLengths(projHeme2)
#      chr1      chr2      chr3      chr4      chr5      chr6      chr7      chr8 
# 234122563 203129947 200656507 173030664 194799049 188935743 176059691 154933826 
#      chr9     chr10     chr11     chr12     chr13     chr14     chr15     chr16 
# 142802635 120451489 140374497 138599521 122129872 133402803 128608738  89758255 
#     chr17     chr18     chr19     chr20      chrX 
# 104623588  81929616  66382070  86056530 162126771 
sum(getChromLengths(projHeme2))
# [1] 3042914375

projHeme3 <- addReproduciblePeakSet(
  ArchRProj = projHeme3,
  groupBy = "CellType",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  genomeSize = 3042914375,
  pathToMacs2 = pathToMacs2
)

# peak summary
summary <- projHeme3@peakSet@metadata$PeakCallSummary
summary
#               Group     Var1   Freq
#1         UnionPeaks   Distal 30.444
#2         UnionPeaks   Exonic  6.714
#3         UnionPeaks Intronic 50.491
#4         UnionPeaks Promoter 17.679

peakSets <- getPeakSet(projHeme3)
projHeme4 <- addPeakMatrix(projHeme3)

projHeme4 <- addIterativeLSI(
  ArchRProj = projHeme4,
  useMatrix = "PeakMatrix",
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  force = T,
  varFeatures = 50000,
  dimsToUse = 1:20
)

projHeme4 <- addHarmony(
    ArchRProj = projHeme4,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    force = TRUE,
    groupBy = "Sample"
)
projHeme4 <- addUMAP(
    ArchRProj = projHeme4, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    force = T,
    metric = "cosine"
)
projHeme4 <- addClusters(
    input = projHeme4,
    # reducedDims = "IterativeLSI",
    reducedDims = "Harmony", 
    method = "Seurat",
    name = "Clusters",
    force = TRUE,
    resolution = 0.4
)

p1 <- plotEmbedding(ArchRProj = projHeme4, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    size = 1,
                    baseSize = 12,
                    embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = projHeme4, 
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    size = 1,
                    baseSize = 15,
                    discreteSet = "stallion",
                    embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = projHeme4, 
                    colorBy = "cellColData", 
                    name = "CellType", 
                    size = 1,
                    baseSize = 15,
                    discreteSet = "stallion",
                    embedding = "UMAPHarmony")
                    
plotPDF(p1,p2,p3, name = "Plot-UMAPHarmony-Sample-Clusters.pdf", 
        ArchRProj = projHeme4, addDOC = FALSE, 
        width = 5, height = 5)

projHeme4 <- addImputeWeights(projHeme4)

markerGenes <- c("AMH","CLU","MYH11","TNP2")

p <- plotEmbedding(
  ArchRProj = projHeme4, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPHarmony",
  size = 0.6,
  imputeWeights = getImputeWeights(projHeme4)
)
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = projHeme4, 
        addDOC = FALSE, width = 5, height = 5)

# track plotting
p <- plotBrowserTrack(
  ArchRProj = projHeme4, 
  groupBy = "CellType", 
  geneSymbol = markerGenes, 
  geneAnnotation= geneAnnotation,
  upstream = 50000,
  downstream = 30000
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projHeme4, 
        addDOC = FALSE, width = 6, height = 6)
        
