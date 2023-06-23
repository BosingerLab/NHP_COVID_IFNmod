
setwd("/data/runs/Analysis/CovTen_Arun/")

"Load all libraries for Analysis"

library(stringr)
library(Seurat)
library(ggplot2)
library(SingleR)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(plotly)
library(ggExtra)
library(kableExtra)
library(knitr)
library(Matrix)
library(biomaRt)
library(mltools)

ref <- BlueprintEncodeData()
y_genes <- as.character(unlist(read.table("/data/runs/Genome_references/macaca_mulatta/mmul10/ensembl_100/Y_genes.txt"),use.names = F))
mt_genes <-  as.character(unlist(read.table("/data/runs/Genome_references/macaca_mulatta/mmul10/ensembl_100/MT_genes.txt"),use.names = F))
ig_genes <- as.character(unlist(read.table("/data/runs/Genome_references/macaca_mulatta/mmul10/ensembl_100/IG_genes.txt"),use.names = F))
tr_genes <- as.character(unlist(read.table("/data/runs/Genome_references/macaca_mulatta/mmul10/ensembl_100/TR_genes.txt"),use.names = F))
prot_coding <- as.character(unlist(read.table("/data/runs/Genome_references/macaca_mulatta/mmul10/ensembl_100/protein_coding_seurat.txt"),
                                 use.names = F))

##Insert Working directories for different groups (Optional, You can give your working directory if you have only 1 super-folder)
group1.dir = "/data/runs/210524_A00945_0114_AHYWLMDSXY/HYWLMDSXY/outs/fastq_path/p21061_Tim"


#Get a list of Sub-Directoties
projectDIR = group1.dir
project.subDIRs = lapply( projectDIR , function(dir){ setNames( list.dirs( dir , recursive = F , full.names = T ) , 
                                                            list.dirs( dir , recursive = F ) ) } )
project.subDIRs = unlist( project.subDIRs )

BALs = project.subDIRs[ !grepl("Lung",names(project.subDIRs) ) ]

#Get Directories with matching pattern (here, CITE for CITE-Seq Counts and GEX for Gene expression matrices)
#BALs_CITE = BALs[ grepl("CITE$",mes(BALs)) ]
BALs_GEX = BALs[ !grepl("Lung$",names(BALs)) ]

#Function to get Gene expression matrices
gex <- function( gex_dir )
{
gex_raw = grep("filtered_feature_bc_matrix",list.dirs( gex_dir , recursive = T , full.names = T ),value=T)

barcode.fn = file.path( gex_raw , "barcodes.tsv.gz" )
features.fn = file.path( gex_raw , "features.tsv.gz" )
matrix.fn = file.path( gex_raw , "matrix.mtx.gz" )

matrx = readMM(file = matrix.fn)
feature.names = read.delim(features.fn, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.fn, header = FALSE, stringsAsFactors = FALSE)
colnames(matrx) = barcode.names$V1
rownames(matrx) = feature.names$V2

return( matrx )
}
gex_matrixList_GEX = lapply( BALs_GEX , gex )

##This part of the code is to intersect Gene expression and HTO/ADT Capture
BAL.groups = list()
Classificationtable = list()


for(i in 4:length(gex_matrixList_GEX)){
  BAL_group1_gex = gex_matrixList_GEX[[i]][grep("^Hash", rownames(gex_matrixList_GEX[[i]]),
                                               value =T,invert=TRUE), ]
  
  HTO_matrix = gex_matrixList_GEX[[i]][grep("^Hash", rownames(gex_matrixList_GEX[[i]]), value =
                                           T), ]
  rowhtos = do.call("rbind" , str_split(rownames(HTO_matrix), pattern = "_"))
  rownames(HTO_matrix) = rowhtos[,1]
  colnames(BAL_group1_gex) = gsub("-.*", "", colnames(BAL_group1_gex))
  colnames(HTO_matrix) = gsub("-.*", "", colnames(HTO_matrix))
  joint.bcs = intersect(colnames(BAL_group1_gex), colnames(HTO_matrix))
  
  #Subset R and HTO counts by joint cell barcodes
  BAL_group1_gex <- BAL_group1_gex[, joint.bcs]
  HTO_matrix <- HTO_matrix[, joint.bcs]
  BAL.hastag <- CreateSeuratObject(counts = BAL_group1_gex,
                                   project = str_split(basename(names(gex_matrixList_GEX[i])),
                                                       pattern = "_count_",n = 2)[[1]][1])
  
  
  cov2_genes <- row.names(BAL.hastag)[grepl("SARS-CoV2",row.names(BAL.hastag))]
  rps_l_genes <- row.names(BAL.hastag)[grepl("^RP[S,L]",row.names(BAL.hastag))]
  BAL.hastag[["percent.cov2"]] <- PercentageFeatureSet(BAL.hastag, features = cov2_genes)
  BAL.hastag[["percent.ig"]] <- PercentageFeatureSet(BAL.hastag, features = ig_genes)
  BAL.hastag[["percent.trgenes"]] <- PercentageFeatureSet(BAL.hastag, features = tr_genes)
  BAL.hastag[["percent.hbb"]]  = PercentageFeatureSet(BAL.hastag, pattern = "^HBB")
  BAL.hastag[["percent.rps"]]  = PercentageFeatureSet(BAL.hastag, pattern = "^RPS")
  BAL.hastag[["percent.rpl"]]  = PercentageFeatureSet(BAL.hastag, pattern = "^RPL")
  BAL.hastag[["percent.mt"]]  = PercentageFeatureSet(BAL.hastag, features = mt_genes)
  BAL.hastag$log10GenesPerUMI <- log10(BAL.hastag$nFeature_RNA) / log10(BAL.hastag$nCount_RNA)
  
  meta <- BAL.hastag@meta.data[,c("percent.hbb","percent.mt","percent.rps",
                                  "percent.rpl","percent.ig","percent.cov2","percent.trgenes")]
  
  
  is_coding <- row.names(BAL.hastag) %in% prot_coding
  non_coding <- row.names(BAL.hastag)[!is_coding]
  filt_genes <- c("HBB",y_genes, mt_genes, rps_l_genes, 
                  ig_genes, tr_genes, cov2_genes, non_coding)
  is_filt <- row.names(BAL.hastag) %in% filt_genes
  counts <- GetAssayData(object = BAL.hastag, slot = "counts")
  Sample <- CreateSeuratObject(counts=counts[!is_filt,],project = str_split(basename(names(gex_matrixList_GEX[i])),
                                                                            pattern = "_count_",n = 2)[[1]][1])
  
  Sample@meta.data <- transform(merge(Sample@meta.data, meta, by=0, Sample.x=TRUE), 
                                row.names=Row.names, Row.names=NULL)
  
  Sample[["Animal_Day"]]  = str_split(Sample@project.name,
                               pattern = "_",n = 3)[[1]][3]
  
  HTO_matrix_subset = HTO_matrix[which(rowSums2(HTO_matrix) > 10000),]
  
  if( is.vector(HTO_matrix_subset) )
    {
    HTO_matrix_subset = t( as.data.frame( HTO_matrix_subset ) )
    hash_name = row.names(HTO_matrix)[which(rowSums2(HTO_matrix) > 10000)]
    rownames(HTO_matrix_subset) = hash_name
    }
  # Add HTO data as a new assay independent from R
  
  Sample[["HTO"]] <- CreateAssayObject(counts = HTO_matrix_subset)
  
    # Normalize R data with log normalization
  Sample <- NormalizeData(Sample)
  # Find and scale variable features
  Sample <-
    FindVariableFeatures(Sample, selection.method = "mean.var.plot")
  Sample <-
    ScaleData(Sample, features = VariableFeatures(Sample))
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  Sample <-
    NormalizeData(Sample,
                  assay = "HTO",
                  normalization.method = "CLR")
  # Global classification results
  # table(BAL.hashtag_C1$ADT_classification.global)
  
  if(length(rownames(Sample[["HTO"]]) >= 2)){
    
    
    Sample <-
      HTODemux(Sample, assay = "HTO", positive.quantile = 0.99 )
    
    
    Idents(Sample) <- "HTO_classification.global"
    # First, we will remove negative cells from the object
    Classificationtable[[i]] = table(Sample$HTO_classification.global)
    
    Sample <-
      subset(Sample, idents = "Singlet", invert = FALSE)
    
    Sample <-
      ScaleData(Sample,
                features = rownames(Sample),
                verbose = FALSE)
    Sample <-
      RunPCA(Sample,
             features = rownames(Sample),
             approx = FALSE)
    Sample <- FindNeighbors(Sample, reduction = "pca", dims = 1:10)
    Sample <- FindClusters(Sample, resolution = 0.6, verbose = FALSE)
    
    Sample <-
      RunTSNE(Sample,
              dims = 1:10, check_duplicates = F)
    Sample <-
      RunUMAP(Sample,
              dims = 1:10)
    
    Sample.sce <- as.SingleCellExperiment(Sample)
    
    "Convert the seurat object to 
          singlecellexperiment"
    
    Sample_singler_bp = SingleR(test = Sample.sce, ref = ref,
                                labels = ref$label.main, fine.tune = T,
                                prune = T)
    Sample$singleR_main = Sample_singler_bp$labels
    Sample$singleR_pruned = Sample_singler_bp$pruned.labels
    
    
    BAL.groups[[i]] = Sample
  }
}

SingleSamples = list()
##Single Samples
gex_matrixList_GEX_v1 = gex_matrixList_GEX[c(1,3)]
for(i in 1:length(gex_matrixList_GEX_v1)){
  BAL_group1_gex = gex_matrixList_GEX_v1[[i]][grep("^Hash", rownames(gex_matrixList_GEX_v1[[i]]),
                                                   value =T,invert=TRUE), ]
  
  HTO_matrix = gex_matrixList_GEX_v1[[i]][grep("^Hash", rownames(gex_matrixList_GEX_v1[[i]]), value =
                                                 T), ]
  rowhtos = do.call("rbind" , str_split(rownames(HTO_matrix), pattern = "_"))
  rownames(HTO_matrix) = rowhtos[,1]
  colnames(BAL_group1_gex) = gsub("-.*", "", colnames(BAL_group1_gex))
  colnames(HTO_matrix) = gsub("-.*", "", colnames(HTO_matrix))
  joint.bcs = intersect(colnames(BAL_group1_gex), colnames(HTO_matrix))
  
  #Subset R and HTO counts by joint cell barcodes
  BAL_group1_gex <- BAL_group1_gex[, joint.bcs]
  HTO_matrix <- HTO_matrix[, joint.bcs]
  BAL.hastag <- CreateSeuratObject(counts = BAL_group1_gex,
                                   project = str_split(basename(names(gex_matrixList_GEX_v1[i])),
                                                       pattern = "_count_",n = 2)[[1]][1])
  
  
  cov2_genes <- row.names(BAL.hastag)[grepl("SARS-CoV2",row.names(BAL.hastag))]
  rps_l_genes <- row.names(BAL.hastag)[grepl("^RP[S,L]",row.names(BAL.hastag))]
  BAL.hastag[["percent.cov2"]] <- PercentageFeatureSet(BAL.hastag, features = cov2_genes)
  BAL.hastag[["percent.ig"]] <- PercentageFeatureSet(BAL.hastag, features = ig_genes)
  BAL.hastag[["percent.trgenes"]] <- PercentageFeatureSet(BAL.hastag, features = tr_genes)
  BAL.hastag[["percent.hbb"]]  = PercentageFeatureSet(BAL.hastag, pattern = "^HBB")
  BAL.hastag[["percent.rps"]]  = PercentageFeatureSet(BAL.hastag, pattern = "^RPS")
  BAL.hastag[["percent.rpl"]]  = PercentageFeatureSet(BAL.hastag, pattern = "^RPL")
  BAL.hastag[["percent.mt"]]  = PercentageFeatureSet(BAL.hastag, features = mt_genes)
  BAL.hastag$log10GenesPerUMI <- log10(BAL.hastag$nFeature_RNA) / log10(BAL.hastag$nCount_RNA)
  
  meta <- BAL.hastag@meta.data[,c("percent.hbb","percent.mt","percent.rps",
                                  "percent.rpl","percent.ig","percent.cov2","percent.trgenes")]
  
  
  is_coding <- row.names(BAL.hastag) %in% prot_coding
  non_coding <- row.names(BAL.hastag)[!is_coding]
  filt_genes <- c("HBB",y_genes, mt_genes, rps_l_genes, 
                  ig_genes, tr_genes, cov2_genes, non_coding)
  is_filt <- row.names(BAL.hastag) %in% filt_genes
  counts <- GetAssayData(object = BAL.hastag, slot = "counts")
  Sample <- CreateSeuratObject(counts=counts[!is_filt,],project = str_split(basename(names(gex_matrixList_GEX_v1[i])),
                                                                            pattern = "_count_",n = 2)[[1]][1])
  
  Sample@meta.data <- transform(merge(Sample@meta.data, meta, by=0, Sample.x=TRUE), 
                                row.names=Row.names, Row.names=NULL)
  
  Sample[["Animal_Day"]]  = str_split(Sample@project.name,
                                      pattern = "_",n = 3)[[1]][3]
  
  HTO_matrix_subset = HTO_matrix[which(rowSums2(HTO_matrix) > 10000),]
  
  if( is.vector(HTO_matrix_subset) )
  {
    HTO_matrix_subset = t( as.data.frame( HTO_matrix_subset ) )
    hash_name = row.names(HTO_matrix)[which(rowSums2(HTO_matrix) > 10000)]
    rownames(HTO_matrix_subset) = hash_name
  }
  # Add HTO data as a new assay independent from R
  
  Sample[["HTO"]] <- CreateAssayObject(counts = HTO_matrix_subset)
  
  # Normalize R data with log normalization
  Sample <- NormalizeData(Sample)
  # Find and scale variable features
  Sample <-
    FindVariableFeatures(Sample, selection.method = "mean.var.plot")
  Sample <-
    ScaleData(Sample, features = VariableFeatures(Sample))
  # Normalize HTO data, here we use centered log-ratio (CLR) transformation
  Sample <-
    NormalizeData(Sample,
                  assay = "HTO",
                  normalization.method = "CLR")
  # Global classification results
  # table(BAL.hashtag_C1$ADT_classification.global)
  
    Sample <-
      ScaleData(Sample,
                features = rownames(Sample),
                verbose = FALSE)
    Sample <-
      RunPCA(Sample,
             features = rownames(Sample),
             approx = FALSE)
    Sample <- FindNeighbors(Sample, reduction = "pca", dims = 1:10)
    Sample <- FindClusters(Sample, resolution = 0.6, verbose = FALSE)
    
    Sample <-
      RunTSNE(Sample,
              dims = 1:10, check_duplicates = F)
    Sample <-
      RunUMAP(Sample,
              dims = 1:10)
    
    Sample.sce <- as.SingleCellExperiment(Sample)
    
    "Convert the seurat object to 
          singlecellexperiment"
    
    Sample_singler_bp = SingleR(test = Sample.sce, ref = ref,
                                labels = ref$label.main, fine.tune = T,
                                prune = T)
    Sample$singleR_main = Sample_singler_bp$labels
    Sample$singleR_pruned = Sample_singler_bp$pruned.labels
    
    SingleSamples[[i]] = Sample
}

BAL.groups[[1]] = SingleSamples[[1]]

BAL.groups[[1]]$hash.ID = "Hash1Human"

BAL.groups[[2]] = SingleSamples[[2]]

BAL.groups[[2]]$hash.ID = "Hash2Human"

BAL.groups_new = BAL.groups[c(1,2,4:23)]
#write_rds(BAL.groups,file = "AllSample_withSingleR.rds",version = 2 )



# ##---MERGE THE SAMPLE LIST----## IMPORTANT STEP 2 
# 
BAL_merge <- merge(BAL.groups_new[[1]], y=BAL.groups_new[2:length(BAL.groups_new)],
                   add.cell.ids=c("p21061_s001_RMs8-GR1_BL","p21061_s002_CD68-GR1_BL",
                                  "p21061_s003_RMj8-KV31_BL","p21061_s004_RMs8-CD68_D2",
                                  "p21061_s005_RMj8-KV31_D2", "p21061_s006_RMs8-CD68_D4",
                                  "p21061_s007_RMj8-KV31_D4", "p21061_s008_CD68-KV31_D7",
                                  "p21061_s009_RRc17-RGo9_BL","p21061_s010_RCk17-RRj11_BL",
                                  "p21061_s011_RRc17-RGo9_D2","p21061_s012_RCk17-RRj11_D2",
                                  "p21061_s013_RRc17-RGo9_D4", "p21061_s014_RCk17-RRj11_D4",
                                  "p21061_s015_RGo9-RRj11_D7", "p21061_s016_LG92-LC10_BL",
                                  "p21061_s017_RNi17-RZs14_BL", "p21061_s018_LG92-LC10_D2",
                                  "p21061_s019_RNi17-RZs14_D2","p21061_s020_LG92-LC10_D4",
                                  "p21061_s021_RNi17-RZs14_D4","p21061_s024_LC10-RZs14_D7"))
DefaultAssay(BAL_merge) = "RNA"

BAL_merge$log10GenesPerUMI <- log10(BAL_merge$nFeature_RNA) / log10(BAL_merge$nCount_RNA)

BAL_merge$hashing <- str_match(row.names(BAL_merge@meta.data),
                              "[[:alnum:]]+-[[:alnum:]]+")

BAL_merge@meta.data <- BAL_merge@meta.data %>% separate(hashing,c("Animal1","Animal2"),"-",
                                                        remove = FALSE)

BAL_merge@meta.data <- BAL_merge@meta.data %>% separate(Animal_Day,c("hashing","Day"),"_",
                                                        remove = FALSE)

#######

# Group	Animal	Hash#
# 1	RMs8	1
# 1	CD68	2
# 1	RMj8	3
# 1	KV31	4
# 2	RRc17	1
# 2	RGo9	2
# 2	RCk17	3
# 2	RRj11	4
# 3	LG92	1
# 3	LC10	2
# 3	RNi17	3
# 3	RZs14	4

# Control:
#   RMs8 , CD68 , RRc17 , RGo9 , LG92 , LC10

# IFNant:
#   RMj8, KV31, RCk17, RRj11, RNi17, RZs14

#######

BAL_merge$Treatment = "NA"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash1Human" & BAL_merge$Animal1 == "RMs8" & BAL_merge$Animal2 == "GR1")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal1 == "CD68" & BAL_merge$Animal2 == "GR1")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash3Human" & BAL_merge$Animal1 == "RMj8")] = "IFnant"
#BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal1 == "CD68")] = "IFnant"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash1Human" & BAL_merge$Animal1 == "RMs8")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal1 == "CD68")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash1Human" & BAL_merge$Animal1 == "RRc17")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal1 == "RGo9")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash1Human" & BAL_merge$Animal1 == "LG92")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal1 == "LC10")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash3Human" & BAL_merge$Animal1 == "RCk17")] = "IFnant"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash4Human" & BAL_merge$Animal1 == "RRj11")] = "IFnant"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash3Human" & BAL_merge$Animal1 == "RNi17")] = "IFnant"

#Aniaml2
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash4Human" & BAL_merge$Animal2 == "KV31")] = "IFnant"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal2 == "CD68")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal2 == "RGo9")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash4Human" & BAL_merge$Animal2 == "RRj11")] = "IFnant"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal2 == "LC10")] = "Control"
BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash4Human" & BAL_merge$Animal2 == "RZs14")] = "IFnant"
#BAL_merge$Treatment[which(BAL_merge$hash.ID=="Hash2Human" & BAL_merge$Animal2 == "CD68")] = "Control"



# BAL_merge$Sample_m = paste0(BAL_merge@meta.data$Sample,
#                              "_", BAL_merge@meta.data$Day,"_",
#                              BAL_merge$Type)

BAL_merged_Filtered <- subset(x = BAL_merge,
                              subset= (nFeature_RNA >= 500) & (nFeature_RNA <=3500) &
                                (nCount_RNA >= 250) & (log10GenesPerUMI > 0.8))

Idents(BAL_merged_Filtered) = "Sample_me"
pdf("/data/runs/Alysis/RAMAS_Data_Alysis/BAL_Alysis/BALs_Merged_VlnPlot.pdf", width =40, height =5)
VlnPlot(object = BAL_merged_Filtered, features = c("nFeature_RNA", "nCount_RNA","percent.hbb","log10GenesPerUMI",
                                                   "percent.mt", "percent.rps",
                                                   "percent.rpl"), ncol = 3)
dev.off()

BAL_merged_Filtered$Animal_Day[which(BAL_merged_Filtered$Animal_Day == "CD68-KV31_D4")] = "RMs8-CD68_D4"

BAL.list_0 <- SplitObject(BAL_merged_Filtered, split.by = "Animal_Day")

for (i in names(BAL.list_0)) {
  BAL.list_0[[i]] <- SCTransform(BAL.list_0[[i]], verbose = FALSE)
}

BAL_sctransformed = BAL.list_0

"Select Integration features"
BAL.features <- SelectIntegrationFeatures(object.list = BAL_sctransformed, nfeatures = 2000)
"Prepare 'SCT' based Integration"
BAL_sctransformed <- PrepSCTIntegration(object.list = BAL_sctransformed,
                                        anchor.features = BAL.features)
"Find Anchors within datasets(samples)
      to integrate using 'SCT'"
BAL.anchors <- FindIntegrationAnchors(object.list = BAL_sctransformed, normalization.method = "SCT",
                                    anchor.features = BAL.features)

BAL.integrated <- IntegrateData(anchorset = BAL.anchors, normalization.method = "SCT")

BAL.integrated <- RunPCA(object = BAL.integrated, verbose = FALSE)
ElbowPlot(BAL.integrated, ndims = 30)
BAL.integrated <- RunUMAP(object = BAL.integrated, dims = 1:30)
BAL.integrated = FindNeighbors(BAL.integrated, reduction = "pca")
BAL.integrated <- FindClusters(object = BAL.integrated, graph.me = "integrated_snn" ,resolution = 0.6)
DimPlot(BAL.integrated, group.by  = "seurat_clusters", split.by = "Treatment")

BAL.integrated$Day_Treatment = paste0(BAL.integrated$Day,"_", BAL.integrated$Treatment)


BAL_Integrated.sce <- as.SingleCellExperiment(BAL.integrated)

"Convert the seurat object to
          singlecellexperiment"

singleR_BAL_intergrated_bp = SingleR(test = BAL_Integrated.sce, ref = ref,
                                     labels = ref$label.main, fine.tune = T,
                                     prune = T)

"Map co-ordites of UMAP onto
      singleR object"


BAL.integrated$singeR_BP = singleR_BAL_intergrated_bp$labels

