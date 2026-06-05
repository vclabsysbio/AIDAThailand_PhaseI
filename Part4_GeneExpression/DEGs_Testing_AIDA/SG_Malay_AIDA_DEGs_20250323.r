library(limma)
library(edgeR)
library(variancePartition)
library(Seurat)
library(dplyr)

#Output path
AIDA_DEGs_Output_PATH <- "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202603/Part4_AIDA_Age_Sex_DEGs_Output/SG_AIDA/"

#Myeloid + pDC
Myeloid_pDC_Output_PATH <- "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2025/Analysis202510/myeloid_pDC_integration"
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered <- readRDS(file.path(Myeloid_pDC_Output_PATH, "/AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered.rds"))
#NK
NK_ILC_Output_PATH <- "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2025/Analysis202511/NK_ILC_integration"
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered <- readRDS(file.path(NK_ILC_Output_PATH, "/AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered.rds"))
#T
T_Output_PATH <- "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2025/Analysis202511/T_integration"
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered <- readRDS(file.path(T_Output_PATH, "/AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered.rds"))
#B and Plasma Cells
B_PlasmaCell_Output_PATH <- "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2025/Analysis202511/B_PlasmaCell_integration"
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered <- readRDS(file.path(B_PlasmaCell_Output_PATH, "/AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered.rds"))

#Myeloid_pDC cells in SG 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_SG_Malay <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered, subset = Ethnicity == "Malay")

#Author Annotation (Prepare count data)
#Using limma and voom for linear mixed model
#Only same cell sets with proportion analysis
Myeloid_pDC_Author_Annotation_celltype <- c("Myeloid_unknown", "CD14+_Monocyte", "CD16+_Monocyte", "cDC2")

#Metadata
#perLibrary
SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary <- select(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_SG_Malay@meta.data, Library, Genotyping_ID, DCP_ID, Age, Sex, Author_Annotation, Blue_ref, Red_ref, Orange_ref) %>%
distinct(DCP_ID, Library, .keep_all = TRUE)  
row.names(SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary) <- NULL

#Format metadata
SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary$DCP_ID <- gsub("_", "-", SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary$DCP_ID)
SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary$Library <-gsub("_", "-", SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary$Library)
rownames(SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary) <- paste(SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary$Library, SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary$DCP_ID, sep = "_")

#LMM by dream 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_SG_Malay_Author_Annotation <- list()
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP <- list()
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_countdata <- list()
SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary_Author_Annotation <- list()
Myeloid_pDC_gene_notlowlyExpressedGenes <- list()
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE <- list()
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE_TMM <- list()
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_design <- list()
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE_vobjGenes_blue <- list()
SG_Malay_AIDA_Myeloid_pDC_dream_modelFit <- list()
SG_Malay_AIDA_Myeloid_pDC_dream_modelFit_eBayes <- list()

for (i in Myeloid_pDC_Author_Annotation_celltype) {
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]] <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_SG_Malay, subset = Author_Annotation == i)
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP[[i]] <- AggregateExpression(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]], group.by = c("Library", "DCP_ID"), return.seurat = F)
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_countdata[[i]] <- as.matrix(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP[[i]][["RNA"]])
SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary_Author_Annotation[[i]] <-  SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary[colnames(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_countdata[[i]]), , drop=FALSE]
Myeloid_pDC_gene_notlowlyExpressedGenes[[i]] <- rowSums(edgeR::cpm(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_countdata[[i]]) > 0.5) >= 10
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE[[i]] <- DGEList(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_countdata[[i]][Myeloid_pDC_gene_notlowlyExpressedGenes[[i]], ])
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE_TMM[[i]] <- calcNormFactors(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE[[i]])
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_design[[i]] <- ~ Sex + Blue_ref + Age + (1 | Library)
bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]] <- voomWithDreamWeights(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE_TMM[[i]], bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_Myeloid_pDC_dream_modelFit[[i]] <- dream(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]], bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_Myeloid_pDC_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_Myeloid_pDC_dream_modelFit_eBayes[[i]] <- eBayes(SG_Malay_AIDA_Myeloid_pDC_dream_modelFit[[i]])
}
saveRDS(bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE, file =  file.path(AIDA_DEGs_Output_PATH, "bulk_SG_Malay_AIDA_Myeloid_pDC_integrated_passQC_Library_DCP_DGE.rds"))
saveRDS(SG_Malay_AIDA_Myeloid_pDC_dream_modelFit_eBayes, file =  file.path(AIDA_DEGs_Output_PATH, "SG_Malay_AIDA_Myeloid_pDC_dream_modelFit_eBayes.rds"))

#NK cells in SG_Malay 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_SG_Malay <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered, subset = Ethnicity == "Malay")

#Author Annotation (Prepare count data)
#Using limma and voom for linear mixed model
#Only same cell sets with proportion analysis
NK_Author_Annotation_celltype <- c("NK_unknown", "CD16+_NK", "CD56+_NK")

#Metadata
#perLibrary
SG_Malay_AIDA_NK_integrated_metadata_perLibrary <- select(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_SG_Malay@meta.data, Library, Genotyping_ID, DCP_ID, Age, Sex, Author_Annotation, Blue_ref, Red_ref, Orange_ref) %>%
distinct(DCP_ID, Library, .keep_all = TRUE)  
row.names(SG_Malay_AIDA_NK_integrated_metadata_perLibrary) <- NULL

#Format metadata
SG_Malay_AIDA_NK_integrated_metadata_perLibrary$DCP_ID <- gsub("_", "-", SG_Malay_AIDA_NK_integrated_metadata_perLibrary$DCP_ID)
SG_Malay_AIDA_NK_integrated_metadata_perLibrary$Library <-gsub("_", "-", SG_Malay_AIDA_NK_integrated_metadata_perLibrary$Library)
rownames(SG_Malay_AIDA_NK_integrated_metadata_perLibrary) <- paste(SG_Malay_AIDA_NK_integrated_metadata_perLibrary$Library, SG_Malay_AIDA_NK_integrated_metadata_perLibrary$DCP_ID, sep = "_")

#LMM by dream 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation <- list()
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP <- list()
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_countdata <- list()
SG_Malay_AIDA_NK_integrated_metadata_perLibrary_Author_Annotation <- list()
NK_gene_notlowlyExpressedGenes <- list()
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE <- list()
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE_TMM <- list()
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_design <- list()
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE_vobjGenes_blue <- list()
SG_Malay_AIDA_NK_dream_modelFit <- list()
SG_Malay_AIDA_NK_dream_modelFit_eBayes <- list()

for (i in NK_Author_Annotation_celltype) {
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]] <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_SG_Malay, subset = Author_Annotation == i)
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP[[i]] <- AggregateExpression(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]], group.by = c("Library", "DCP_ID"), return.seurat = F)
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_countdata[[i]] <- as.matrix(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP[[i]][["RNA"]])
SG_Malay_AIDA_NK_integrated_metadata_perLibrary_Author_Annotation[[i]] <-  SG_Malay_AIDA_NK_integrated_metadata_perLibrary[colnames(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_countdata[[i]]), , drop=FALSE]
NK_gene_notlowlyExpressedGenes[[i]] <- rowSums(edgeR::cpm(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_countdata[[i]]) > 0.5) >= 10
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE[[i]] <- DGEList(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_countdata[[i]][NK_gene_notlowlyExpressedGenes[[i]], ])
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE_TMM[[i]] <- calcNormFactors(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE[[i]])
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_design[[i]] <- ~ Sex + Blue_ref + Age + (1 | Library)
bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]] <- voomWithDreamWeights(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE_TMM[[i]], bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_NK_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_NK_dream_modelFit[[i]] <- dream(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]], bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_NK_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_NK_dream_modelFit_eBayes[[i]] <- eBayes(SG_Malay_AIDA_NK_dream_modelFit[[i]])
}
saveRDS(bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE, file =  file.path(AIDA_DEGs_Output_PATH, "bulk_SG_Malay_AIDA_NK_integrated_passQC_Library_DCP_DGE.rds"))
saveRDS(SG_Malay_AIDA_NK_dream_modelFit_eBayes, file =  file.path(AIDA_DEGs_Output_PATH, "SG_Malay_AIDA_NK_dream_modelFit_eBayes.rds"))

#T cells in SG_Malay 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_SG_Malay <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered, subset = Ethnicity == "Malay")

#Author Annotation (Prepare count data)
#Using limma and voom for linear mixed model
#Only same cell sets with proportion analysis
T_Author_Annotation_celltype <- c("T_unknown", "CD4+_T_unknown", "CD4+_T_naive", "CD4+_T_cm",  "CD4+_T_cyt", "CD4+_T_em", "CD8+_T_unknown", "CD8+_T_naive", "CD8+_T_GZMBhi", "CD8+_T_GZMKhi", "MAIT",   "gdT_GZMBhi", "gdT_GZMKhi",  "Treg")

#Metadata
#perLibrary
SG_Malay_AIDA_T_integrated_metadata_perLibrary <- select(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_SG_Malay@meta.data, Library, Genotyping_ID, DCP_ID, Age, Sex, Author_Annotation, Blue_ref, Red_ref, Orange_ref) %>%
distinct(DCP_ID, Library, .keep_all = TRUE)  
row.names(SG_Malay_AIDA_T_integrated_metadata_perLibrary) <- NULL

#Format metadata
SG_Malay_AIDA_T_integrated_metadata_perLibrary$DCP_ID <- gsub("_", "-", SG_Malay_AIDA_T_integrated_metadata_perLibrary$DCP_ID)
SG_Malay_AIDA_T_integrated_metadata_perLibrary$Library <-gsub("_", "-", SG_Malay_AIDA_T_integrated_metadata_perLibrary$Library)
rownames(SG_Malay_AIDA_T_integrated_metadata_perLibrary) <- paste(SG_Malay_AIDA_T_integrated_metadata_perLibrary$Library, SG_Malay_AIDA_T_integrated_metadata_perLibrary$DCP_ID, sep = "_")

#LMM by dream 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation <- list()
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP <- list()
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_countdata <- list()
SG_Malay_AIDA_T_integrated_metadata_perLibrary_Author_Annotation <- list()
T_gene_notlowlyExpressedGenes <- list()
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE <- list()
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE_TMM <- list()
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_design <- list()
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE_vobjGenes_blue <- list()
SG_Malay_AIDA_T_dream_modelFit <- list()
SG_Malay_AIDA_T_dream_modelFit_eBayes <- list()

for (i in T_Author_Annotation_celltype) {
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]] <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_SG_Malay, subset = Author_Annotation == i)
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP[[i]] <- AggregateExpression(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]], group.by = c("Library", "DCP_ID"), return.seurat = F)
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_countdata[[i]] <- as.matrix(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP[[i]][["RNA"]])
SG_Malay_AIDA_T_integrated_metadata_perLibrary_Author_Annotation[[i]] <-  SG_Malay_AIDA_T_integrated_metadata_perLibrary[colnames(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_countdata[[i]]), , drop=FALSE]
T_gene_notlowlyExpressedGenes[[i]] <- rowSums(edgeR::cpm(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_countdata[[i]]) > 0.5) >= 10
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE[[i]] <- DGEList(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_countdata[[i]][T_gene_notlowlyExpressedGenes[[i]], ])
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE_TMM[[i]] <- calcNormFactors(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE[[i]])
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_design[[i]] <- ~ Sex + Blue_ref + Age + (1 | Library)
bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]] <- voomWithDreamWeights(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE_TMM[[i]], bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_T_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_T_dream_modelFit[[i]] <- dream(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]], bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_T_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_T_dream_modelFit_eBayes[[i]] <- eBayes(SG_Malay_AIDA_T_dream_modelFit[[i]])
}
saveRDS(bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE, file =  file.path(AIDA_DEGs_Output_PATH, "bulk_SG_Malay_AIDA_T_integrated_passQC_Library_DCP_DGE.rds"))
saveRDS(SG_Malay_AIDA_T_dream_modelFit_eBayes, file =  file.path(AIDA_DEGs_Output_PATH, "SG_Malay_AIDA_T_dream_modelFit_eBayes.rds"))

#B cells in SG_Malay 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_QC_RPCA_joined_clustered_SG_Malay <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered, subset = Ethnicity == "Malay")

#Author Annotation (Prepare count data)
#Using limma and voom for linear mixed model
#Only same cell sets with proportion analysis
B_Author_Annotation_celltype <- c("B_unknown", "naive_B", "memory_B_IGHMlo", "memory_B_IGHMhi")

#Metadata
#perLibrary
SG_Malay_AIDA_B_integrated_metadata_perLibrary <- select(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_QC_RPCA_joined_clustered_SG_Malay@meta.data, Library, Genotyping_ID, DCP_ID, Age, Sex, Author_Annotation, Blue_ref, Red_ref, Orange_ref) %>%
distinct(DCP_ID, Library, .keep_all = TRUE)  
row.names(SG_Malay_AIDA_B_integrated_metadata_perLibrary) <- NULL

#Format metadata
SG_Malay_AIDA_B_integrated_metadata_perLibrary$DCP_ID <- gsub("_", "-", SG_Malay_AIDA_B_integrated_metadata_perLibrary$DCP_ID)
SG_Malay_AIDA_B_integrated_metadata_perLibrary$Library <-gsub("_", "-", SG_Malay_AIDA_B_integrated_metadata_perLibrary$Library)
rownames(SG_Malay_AIDA_B_integrated_metadata_perLibrary) <- paste(SG_Malay_AIDA_B_integrated_metadata_perLibrary$Library, SG_Malay_AIDA_B_integrated_metadata_perLibrary$DCP_ID, sep = "_")

#LMM by dream 
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation <- list()
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP <- list()
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_countdata <- list()
SG_Malay_AIDA_B_integrated_metadata_perLibrary_Author_Annotation <- list()
B_gene_notlowlyExpressedGenes <- list()
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE <- list()
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE_TMM <- list()
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_design <- list()
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE_vobjGenes_blue <- list()
SG_Malay_AIDA_B_dream_modelFit <- list()
SG_Malay_AIDA_B_dream_modelFit_eBayes <- list()

for (i in B_Author_Annotation_celltype) {
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]] <- subset(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_QC_RPCA_joined_clustered_SG_Malay, subset = Author_Annotation == i)
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP[[i]] <- AggregateExpression(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_QC_RPCA_joined_clustered_SG_Malay_Author_Annotation[[i]], group.by = c("Library", "DCP_ID"), return.seurat = F)
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_countdata[[i]] <- as.matrix(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP[[i]][["RNA"]])
SG_Malay_AIDA_B_integrated_metadata_perLibrary_Author_Annotation[[i]] <-  SG_Malay_AIDA_B_integrated_metadata_perLibrary[colnames(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_countdata[[i]]), , drop=FALSE]
B_gene_notlowlyExpressedGenes[[i]] <- rowSums(edgeR::cpm(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_countdata[[i]]) > 0.5) >= 10
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE[[i]] <- DGEList(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_countdata[[i]][B_gene_notlowlyExpressedGenes[[i]], ])
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE_TMM[[i]] <- calcNormFactors(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE[[i]])
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_design[[i]] <- ~ Sex + Blue_ref + Age + (1 | Library)
bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]] <- voomWithDreamWeights(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE_TMM[[i]], bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_B_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_B_dream_modelFit[[i]] <- dream(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE_vobjGenes_blue[[i]], bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_design[[i]], SG_Malay_AIDA_B_integrated_metadata_perLibrary_Author_Annotation[[i]])
SG_Malay_AIDA_B_dream_modelFit_eBayes[[i]] <- eBayes(SG_Malay_AIDA_B_dream_modelFit[[i]])
}
saveRDS(bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE, file =  file.path(AIDA_DEGs_Output_PATH, "bulk_SG_Malay_AIDA_B_integrated_passQC_Library_DCP_DGE.rds"))
saveRDS(SG_Malay_AIDA_B_dream_modelFit_eBayes, file =  file.path(AIDA_DEGs_Output_PATH, "SG_Malay_AIDA_B_dream_modelFit_eBayes.rds"))

