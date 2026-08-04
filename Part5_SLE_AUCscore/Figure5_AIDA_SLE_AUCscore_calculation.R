#AUCell analysis in cell types
library(AUCell)
library(Seurat)
library(SingleCellExperiment) 

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

#SLE marker geneSets obtained from https://pubmed.ncbi.nlm.nih.gov/31971918/
##https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
SLE_marker_geneSets <-  list(geneSet1 = c("IFI44L", "EPSTI1", "HERC5", "IFI44", "IFI27", "RSAD2", "ISG15", "OASL", "IFIT3", "CMPK2", "IFIT1", "USP18", "SIGLEC1", "MX2", "MX1", "LY6E", "PLSCR1", "OAS2", "SPATS2L", "OAS1", "OAS3", "PARP12", "SCO2", "IFIH1", "IFITM3", "IFITM1", "DDX60", "TRIM22", "RTP4", "SAMD9L", "XAF1", "IFI35", "TNFAIP6", "MT2A", "LAP3", "HERC6", "FBXO6", "TDRD7", "IFIT5", "PHF11", "IFIT2", "TAP1", "TOR1B", "TNFSF13B", "ELANE", "IRF7", "DHX58", "ZBP1", "TYMP", "GRN", "LAMP3", "IFI6", "NTNG2", "SAT1", "SERPING1", "DEFA4", "DDX58", "SAMD9", "PARP9", "MT1E", "IRF9", "TCN2", "MT1HL1", "HSH2D", "ISG20", "ZCCHC2", "SP100", "LMO2", "GBP1", "IFITM2", "LGALS3BP", "REM2", "TNFSF10", "HESX1", "MT1A", "CCR1", "CEACAM1", "MYD88", "BST2", "LHFPL2", "FFAR2", "MT1F"))

#Other cell types 
#single cell experiment
#Myeloid_pDC
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce <- as.SingleCellExperiment(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered, assay = "RNA")
#NK
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce <- as.SingleCellExperiment(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered, assay = "RNA")
#T
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce <- as.SingleCellExperiment(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered, assay = "RNA")
#B
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce <- as.SingleCellExperiment(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered, assay = "RNA")

#obtain expression matrix and convert to sparse
#Myeloid_pDC
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce_exprMatrix <- assay(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce)
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce_exprMatrix <- as(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce_exprMatrix, "dgCMatrix")
#NK
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce_exprMatrix <- assay(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce)
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce_exprMatrix <- as(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce_exprMatrix, "dgCMatrix")
#T
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce_exprMatrix <- assay(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce)
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce_exprMatrix <- as(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce_exprMatrix, "dgCMatrix")
rm(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered)
#B
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce_exprMatrix <- assay(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce)
AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce_exprMatrix <- as(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce_exprMatrix, "dgCMatrix")

#Calculate enrichment scores from cell rangings
#Ranking cells
Myeloid_pDC_cells_rankings <- AUCell_buildRankings(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce_exprMatrix, plotStats=TRUE)
NK_ILC_cells_rankings <- AUCell_buildRankings(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce_exprMatrix, plotStats=TRUE)
T_cells_rankings <- AUCell_buildRankings(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce_exprMatrix, plotStats=TRUE)
B_Plasma_Cell_cells_rankings <- AUCell_buildRankings(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce_exprMatrix, plotStats=TRUE)

#Calculate enrichment scores
AIDA_myeloid_SLE_AUC_rankings <- AUCell_run(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_Myeloid_pDC_merged_RPCA_joined_clustered_sce_exprMatrix, SLE_marker_geneSets, aucMaxRank=nrow(cells_rankings)*0.05)
AIDA_NK_ILC_SLE_AUC_rankings <- AUCell_run(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_NK_ILC_QC_RPCA_joined_clustered_sce_exprMatrix, SLE_marker_geneSets, aucMaxRank=nrow(NK_ILC_cells_rankings)*0.05)
AIDA_T_SLE_AUC_rankings <- AUCell_run(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_T_QC_RPCA_joined_clustered_sce_exprMatrix, SLE_marker_geneSets, aucMaxRank=nrow(T_cells_rankings)*0.05)
AIDA_B_Plasma_Cell_SLE_AUC_rankings <- AUCell_run(AIDA_Phase1_DataFreeze_v2_Step03_Seurat_JPKRSGTH_B_Plasma_Cell_QC_RPCA_joined_clustered_sce_exprMatrix, SLE_marker_geneSets, aucMaxRank=nrow(B_Plasma_Cell_cells_rankings)*0.05)

#Set the assignment thresholds
set.seed(123)
myeloid_cells_assignment_rankings <- AUCell_exploreThresholds(AIDA_myeloid_SLE_AUC_rankings, plotHist=TRUE, nCores=1, assign=TRUE)
NK_ILC_cells_assignment_rankings <- AUCell_exploreThresholds(AIDA_NK_ILC_SLE_AUC_rankings, plotHist=TRUE, nCores=1, assign=TRUE)
T_cells_assignment_rankings <- AUCell_exploreThresholds(AIDA_T_SLE_AUC_rankings, plotHist=TRUE, nCores=1, assign=TRUE)
B_Plasma_Cell_cells_assignment_rankings <- AUCell_exploreThresholds(AIDA_B_Plasma_Cell_SLE_AUC_rankings, plotHist=TRUE, nCores=1, assign=TRUE)

#Save outputs
saveRDS(AIDA_myeloid_SLE_AUC_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_myeloid_SLE_AUC_rankings.rds")
saveRDS(myeloid_cells_assignment_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_myeloid_SLE_AUC_rankings_cells_assignment.rds")
saveRDS(AIDA_NK_ILC_SLE_AUC_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_NK_ILC_SLE_AUC_rankings.rds")
saveRDS(NK_ILC_cells_assignment_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_NK_ILC_SLE_AUC_rankings_cells_assignment.rds")
saveRDS(AIDA_T_SLE_AUC_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_T_SLE_AUC_rankings.rds")
saveRDS(T_cells_assignment_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_T_SLE_AUC_rankings_cells_assignment.rds")
saveRDS(AIDA_B_Plasma_Cell_SLE_AUC_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_B_Plasma_Cell_SLE_AUC_rankings.rds")
saveRDS(B_Plasma_Cell_cells_assignment_rankings, "/mnt/icbs_shared_storage_general_2022/Benjamaporn/Analysis2026/Analysis202607/SLE/AUCell/AIDA_B_Plasma_Cell_SLE_AUC_rankings_cells_assignment.rds")
