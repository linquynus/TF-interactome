library(TFregulomeR)

# ATAC-seq signal in cistrome
# available in TF-interactome website
atac <- read.table("K562_ATACseq_in_K562_cistrome.txt")

K562 <- dataBrowser(cell_tissue_name = "K562")

# TF interactome coupled with CpG methylation and DNase signal
K562_intersect <- intersectPeakMatrix(peak_id_x = K562$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = K562$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = atac)
saveRDS(K562_intersect, "K562_intersectPeakMatrix.rds")



K562_intersect <- readRDS("K562_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = K562_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


