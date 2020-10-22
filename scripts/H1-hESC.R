library(TFregulomeR)

#  DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("H1-hESC_DNase-seq-DUKE_in_H1-hESC_cistrome.txt")

H1 <- dataBrowser(cell_tissue_name = "H1-hESC")

# TF interactome coupled with CpG methylation and DNase signal
H1_intersect <- intersectPeakMatrix(peak_id_x = H1$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = H1$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(H1_intersect, "H1_intersectPeakMatrix.rds")




H1_intersect <- readRDS("H1_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = H1_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


