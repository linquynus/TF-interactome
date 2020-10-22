library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("HepG2_DNase-seq-UW_in_HepG2_cistrome.txt")

HepG2 <- dataBrowser(cell_tissue_name = "HepG2")

# TF interactome coupled with CpG methylation and DNase signal
HepG2_intersect <- intersectPeakMatrix(peak_id_x = HepG2$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = HepG2$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(HepG2_intersect, "HepG2_intersectPeakMatrix.rds")




HepG2_intersect <- readRDS("HepG2-S3_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = HepG2_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


