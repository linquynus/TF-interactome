library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("HEK293_DNase-seq-UW_in_HEK293_cistrome.txt")

HEK293 <- dataBrowser(cell_tissue_name = "HEK293")

# TF interactome coupled with CpG methylation and DNase signal
HEK293_intersect <- intersectPeakMatrix(peak_id_x = HEK293$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = HEK293$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(HEK293_intersect, "HEK293_intersectPeakMatrix.rds")




HEK293_intersect <- readRDS("HEK293_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = HEK293_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


