library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("IMR-90_DNase-seq_in_IMR-90_cistrome.txt")

IMR90 <- dataBrowser(cell_tissue_name = "IMR-90")

# TF interactome coupled with CpG methylation and DNase signal
IMR90_intersect <- intersectPeakMatrix(peak_id_x = IMR90$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = IMR90$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(IMR90_intersect, "IMR90_intersectPeakMatrix.rds")




IMR90_intersect <- readRDS("IMR90_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = IMR90_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


