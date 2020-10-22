library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("HEK293T_DNase-seq_in_HEK293T_cistrome.txt")

HEK293T <- dataBrowser(cell_tissue_name = "HEK293T")

# TF interactome coupled with CpG methylation and DNase signal
HEK293T_intersect <- intersectPeakMatrix(peak_id_x = HEK293T$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = HEK293T$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(HEK293T_intersect, "HEK293T_intersectPeakMatrix.rds")




HEK293T_intersect <- readRDS("HEK293T_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = HEK293T_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


