library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("SK-N-SH_DNase-seq_in_SK-N-SH_cistrome.txt")

SKNSH <- dataBrowser(cell_tissue_name = "SK-N-SH")

# TF interactome coupled with CpG methylation and DNase signal
SKNSH_intersect <- intersectPeakMatrix(peak_id_x = SKNSH$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = SKNSH$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(SKNSH_intersect, "SKNSH_intersectPeakMatrix.rds")




SKNSH_intersect <- readRDS("SKNSH_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = SKNSH_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


