library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("HeLa-S3_DNase-seq-DUKE_in_HeLa-S3_cistrome.txt")

HeLa <- dataBrowser(cell_tissue_name = "HeLa-S3")

# TF interactome coupled with CpG methylation and DNase signal
HeLa_intersect <- intersectPeakMatrix(peak_id_x = HeLa$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = HeLa$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(HeLa_intersect, "HeLa-S3_intersectPeakMatrix.rds")




HeLa_intersect <- readRDS("HeLa-S3_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = HeLa_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


