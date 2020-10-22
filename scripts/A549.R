library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("A549_DNase-seq-UW_in_A549_cistrome.txt")

A549 <- dataBrowser(cell_tissue_name = "A549")

# TF interactome coupled with CpG methylation and DNase signal
A549_intersect <- intersectPeakMatrix(peak_id_x = A549$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = A549$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(A549_intersect, "A549_intersectPeakMatrix.rds")


A549_intersect <- readRDS("A549_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = A549_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


