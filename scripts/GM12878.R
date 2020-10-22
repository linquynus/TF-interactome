library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("GM12878_DNase-seq-DUKE_in_GM12878_cistrome.txt")

GM12878 <- dataBrowser(cell_tissue_name = "GM12878")

# TF interactome coupled with CpG methylation and DNase signal
GM12878_intersect <- intersectPeakMatrix(peak_id_x = GM12878$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = GM12878$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(GM12878_intersect, "GM12878_intersectPeakMatrix.rds")




GM12878_intersect <- readRDS("GM12878_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = GM12878_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


