library(TFregulomeR)

# DNase-seq signal in cistrome
# available in TF-interactome website
dnase <- read.table("MCF-7_DNase-seq-DUKE_in_MCF-7_cistrome.txt")

MCF7 <- dataBrowser(cell_tissue_name = "MCF-7")

# TF interactome coupled with CpG methylation and DNase signal
MCF7_intersect <- intersectPeakMatrix(peak_id_x = MCF7$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = MCF7$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = dnase)
saveRDS(MCF7_intersect, "MCF7_intersectPeakMatrix.rds")




MCF7_intersect <- readRDS("MCF7_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = MCF7_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


