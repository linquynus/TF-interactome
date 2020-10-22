library(TFregulomeR)

# ATAC-seq signal in cistrome
# available in TF-interactome website
atac <- read.table("HCT116_ATACseq_in_HCT116_cistrome.txt")

HCT116 <- dataBrowser(cell_tissue_name = "HCT116")

# TF interactome coupled with CpG methylation and DNase signal
HCT116_intersect <- intersectPeakMatrix(peak_id_x = HCT116$ID, 
                                      motif_only_for_id_x = TRUE,
                                      peak_id_y = HCT116$ID, 
                                      motif_only_for_id_y = TRUE, 
                                      methylation_profile_in_narrow_region = TRUE, 
                                      external_source = atac)
saveRDS(HCT116_intersect, "HCT116_intersectPeakMatrix.rds")




HCT116_intersect <- readRDS("HCT116_intersectPeakMatrix.rds")

interactome3D(intersectPeakMatrix = HCT116_intersect,
              return_interactome_with_mCpG = TRUE,
              return_interactome_with_external_source = TRUE)


