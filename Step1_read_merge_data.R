#Step1_read_merge_data.R

options(rgl.useNULL=TRUE)

#Install Voltron from github
devtools::install_github("Artur-man/VoltRon")

library(Seurat)
library(VoltRon)


setwd("/path/to/working/directory/")

#Read all Xenium output folders as they come from the machine
d0_WT1 <- importXenium("d0_WT1/outs/", resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_WT2 <- importXenium("d0_WT2/outs/", resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_WT3 <- importXenium("d0_WT3/outs/", resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_KO1_slice1 <- importXenium("d0_KO1_slice1/outs/", resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_KO2_slice2 <- importXenium("d0_KO2_slice2/outs/", resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_KO3 <- importXenium("d0_KO3/outs/", resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_WT1 <- importXenium("d15_WT1/outs/", resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_WT2 <- importXenium("d15_WT2/outs/", resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_WT3 <- importXenium("d15_WT3/outs/", resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_KO1 <- importXenium("d15_KO1/outs/", resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_KO2 <- importXenium("d15_KO2/outs/", resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_KO3_slice2 <- importXenium("d15_KO3_slice2/outs/", resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

#merge day 0 and d15 objects
expr <- list(d0_WT1, d0_WT2, d0_WT3, d0_KO1_slice1, d0_KO2_slice2, d0_KO3)
all_d0 <- merge(expr[[1]], expr[-1])

expr <- list(d15_WT1, d15_WT2, d15_WT3, d15_KO1, d15_KO2, d15_KO3_slice2)
all_d15 <- merge(expr[[1]], expr[-1])

#save for usage in downstream scripts
saveRDS(all_d15, "all_d15_mols.rds")
saveRDS(all_d0, "all_d0_mols.rds")
