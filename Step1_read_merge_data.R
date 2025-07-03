#Step 1

options(rgl.useNULL=TRUE)
options(java.parameters = "-Xmx28g")

#Install Voltron from github if not done previously
#devtools::install_github("Artur-man/VoltRon")

library(Seurat)
library(VoltRon)
library(rJava)
library(RBioFormats)
library(arrow)

working_dir <- "/path/to/working/directory/"
data_dir <- "/path/to/data/dir/"

setwd(working_dir)

#Read all Xenium output folders as they come from the machine
d0_WT1 <- importXenium(paste0(data_dir, "dist10_d0_WT1/outs/"), resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_WT2 <- importXenium(paste0(data_dir, "dist10_d0_WT2/outs/"), resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_WT3 <- importXenium(paste0(data_dir, "dist10_d0_WT3/outs/"), resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_KO1_slice1 <- importXenium(paste0(data_dir, "dist10_d0_KO1_slice1/outs/"), resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_KO2_slice2 <- importXenium(paste0(data_dir, "dist10_d0_KO2_slice2/outs/"), resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_KO3 <- importXenium(paste0(data_dir, "dist10_d0_KO3/outs/"), resolution_level = 3, overwrite_resolution = TRUE, import_molecules = TRUE)

d0_WT1$Sample <- "d0_WT1"
d0_WT2$Sample <- "d0_WT2"
d0_WT3$Sample <- "d0_WT3"
d0_KO1_slice1$Sample <- "d0_KO1_slice1"
d0_KO2_slice2$Sample <- "d0_KO2_slice2"
d0_KO3$Sample <- "d0_KO3"


#merge day 0 objects
expr <- list(d0_WT1, d0_WT2, d0_WT3, d0_KO1_slice1, d0_KO2_slice2, d0_KO3)
all_d0 <- merge(expr[[1]], expr[-1])

#save for usage in downstream scripts
saveRDS(all_d0, "all_d0_mols.rds")

#Remove objects from memory
rm(list=ls(pattern="d0"))


d15_WT1 <- importXenium(paste0(data_dir, "dist10_d15_WT1/outs/"), resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_WT2 <- importXenium(paste0(data_dir, "dist10_d15_WT2/outs/"), resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_WT3 <- importXenium(paste0(data_dir, "dist10_d15_WT3/outs/"), resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_KO1 <- importXenium(paste0(data_dir, "dist10_d15_KO1/outs/"), resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_KO2 <- importXenium(paste0(data_dir, "dist10_d15_KO2/outs/"), resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_KO3_slice2 <- importXenium(paste0(data_dir, "dist10_d15_KO3_slice2/outs/"), resolution_level = 2, overwrite_resolution = TRUE, import_molecules = TRUE)

d15_WT1$Sample <- "d15_WT1"
d15_WT2$Sample <- "d15_WT2"
d15_WT3$Sample <- "d15_WT3"
d15_KO1$Sample <- "d15_KO1"
d15_KO2$Sample <- "d15_KO2"
d15_KO3_slice2$Sample <- "d15_KO3_slice2"

expr <- list(d15_WT1, d15_WT2, d15_WT3, d15_KO1, d15_KO2, d15_KO3_slice2)
all_d15 <- merge(expr[[1]], expr[-1])

saveRDS(all_d15, "all_d15_mols.rds")

rm(list=ls(pattern="d15"))
