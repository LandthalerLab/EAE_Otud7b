#Step 4

options(rgl.useNULL=TRUE)

#Install Voltron from github
devtools::install_github("Artur-man/VoltRon")

library(Seurat)
library(VoltRon)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(tibble)
library(igraph)

#From https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


setwd("/path/to/working/directory/")


all_d0 <- readRDS("all_d0_mols_celltypeann_reg.rds")
all_d15 <- readRDS("all_d15_mols_celltypeann_reg.rds")


#Annotate regions: first dorsal/ventral at both d0 and d15, then lesions at d15

#Dorsal/ventral annotations  at day 0 samples are done manually. The region annotations used in the paper can be transferred from the object GSE286422_all_d15_mols_reg_ann_lesions_v4.rds
#downloadable at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286422

#Adjust alpha: if 0, do not show cells. Increase up to 1 to increase opacity of cells


#Annotate day 0 sections
all_d0 <- annotateSpatialData(all_d0, label = "annotation_d0_WT1", use.image = FALSE, assay = "Assay1",
                               group.by = "celltype", pt.size = 1,
                               image_name = "main", channel = "H&E", alpha=0.2)
all_d0 <- annotateSpatialData(all_d0, label = "annotation_d0_WT2", use.image = FALSE, assay = "Assay3",
                              group.by = "celltype", pt.size = 1,
                              image_name = "main", channel = "H&E", alpha=0.2)
all_d0 <- annotateSpatialData(all_d0, label = "annotation_d0_WT3", use.image = FALSE, assay = "Assay5",
                              group.by = "celltype", pt.size = 1,
                              image_name = "main", channel = "H&E", alpha=0.2)
all_d0 <- annotateSpatialData(all_d0, label = "annotation_d0_KO1_slice1", use.image = FALSE, assay = "Assay7",
                              group.by = "celltype", pt.size = 1,
                              image_name = "main", channel = "H&E", alpha=0.2)
all_d0 <- annotateSpatialData(all_d0, label = "annotation_d0_KO2_slice2", use.image = FALSE, assay = "Assay9",
                              group.by = "celltype", pt.size = 1,
                              image_name = "main", channel = "H&E", alpha=0.2)
all_d0 <- annotateSpatialData(all_d0, label = "annotation_d0_KO3", use.image = FALSE, assay = "Assay11",
                              group.by = "celltype", pt.size = 1,
                              image_name = "main", channel = "H&E", alpha=0.2)

for (the_sample in unique(all_d0@metadata@cell$Sample)) {
  the_assay=subset(all_d0@metadata@cell, Sample==the_sample) %>% rownames() %>% gsub(".*\\_([^\\_]*)$", "\\1", ., perl=TRUE) %>% unique()
  vrSpatialPlot(all_d0, group.by = paste0("annotation_", the_sample), alpha = 0.3, plot.segments = TRUE, background = c("main", "H&E"), assay=the_assay)
  ggsave(paste0("annotation_", the_sample, "_HE.pdf"))
  vrSpatialPlot(all_d0, group.by = paste0("annotation_", the_sample), alpha = 1, plot.segments = TRUE, background.color = "black", assay=the_assay)
  ggsave(paste0("annotation_", the_sample, ".pdf"))
}


#For day 15, we start out with defining lesions based on T cells. This is done as following:
#First, groups of T cells are defined, however not based on the celltype definition from above but based on whether they have a least three counts of a typical T cell marker gene
#In order to form a group, at least 5 cells have to be within 100 µm distances (variables intra_lesion_distance and min_group_size). T cells outsize of such groups are in small cluster
#and not taken into account.
#Second, the T cell groups are spatially expanded, i.e. proximal and distal regions are defined. At the end of this script, plots are generated with cells colored by lesion/region.
#However, these are not the lesions/regions that are used in the paper! We have subsequently manually annotated lesions based on visual examination of the H&E staining.


#T cell genes
lesion_genes <- c("Cd3e", "Cd8a", "Cd4")

#Dataframes with expressions of lesion genes
counts_lesion_genes <- vrData(all_d15) %>% as.data.frame() %>% rownames_to_column(var="gene") %>% filter(gene %in% lesion_genes) %>%
  column_to_rownames(var="gene") %>% t() %>% as.data.frame()
norm_lesion_genes <- vrData(all_d15, norm=TRUE) %>% as.data.frame() %>% rownames_to_column(var="gene") %>% filter(gene %in% lesion_genes) %>%
  column_to_rownames(var="gene") %>% t() %>% as.data.frame()
df_norm <- cbind.data.frame(all_d15@metadata@cell, vrCoordinates(all_d15), norm_lesion_genes)
df_counts <- cbind.data.frame(all_d15@metadata@cell, vrCoordinates(all_d15), counts_lesion_genes)
df2 <- df_norm %>% dplyr::select(all_of(lesion_genes))
df2 <- reshape2::melt(df2)
df_norm$cell <- rownames(df_norm)
df_counts$cell <- rownames(df_counts)

#In this loop, we identify the lesion-defining cells, and group them spatially.
lesion_t <- 2
#lesion_distance is 100 for T cells
intra_lesion_distance <- 100
min_group_size <- 5
expr <- list()

for (the_sample in unique(all_d15@metadata@cell$Sample)) {
  df_counts_temp <- df_norm %>% filter(Sample == the_sample)
  metadata_temp <- all_d15@metadata@cell %>% filter(Sample == the_sample)
  metadata_temp$cell <- rownames(metadata_temp)
  #Get coordinates of lesion-defining cells, i.e. which have at least lesion_t+1 of the lesion defining genes (e.g. Cd3/Cd4/Cd8 )
  points_coords <- df_counts_temp %>%
    filter(if_any(contains(lesion_genes), ~ . > lesion_t)) %>% dplyr::select(c("x", "y", "cell"))
  #Calculate adjacency matrix, everything that is closer of <intra_lesion_distance> µm to each other is put into a group
  adj2 <- matrix(as.numeric(as.matrix(dist(df_counts_temp %>% filter(if_any(contains(lesion_genes), ~ . > 0)) %>%
                                             dplyr::select(c("x", "y"))))) < intra_lesion_distance, nrow = nrow(df_counts_temp %>% filter(if_any(contains(lesion_genes), ~ . > 0))))
  g <- graph_from_adjacency_matrix(adj2)
  points_coords$lesiongroup <- components(g)$membership
  metadata_temp <- left_join(metadata_temp, points_coords %>%
                               dplyr::select(c("cell", "lesiongroup")), by="cell") %>% replace(., is.na(.), 0)
  rownames(metadata_temp) <- metadata_temp$cell
  #Groups with less than five cells are filtered out
  small_clusts <- metadata_temp %>% group_by(lesiongroup) %>% tally() %>% filter(n<min_group_size) %>% pull(lesiongroup)
  metadata_temp <- metadata_temp %>% mutate(lesiongroup=if_else(lesiongroup %in% small_clusts, 10000, lesiongroup))
  metadata_temp$lesiongroup <- paste0(metadata_temp$lesiongroup, "_", the_sample)
  expr[[the_sample]] <- metadata_temp
}
df <- do.call(rbind.data.frame, expr)

#Add lesiongroup information to metadata frame
all_d15@metadata@cell <- df

#Plot lesiongroups, i.e. groups of T cells that define lesions

assay_names <- as.list(as.vector(all_d15@sample.metadata %>% filter(Assay=="Xenium") %>% rownames()))
names(assay_names) <- all_d15@sample.metadata %>% filter(Assay=="Xenium") %>% pull(var="Sample")

#Plots for individual sections
for (the_sample in unique(all_d15@sample.metadata$Sample)) {
  df <- subset(all_d15@metadata@cell, Sample==the_sample)
  dist_names <- df %>% pull(var="lesiongroup") %>% unique()
  the_cols <- gg_color_hue(length(dist_names))
  names(the_cols) <- dist_names
  the_cols[grep("^0_", names(the_cols), value=TRUE, perl=TRUE)] <- "gray95"
  the_cols[grep("^10000_", names(the_cols), value=TRUE, perl=TRUE)] <- "gray50"
  vrSpatialPlot(all_d15, assay=assay_names[[the_sample]], colors=the_cols, group.by = "lesiongroup", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")
  ggsave(paste(the_sample, "colored_by_lesiongroup", "pdf", sep="."), useDingbats=FALSE, height=6, width=8)
}

#Now, we expand the groups of T cells that define a lesion (lesiongroups)
expr1 <- list()
for (the_sample in unique(all_d15@metadata@cell$Sample)) {
  df_counts_temp <- df_counts %>% filter(Sample == the_sample)
  metadata_temp <- all_d15@metadata@cell %>% filter(Sample == the_sample)
  #Get cells that are in lesion cells, i.e. lesion-defining cells that are however not in too small clusters (those starting with 10000)
  only_lesiongroup <- metadata_temp %>% dplyr::filter(!(grepl("^0_|^10000_", lesiongroup, perl=TRUE, ))) %>%
    dplyr::select(all_of(c("cell", "lesiongroup")))
  
  #Get distance matrix for the entire sample, i.e. distances of every cell to every other cells
  d <- dist(df_counts_temp %>% dplyr::select(c("x", "y")))
  distance_table <- as.data.frame(as.matrix(d))
  #For every cell that is not a lesion-defining cell (i.e. starting with 0), calculate distance to the closest lesion-defining cell in a lesion group
  expr2 <- list()
  for (the_cell in metadata_temp %>% dplyr::filter(grepl("^0_", lesiongroup, perl=TRUE, )) %>% pull(var="cell")) {
    df <- distance_table %>% dplyr::select(all_of(the_cell)) %>% rownames_to_column(var="cell")
    colnames(df) <- c("cell", "distance")
    expr2[[the_cell]] <- merge(df, only_lesiongroup, by="cell") %>% slice_min(order_by=distance, n=1)
  }
  distances <- do.call(rbind,expr2) %>% replace(., is.na(.), 0)
  #Define distance group (measured in micrometers)
  distances_2 <- distances %>% mutate(distance = case_when(distance < 50 ~ "proximal",
                                                           distance >= 50 & distance < 100 ~ "intermediate",
                                                           distance >= 100 & distance < 200 ~ "distal",
                                                           .default = "outside"))
  distances_2$cell <- rownames(distances_2)
  df4 <- metadata_temp %>% dplyr::select(c("cell", "celltype"))
  #distances_3 contains for every non-lesion defining cell the distance group and to which lesion group it belongs
  distances_3 <- left_join(distances_2, df4, by="cell")
  distances_3$distgroup <- distances_3$distance
  metadata_temp$distance <- NULL
  #In the column distance, the cells are classified as either "smalllesiondefininggroup" (lesion defining cells but in too small grousps),
  #or the distance group (which is in column distgroup) and the name of the lesiongroup
  expr1[[the_sample]] <- left_join(metadata_temp, distances_3 %>% mutate(distance=paste(distance, lesiongroup, sep="_")) %>% dplyr::select(c("cell", "distance", "distgroup")), by="cell") %>%
    mutate(distance=ifelse(grepl("^10000_", lesiongroup), "smalllesiondefininggroup", ifelse(grepl("^0_", lesiongroup), distance, paste0("lesiondefining_", lesiongroup)))) %>%
    mutate(distgroup=ifelse(grepl("^10000_", lesiongroup), "smalllesiondefininggroup", ifelse(grepl("^0_", lesiongroup), distgroup, "lesiondefining")))
  rownames(expr1[[the_sample]]) <- expr1[[the_sample]]$cell
}
df <- do.call(rbind.data.frame, expr1)
all_d15@metadata@cell <- df

#Plot lesions and the expanded segments
for (the_sample in unique(all_d15@sample.metadata$Sample)) {
  df <- subset(all_d15@metadata@cell, Sample==the_sample) %>% mutate(lesioncluster=gsub("^[^\\_]*_(.*)", "\\1", distance)) %>% arrange(distgroup, lesioncluster)
  dist_names <- df %>% pull(var="distance") %>% unique()
  the_cols <- gg_color_hue(length(dist_names))
  names(the_cols) <- dist_names
  the_cols[grep("^lesiondefining_", names(the_cols), value=TRUE, perl=TRUE)] <- "#fbff02"
  the_cols["smalllesiondefininggroup"] <- "gray50"
  vrSpatialPlot(all_d15, assay=assay_names[[the_sample]], colors=the_cols, group.by = "distance", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")
  ggsave(paste(the_sample, "by_lesiongroup", "pdf", sep="."), useDingbats=FALSE, height=6, width=24)
  vrSpatialPlot(all_d15, assay=assay_names[[the_sample]], group.by = "distgroup", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background.color = "black")
  ggsave(paste(the_sample, "by_distancegroup", "pdf", sep="."), useDingbats=FALSE)
}

#At this stage, lesions were manually annotated. In the final object (downloadable from GEO, see above), lesions ended up in the metadata column Tlesiongroup_ann.
#This then also included the region (dorsal or ventral). 
#If you would want to recapitulate the analysis we did in our paper, simply transfer the medata column Tlesiongroup_ann from our object.
#Otherwise, create your own Tlesiongroup_ann column based applying custom criteria to the lesiongroup metadata column within the all_d15 VoltRon object.

#Read in objected downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286422 and transfer information
d15_final_2 <- readRDS("GSE286422_all_d15_mols_reg_ann_lesions_v4.rds")
all_d15@metadata@cell$Tlesiongroup_ann <- d15_final_2@metadata@cell$Tlesiongroup_ann

#Expand T cell lesions from the Tlesiongroup_ann column
all_d15@metadata@cell <- cbind.data.frame(all_d15@metadata@cell, vrCoordinates(all_d15))
expr1 <- list() 
for (the_sample in unique(all_d15@metadata@cell$Sample)) {
  #Create a temporary dataframe to work on
  metadata_df <- all_d15@metadata@cell %>% filter(Sample == the_sample)
  #Get the cells from which distances are to be calculated (basecells) 
  metadata_df_basecells <- metadata_df %>% filter(grepl("dorsal|ventral", Tlesiongroup_ann, perl=TRUE))
  metadata_df_basecells$cell <- rownames(metadata_df_basecells)
  #Calculate distance matrix
  d <- dist(metadata_df %>% dplyr::select(c("x", "y")))
  distance_table <- as.data.frame(as.matrix(d))
  #For every cell, calculate distance to the closest basecell
  expr2 <- list()
  for (the_cell in rownames(metadata_df)) {
    df <- distance_table %>% dplyr::select(all_of(the_cell)) %>% rownames_to_column(var="cell")
    colnames(df) <- c("cell", "distanceT")
    expr2[[the_cell]] <- merge(df, metadata_df_basecells, by="cell") %>% slice_min(order_by=distanceT, n=1)
  }
  distances <- do.call(rbind,expr2) %>% replace(., is.na(.), 0)
  #Define distance group (measured in micrometers)
  distances_2 <- distances %>% mutate(distgroupT = case_when(distanceT==0 ~"base",
                                                             distanceT > 0 & distanceT < 100 ~ "lesion",
                                                             distanceT >= 100 & distanceT < 200 ~ "proximal",
                                                             .default = "outside")) %>%
    mutate(closestbaseT = Tlesiongroup_ann) %>%
    mutate(distgroup_closestbaseT = paste(distgroupT, closestbaseT, sep=" ")) %>%
    select(c("distanceT", "distgroupT", "closestbaseT", "distgroup_closestbaseT"))
  
  metadata_df <- cbind.data.frame(metadata_df, distances_2)
  metadata_df$cell <- rownames(metadata_df)
  expr1[[the_sample]]  <- metadata_df
}
df <- do.call(rbind.data.frame, expr1)
rownames(df) <- df$cell

#Briefly check whether nothing was messed up
identical(rownames(df), rownames(all_d15@metadata@cell))

#Transfer new annotations to the main object
all_d15@metadata@cell$Tdistance <- df$distanceT
all_d15@metadata@cell$Tdistgroup <- df$distgroupT
all_d15@metadata@cell$Tclosestbase <- df$closestbaseT
all_d15@metadata@cell$Tdistgroup_closestbase <- df$distgroup_closestbaseT

vrSpatialPlot(all_d15, group.by = "Tdistgroup", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background = "black")
ggsave("all_d15_Tdistgroup.pdf", height=8, width=7, useDingbats=FALSE)

vrSpatialPlot(all_d15, group.by = "Tdistgroup_closestbase", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay1")
ggsave("d15_WT1_Tdistgroup_closestbase.pdf")
vrSpatialPlot(all_d15, group.by = "Tdistgroup_closestbase", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay3")
ggsave("d15_WT2_Tdistgroup_closestbase.pdf")
vrSpatialPlot(all_d15, group.by = "Tdistgroup_closestbase", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay5")
ggsave("d15_WT3_Tdistgroup_closestbase.pdf")
vrSpatialPlot(all_d15, group.by = "Tdistgroup_closestbase", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay7")
ggsave("d15_KO1_Tdistgroup_closestbase.pdf")
vrSpatialPlot(all_d15, group.by = "Tdistgroup_closestbase", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay9")
ggsave("d15_KO2_Tdistgroup_closestbase.pdf")
vrSpatialPlot(all_d15, group.by = "Tdistgroup_closestbase", alpha = 1, plot.segments = TRUE, background.color = "black", assay="Assay11")
ggsave("d15_KO3_Tdistgroup_closestbase.pdf")

saveRDS(all_d15, "all_d15_mols_celltypeann_reg_lesions.rds")

