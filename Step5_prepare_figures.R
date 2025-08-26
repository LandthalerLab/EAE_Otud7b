options(rgl.useNULL=TRUE)

#Install Voltron from github
#devtools::install_github("Artur-man/VoltRon")
#install other packages from CRAN or Bioconductor



library(Seurat)
library(harmony)
library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(VoltRon)
library(tibble)
library(tidyr)
library(igraph)
library(reshape2)
library(DESeq2)
library(apeglm)
library(dendextend)
library(stringr)
library(gridExtra)
library(ggbeeswarm)

#Get the summarySE function from http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/ , save it in summarySE.R
source("./summarySE.R")


#From https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#setwd("/Volumes/Storage/MoreApplications/stuff/Xenium/EAE/new_September_2024/")
#Redoing figures in December 2024 to include base in lesion 
setwd("/Volumes/Storage/MoreApplications/stuff/Xenium/EAE/new_December_2024/")

#Read the completely annotated object, resulting from previous processing steps or downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286422
all_d15 <- readRDS("all_d15_mols_reg_ann_lesions_v4.rds")
all_d0 <- readRDS("all_d0_mols_reg_ann_2.rds")

assay_names <- as.list(as.vector(all_d15@sample.metadata %>% filter(Assay=="Xenium") %>% rownames()))
names(assay_names) <- all_d15@sample.metadata %>% filter(Assay=="Xenium") %>% pull(var="Sample")

the_cols = c(neurons = "#CD9600", astrocytes = "#7CAE00", oligodendrocytes = "#00BF75", microglia = "#fc8d86", `T cells` = "#FF689F", `endothelial cells` = "#00BCD6", 
  pericytes = "#35A2FF", `ventral leptomeningeal cells` = "gray78")
the_cols_d0 = c(neurons = "#CD9600", astrocytes = "#7CAE00", oligodendrocytes = "#00BF75", microglia = "#fc8d86", `endothelial cells` = "#00BCD6", 
             pericytes = "#35A2FF", `ventral leptomeningeal cells` = "gray78")



#Fig. 1E: overview of an example WT section
#This is done in a rather complicated way, better do the cropping etc. with VoltRon
vrSpatialPlot(all_d15, assay="Assay5", colors=the_cols, group.by = "celltype", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background = c("main", "H&E"), scale.image = TRUE)
ggsave("d15_WT3_celltypes_legend.pdf")
vrSpatialPlot(all_d15, assay="Assay5", colors=the_cols, group.by = "celltype", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background = c("main", "H&E"), scale.image = FALSE)+NoLegend()
ggsave("d15_WT3_celltypes_nolegend.pdf")
vrSpatialPlot(all_d15, assay="Assay5", group.by = "Tdistgroup", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background = c("main", "H&E"), scale.image = TRUE)+NoLegend()
ggsave("d15_WT3_Tdistgroup_nolegend.pdf")
#Processing in Illustrator
#Rotate 20 degrees
#Remove all clipping mask etc.
#Place shape 240x300 pt
#Select everything without the background H&E image
#Pathfinder > Crop
#Copy remaining cell shapes to new file
#Undo, remove everything except shape and image
#Crop image, using smart guides crop frame should snap to shape, copy cropped image as well
#In the Tdistgroup plot, outline lesion regions, copy them to the file with the cropped material, align using the cells (size and proportion should be the same) 


assay_names_d0 <- as.list(as.vector(all_d0@sample.metadata %>% filter(Assay=="Xenium") %>% rownames()))
names(assay_names_d0) <- all_d0@sample.metadata %>% filter(Assay=="Xenium") %>% pull(var="Sample")

vrSpatialPlot(all_d0, assay="Assay5", colors=the_cols_d0, group.by = "celltype", pt.size = 0.5, alpha = 1, plot.segments = TRUE, background = c("main", "H&E"), scale.image = FALSE)+NoLegend()
ggsave("d0_WT3_celltypes_nolegend.pdf")
#Rotate 270 degrees, process as above (without outlines)



#Fig. 1G left: number of cells in WT sections
df1 <- all_d0@metadata@cell %>% select(Sample) %>% mutate(timepoint="d0") %>% mutate(condition=str_match(Sample, "KO|WT"))
df2 <- all_d15@metadata@cell %>% select(Sample) %>% mutate(timepoint="d15") %>% mutate(condition=str_match(Sample, "KO|WT"))
c = rbind.data.frame(df1, df2) %>% group_by(timepoint, Sample) %>% tally() %>% as.data.frame() %>%
  mutate(condition=str_match(Sample, "KO|WT")) %>% filter(condition=="WT")
tgc <- summarySE(c, measurevar="n", groupvars=c("timepoint", "condition"))
ggplot(tgc, aes(x=timepoint, y=n, fill=condition))+
  geom_bar(position=position_dodge(.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=timepoint, y=n, fill=condition), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=n-sd, ymax=n+sd), width=.2, position=position_dodge(.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("Cells per timepoint and genotype", sep=" "))
ggsave("number_of_cells_onlyWT.pdf")



#Fig. 1G right, stacked barplot with relative cell type composition (WT and KO separately, but generated in parallel as below)
#This code also can be used for side-by-side barplots, which show the same data but less compact
df1 <- all_d0@metadata@cell %>% select(Sample, celltype) %>% mutate(Tdistgroup="all") %>% mutate(timepoint="d0") %>%
  mutate(condition=str_match(Sample, "KO|WT")) %>% filter(condition=="WT")
df2 <- all_d15@metadata@cell %>% select(Sample, celltype, Tdistgroup) %>% mutate(timepoint="d15") %>%
  mutate(condition=str_match(Sample, "KO|WT")) %>% filter(condition=="WT") %>% mutate(Tdistgroup=ifelse(Tdistgroup=="base", "lesion", Tdistgroup))
a = rbind.data.frame(df1, df2) %>% group_by(timepoint, Tdistgroup, Sample) %>% tally(name="tot")  
b = rbind.data.frame(df1, df2) %>% group_by(timepoint, Tdistgroup, Sample, celltype) %>% tally(name="pos") %>% ungroup() %>% as.data.frame()
#Complete rows 
b = b %>% tidyr::complete(tidyr::expand(b, nesting(timepoint, Tdistgroup, Sample)), celltype=unique(b$celltype), fill=list(pos=0))
c =
  left_join(a , b , by = c('timepoint', 'Tdistgroup', 'Sample')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot) %>%
  mutate(celltype = forcats::fct_relevel(celltype, names(the_cols)))
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("timepoint", "Tdistgroup", "celltype"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, names(the_cols)))
#Barplot
p3 <- ggplot(tgc, aes(x=factor(paste(timepoint, Tdistgroup, sep=" "), levels=c("d0 all", "d15 base", "d15 lesion", "d15 proximal", "d15 outside")), y=fraction, fill=celltype))+
  geom_bar(position=position_dodge(0.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=paste(timepoint, Tdistgroup, sep=" "), y=fraction, fill=celltype), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(0.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_cols)+
  ylim(0,0.8)+
  ggtitle("WT")
#Stacked bar plot
tgc_2 <- tgc %>% dplyr::ungroup() %>% dplyr::group_by(Tdistgroup) %>% dplyr::mutate(SDPos = rev(cumsum(rev(fraction))))
p1 <- ggplot(tgc_2, aes(x=factor(paste(timepoint, Tdistgroup, sep=" "), levels=c("d0 all", "d15 base", "d15 lesion", "d15 proximal", "d15 outside")), y=fraction, fill=celltype))+
  geom_bar(stat="identity", size=0.2)+
  geom_errorbar(aes(ymin=SDPos-sd, ymax=SDPos), width=0, position="identity", size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_cols)+
  ggtitle("WT")
#ggsave("percentage_of_cells_onlyWT_bydistgroup.pdf")

#Now with KO
df1 <- all_d0@metadata@cell %>% select(Sample, celltype) %>% mutate(Tdistgroup="all") %>% mutate(timepoint="d0") %>%
  mutate(condition=str_match(Sample, "KO|WT")) %>% filter(condition=="KO")
df2 <- all_d15@metadata@cell %>% select(Sample, celltype, Tdistgroup) %>% mutate(timepoint="d15") %>%
  mutate(condition=str_match(Sample, "KO|WT")) %>% filter(condition=="KO") %>% mutate(Tdistgroup=ifelse(Tdistgroup=="base", "lesion", Tdistgroup))
a = rbind.data.frame(df1, df2) %>% group_by(timepoint, Tdistgroup, Sample) %>% tally(name="tot")  
b = rbind.data.frame(df1, df2) %>% group_by(timepoint, Tdistgroup, Sample, celltype) %>% tally(name="pos") %>% ungroup() %>% as.data.frame()
#Complete rows 
b = b %>% tidyr::complete(tidyr::expand(b, nesting(timepoint, Tdistgroup, Sample)), celltype=unique(b$celltype), fill=list(pos=0))
c =
  left_join(a , b , by = c('timepoint', 'Tdistgroup', 'Sample')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot) %>%
  mutate(celltype = forcats::fct_relevel(celltype, names(the_cols)))
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("timepoint", "Tdistgroup", "celltype"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, names(the_cols)))
#Barplot
p4 <- ggplot(tgc, aes(x=factor(paste(timepoint, Tdistgroup, sep=" "), levels=c("d0 all", "d15 base", "d15 lesion", "d15 proximal", "d15 outside")), y=fraction, fill=celltype))+
  geom_bar(position=position_dodge(0.9), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=paste(timepoint, Tdistgroup, sep=" "), y=fraction, fill=celltype), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(0.9), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_cols)+
  ylim(0,0.8)+
  ggtitle("KO")
#Stacked bar plot
tgc_2 <- tgc %>% dplyr::ungroup() %>% dplyr::group_by(Tdistgroup) %>% dplyr::mutate(SDPos = rev(cumsum(rev(fraction))))
p2 <- ggplot(tgc_2, aes(x=factor(paste(timepoint, Tdistgroup, sep=" "), levels=c("d0 all", "d15 base", "d15 lesion", "d15 proximal", "d15 outside")), y=fraction, fill=celltype))+
  geom_bar(stat="identity", size=0.2)+
  geom_errorbar(aes(ymin=SDPos-sd, ymax=SDPos), width=0, position="identity", size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_cols)+
  ggtitle("KO")

pdf("percentage_of_cells_bydistgroup.pdf", width=10, height =5)
grid.arrange(p1, p2, nrow=1)
dev.off()

pdf("percentage_of_cells_bydistgroup_dodged.pdf", width=15, height =5)
grid.arrange(p3, p4, nrow=1)
dev.off()

#Fig. 1F and 1H: differential expression of Otud7b and other genes in astrocytes in specific areas (Tdistgroups) at day 15 compared to the entire sample at day 0
#count table and column data for d0
celltypes_to_check=c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial cells", "pericytes")
samples_to_check <- unique(all_d0@metadata@cell$Sample)
expr <- list()
count_data <- vrData(all_d0)
for (the_celltype in celltypes_to_check) {
  for (the_sample in samples_to_check) {
    cells <- all_d0@metadata@cell %>% filter(celltype==the_celltype) %>% filter(Sample==the_sample) %>%
      rownames()
    if (length(cells) > 10) { 
      expr[[paste0(the_celltype,'|',the_sample,'|all')]] <- rowSums(count_data[,cells])
    }
  }
}
counts_d0 <- do.call(cbind,expr)
colData_d0 <- data.frame(condition=factor(str_match(names(expr), "KO|WT")),
                         sample=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\2',names(expr), perl=TRUE)),
                         celltype=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\1',names(expr), perl=TRUE)),
                         distgroup="all",
                         group=names(expr),
                         row.names=names(expr))

#d15, all celltypes, not separated by area
celltypes_to_check=c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial cells", "pericytes")
samples_to_check <- unique(all_d15@metadata@cell$Sample)
expr <- list()
count_data <- vrData(all_d15)
for (the_celltype in celltypes_to_check) {
  for (the_sample in samples_to_check) {
    cells <- all_d15@metadata@cell %>% filter(celltype==the_celltype) %>% filter(Sample==the_sample) %>%
      rownames()
    if (length(cells) > 10) { 
      expr[[paste0(the_celltype,'|',the_sample,'|all')]] <- rowSums(count_data[,cells])
    }
  }
}
counts_d15 <- do.call(cbind,expr)
colData_d15 <- data.frame(condition=factor(str_match(names(expr), "KO|WT")),
                          sample=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\2',names(expr), perl=TRUE)),
                          celltype=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\1',names(expr), perl=TRUE)),
                          distgroup="all",
                          group=names(expr),
                          row.names=names(expr))


#count table and column data for d15 separated between the lesion areas (Tdistgroups).
#Note that in the manuscript and figures, naming is a bit different.
#The "lesion core" there is "base" and "lesion" here together, "proximal" is "lesion rim" and "outside" is "peri-lesion"
celltypes_to_check=c("neurons", "oligodendrocytes", "astrocytes", "microglia", "endothelial cells", "pericytes")
samples_to_check <- unique(all_d15@metadata@cell$Sample)
expr <- list()
count_data <- vrData(all_d15)
for (the_celltype in celltypes_to_check) {
  for (the_sample in samples_to_check) {
    for (the_distgroup in c("base", "lesion", "proximal", "outside")) {
      cells <- all_d15@metadata@cell %>% filter(celltype==the_celltype) %>% filter(Sample==the_sample) %>%
        filter(grepl(the_distgroup, Tdistgroup)) %>% rownames()
      if (length(cells) > 10) { 
        expr[[paste0(the_celltype,'|',the_sample,'|',the_distgroup, '|', 'all')]] <- rowSums(count_data[,cells])
      }
    }
    cells <- all_d15@metadata@cell %>% filter(celltype==the_celltype) %>% filter(Sample==the_sample) %>%
      filter(grepl("base|lesion", Tdistgroup)) %>% rownames()
    if (length(cells) > 10) { 
      expr[[paste0(the_celltype,'|',the_sample,'|','baseandlesion', '|', 'all')]] <- rowSums(count_data[,cells])
    }
  }
}
counts_Tdistgroups_d15 <- do.call(cbind,expr)
colData_Tdistgroups_d15 <- data.frame(condition=factor(str_match(names(expr), "KO|WT")),
                                      celltype=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\1',names(expr), perl=TRUE)),
                                      sample=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\2',names(expr), perl=TRUE)),
                                      distgroup=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\3',names(expr), perl=TRUE)),
                                      group=names(expr),
                                      row.names=names(expr))

counts_merged_all <- cbind.data.frame(counts_d0, counts_Tdistgroups_d15)
colData_all <- rbind.data.frame(colData_d0, colData_Tdistgroups_d15)
colData_all$timepoint = factor(str_match(colData_all$sample, "d0|d15"))

res <- list()
for (the_condition in c("WT", "KO")) {
  for (the_celltype in celltypes_to_check) {
    for (the_distgroup in c("base", "baseandlesion", "lesion", "outside", "proximal")) {
      take_row <- rowSums(counts_merged_all) > 0
      take_col <- (colSums(counts_merged_all) > 0) & (colData_all[,'celltype']==the_celltype &
                                                        colData_all[,'condition']==the_condition &
                                                        ((colData_all[,'timepoint']=="d0" & colData_all[,'distgroup'] == "all") |
                                                           (colData_all[,'timepoint']=="d15" & colData_all[,'distgroup'] == the_distgroup)))
      try({
        dds <- DESeqDataSetFromMatrix(countData=counts_merged_all[take_row,take_col],
                                      colData=colData_all[take_col,,drop=FALSE],
                                      design=~timepoint)
        dds <- DESeq(dds)
        for (the_comparison in grep("Intercept", resultsNames(dds), value=TRUE, invert = TRUE)) {
          #Do comparison only when there are three entries in both regions
          nos=table(take_col)[[2]]
          if (nos==6) {
            res[[paste(the_condition,the_celltype,the_distgroup,the_comparison, sep='_')]] <- lfcShrink(dds,
                                                                                                        coef=the_comparison,
                                                                                                        type='apeglm',
                                                                                                        format="DataFrame") %>%
              as.data.frame() %>%
              tibble::rownames_to_column('gene_name') %>%
              dplyr::mutate(celltype=the_celltype,contrast=the_comparison,condition=the_condition,distgroup=the_distgroup)
          }
        }
      })
    }
  }
}

pseudobulk_d0_d15_distgroupcomparison <- do.call(rbind,res)
saveRDS(pseudobulk_d0_d15_distgroupcomparison, "pseudobulk_d0_d15_distgroupcomparison_new.rds")
pseudobulk_d0_d15_distgroupcomparison <- readRDS("pseudobulk_d0_d15_distgroupcomparison_new.rds")

gene_list <- list()
gene_list[["inflmgenes"]][["genes_of_interest"]] <- c("Isg15", "Cxcl9", "Cxcl10", "Ccl2", "Il1b", "Il6")
gene_list[["inflmgenes"]][["w"]] <- 3
gene_list[["inflmgenes"]][["h"]] <- 3
gene_list[["Otud7b"]][["genes_of_interest"]] <- c("Otud7b")
gene_list[["Otud7b"]][["w"]] <- 3
gene_list[["Otud7b"]][["h"]] <- 3

contrasts_to_use <- unique(pseudobulk_d0_d15_distgroupcomparison$contrast)
distgroups_to_use <- c("baseandlesion", "proximal", "outside")
the_celltype="astrocytes"
the_condition = "WT"

for (the_name in names(gene_list)) {
  genes_of_interest <- gene_list[[the_name]][["genes_of_interest"]]
  w <- gene_list[[the_name]][["w"]]
  h <- gene_list[[the_name]][["h"]]
  
  the_filename <- paste0("d15vsd0_pseudobulk_distgroups_", the_name, "_", the_condition, "_", the_celltype ,".pdf")
  #prefilter the pseudobulk data frame
  df <- pseudobulk_d0_d15_distgroupcomparison %>%
    dplyr::filter(celltype == the_celltype) %>%
    dplyr::filter(condition == the_condition) %>%
    dplyr::filter(distgroup %in% distgroups_to_use) %>%
    dplyr::filter(gene_name %in% genes_of_interest) %>%
    dplyr::mutate(padj=ifelse(is.na(padj), 1, padj))
  
  #Clustering to get order of genes and groups
  df2 <- df %>%
    dplyr::mutate(group=paste(celltype,condition,distgroup,contrast,sep='_'))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  df2[is.na(df2)] <- 0
  if (length(genes_of_interest) > 1) {
    hc <- hclust(dist(df2))
    gene.order <- row.names(df2)[order.hclust(hc)]
  } else {
    gene.order <- genes_of_interest
  }
  #plot
  minmax <- as.integer(max(max(df2), abs(min(df2))))+1
  ggplot(df %>%
           dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
         aes(y=gene,x=factor(distgroup, levels=distgroups_to_use),size=-log10(padj),color=log2FoldChange)) +
    geom_point(shape=16) +
    #scale_size_continuous(name='adj. p-value',
    #                      breaks=c(2,4,6),
    #                      labels=c("1e-2","1e-4", "<1e-6"),
    #                      limits=c(0,6)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-minmax,minmax),oob=scales::squish) +
    facet_wrap(~celltype, nrow=1) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')+
    ggtitle(the_celltype)
  ggsave(the_filename, width=w, height = h, units="in", useDingbats=FALSE)
}

#Figure 2A: differential expression between WT and KO at d0 per celltype
res <- list()
for (the_celltype in celltypes_to_check) {
  take_row <- rowSums(counts_d0) > 0
  take_col <- (colSums(counts_d0) > 0) & (colData_d0[,'celltype']==the_celltype)
  try({
    dds <- DESeqDataSetFromMatrix(countData=counts_d0[take_row,take_col],
                                  colData=colData_d0[take_col,,drop=FALSE],
                                  design=~condition)
    dds <- DESeq(dds)
    for (the_comparison in grep("Intercept", resultsNames(dds), value=TRUE, invert = TRUE)) {
      #Do comparison only when there are three entries in both regions
      nos=table(take_col)[[2]]
      if (nos==6) {
        res[[paste(the_condition,the_celltype,the_distgroup,the_comparison, sep='_')]] <- lfcShrink(dds,
                                                                                                    coef=the_comparison,
                                                                                                    type='apeglm',
                                                                                                    format="DataFrame") %>%
          as.data.frame() %>%
          tibble::rownames_to_column('gene_name') %>%
          dplyr::mutate(celltype=the_celltype,contrast=the_comparison)
      }
    }
  })
}

pseudobulk_d0_WTvsKO_comparison <- do.call(rbind,res)
saveRDS(pseudobulk_d0_WTvsKO_comparison, "pseudobulk_d0_WTvsKO_comparison.rds")

res <- list()
for (the_celltype in celltypes_to_check) {
  take_row <- rowSums(counts_d15) > 0
  take_col <- (colSums(counts_d15) > 0) & (colData_d15[,'celltype']==the_celltype)
  try({
    dds <- DESeqDataSetFromMatrix(countData=counts_d15[take_row,take_col],
                                  colData=colData_d15[take_col,,drop=FALSE],
                                  design=~condition)
    dds <- DESeq(dds)
    for (the_comparison in grep("Intercept", resultsNames(dds), value=TRUE, invert = TRUE)) {
      #Do comparison only when there are three entries in both regions
      nos=table(take_col)[[2]]
      if (nos==6) {
        res[[paste(the_condition,the_celltype,the_distgroup,the_comparison, sep='_')]] <- lfcShrink(dds,
                                                                                                    coef=the_comparison,
                                                                                                    type='apeglm',
                                                                                                    format="DataFrame") %>%
          as.data.frame() %>%
          tibble::rownames_to_column('gene_name') %>%
          dplyr::mutate(celltype=the_celltype,contrast=the_comparison)
      }
    }
  })
}

pseudobulk_d15_WTvsKO_comparison <- do.call(rbind,res)
saveRDS(pseudobulk_d15_WTvsKO_comparison, "pseudobulk_d0_WTvsKO_comparison.rds")

#Comparison Otud7b WT vs. KO
pseudobulk_d15_WTvsKO_comparison$timepoint <- "d15"
pseudobulk_d0_WTvsKO_comparison$timepoint <- "d0"
df <- rbind.data.frame(pseudobulk_d0_WTvsKO_comparison, pseudobulk_d15_WTvsKO_comparison) %>% filter (gene_name=="Otud7b")
minmax <- as.integer(max(max(df$log2FoldChange), abs(min(df$log2FoldChange))))+1
ggplot(df,
       aes(y=celltype,x=factor(timepoint, levels=c("d0", "d15")),size=-log10(padj),color=log2FoldChange)) +
  geom_point(shape=16) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-minmax,minmax),oob=scales::squish) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')+
  ggtitle("Otud7b")
ggsave("Otud7b.pdf", width=3.5, height = 2, units="in", useDingbats=FALSE)

#Figure 2I: comparison of areas with d0, for both WT and KO
#We use the dataframe from above
contrasts_to_use <- unique(pseudobulk_d0_d15_distgroupcomparison$contrast)
distgroups_to_use <- c("baseandlesion", "proximal", "outside")
the_celltype="astrocytes"
w=3
h=2.5
for (the_name in c("Gfap")) {
  genes_of_interest <- the_name
  
  the_filename <- paste0("d15vsd0_pseudobulk_distgroups_", the_name, "_", the_celltype ,".pdf")
  #prefilter the pseudobulk data frame
  df <- pseudobulk_d0_d15_distgroupcomparison %>%
    dplyr::filter(celltype == the_celltype) %>%
    dplyr::filter(distgroup %in% distgroups_to_use) %>%
    dplyr::filter(gene_name %in% genes_of_interest) %>%
    dplyr::mutate(padj=ifelse(is.na(padj), 1, padj))
  
  #plot
  minmax <- as.integer(max(max(df$log2FoldChange), abs(min(df$log2FoldChange))))+1
  ggplot(df %>%
           dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
         aes(y=condition,x=factor(distgroup, levels=distgroups_to_use),size=-log10(padj),color=log2FoldChange)) +
    geom_point(shape=16) +
    #scale_size_continuous(name='adj. p-value',
    #                      breaks=c(2,4,6),
    #                      labels=c("1e-2","1e-4", "<1e-6"),
    #                      limits=c(0,6)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-minmax,minmax),oob=scales::squish) +
    facet_wrap(~celltype, nrow=1) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')+
    ggtitle(the_celltype)
  ggsave(the_filename, width=w, height = h, units="in", useDingbats=FALSE)
}

#Figure 7A: expression plot for Tnf (circle around dots was removed in Adobe Illustrator)
all_d0_seu <- VoltRon::as.Seurat(all_d0, cell.assay = "Xenium", type = "image")
all_d15_seu <- VoltRon::as.Seurat(all_d15, cell.assay = "Xenium", type = "image")
all_seu <- merge(all_d0_seu, all_d15_seu)
all_seu <- NormalizeData(all_seu)
all_seu@meta.data$condition=str_match(all_seu@meta.data$Sample, "KO|WT")
all_seu@meta.data$timepoint=str_match(all_seu@meta.data$Sample, "d0|d15")
all_seu@meta.data$Tdistgroupred <- ifelse(all_seu@meta.data$Tdistgroup == "base" | all_seu@meta.data$Tdistgroup == "lesion", "baseandlesion", all_seu@meta.data$Tdistgroup)
all_seu@meta.data$timepoint_condition_celltype_Tdistgroupred = paste(all_seu@meta.data$timepoint, all_seu@meta.data$condition, all_seu@meta.data$celltype, all_seu@meta.data$Tdistgroupred, sep="|")
Idents(all_seu) <- all_seu@meta.data$timepoint_condition_celltype_Tdistgroupred
genes_files <- list()
genes_files[["expr_values_Tnf_2.pdf"]] <- c("Tnf")
celltype_levels <- c("astrocytes", "endothelial cells", "microglia", "neurons", "oligodendrocytes", "pericytes", "T cells", "ventral leptomeningeal cells")
for (the_filename in names(genes_files)) {
  p <- DotPlot(all_seu, group.by='timepoint_condition_celltype_Tdistgroupred', cols='RdYlBu',
               features = genes_files[[the_filename]]) +
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  #ggsave(the_filename, width=8, height=15)
  df <- p$data
  df <- df %>%
    separate_wider_delim(., cols=id, "|", names=c("timepoint", "condition", "celltype", "Tdistgroupred")) %>%
    dplyr::mutate(Tdistgroupred = factor(Tdistgroupred, levels=rev(c("NA", "baseandlesion", "proximal", "outside"))))
  
  ggplot()+geom_point(data=df, aes(x=factor(celltype, levels=celltype_levels), y=Tdistgroupred, size=pct.exp, fill=avg.exp.scaled), pch=21)+
    scale_fill_distiller(palette="RdYlBu", direction=-1)+
    theme_bw()+
    facet_wrap(~condition)+
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10)+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
  ggsave(the_filename, width=7, height=3)
}



#Figure 7B and 8A: compare lesion areas against d0 for all cells together
#d0 first
samples_to_check <- unique(all_d0@metadata@cell$Sample)
expr <- list()
count_data <- vrData(all_d0)
for (the_sample in samples_to_check) {
  cells <- all_d0@metadata@cell %>% filter(Sample==the_sample) %>%
    rownames()
  if (length(cells) > 10) { 
    expr[[paste0(the_sample,'|all')]] <- rowSums(count_data[,cells])
  }
}
counts_d0 <- do.call(cbind,expr)
colData_d0 <- data.frame(condition=factor(str_match(names(expr), "KO|WT")),
                         sample=factor(gsub('^([^\\|]*)\\|([^\\|]*)$','\\1',names(expr), perl=TRUE)),
                         distgroup="all",
                         group=names(expr),
                         row.names=names(expr))

#d15
samples_to_check <- unique(all_d15@metadata@cell$Sample)
expr <- list()
count_data <- vrData(all_d15)
for (the_sample in samples_to_check) {
  for (the_distgroup in c("proximal", "outside")) {
    cells <- all_d15@metadata@cell %>% filter(Sample==the_sample) %>% filter(grepl(the_distgroup, Tdistgroup)) %>% rownames()
    if (length(cells) > 10) { 
      expr[[paste0(the_sample,'|',the_distgroup, '|', 'all')]] <- rowSums(count_data[,cells])
    }
  }
  cells <- all_d15@metadata@cell %>% filter(Sample==the_sample) %>% filter(grepl("base|lesion", Tdistgroup)) %>% rownames()
  if (length(cells) > 10) { 
    expr[[paste0(the_sample,'|','baseandlesion', '|', 'all')]] <- rowSums(count_data[,cells])
  }
}
counts_Tdistgroups_d15 <- do.call(cbind,expr)
colData_Tdistgroups_d15 <- data.frame(condition=factor(str_match(names(expr), "KO|WT")),
                                      sample=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\1',names(expr), perl=TRUE)),
                                      distgroup=factor(gsub('^([^\\|]*)\\|([^\\|]*)\\|([^\\|]*)$','\\2',names(expr), perl=TRUE)),
                                      group=names(expr),
                                      row.names=names(expr))

counts_merged_all <- cbind.data.frame(counts_d0, counts_Tdistgroups_d15)
colData_all <- rbind.data.frame(colData_d0, colData_Tdistgroups_d15)
colData_all$timepoint = factor(str_match(colData_all$sample, "d0|d15"))

res <- list()
for (the_condition in c("WT", "KO")) {
  for (the_distgroup in c("baseandlesion", "outside", "proximal")) {
    take_row <- rowSums(counts_merged_all) > 0
    take_col <- (colSums(counts_merged_all) > 0) & (colData_all[,'condition']==the_condition &
                                                      ((colData_all[,'timepoint']=="d0" & colData_all[,'distgroup'] == "all") |
                                                         (colData_all[,'timepoint']=="d15" & colData_all[,'distgroup'] == the_distgroup)))
    try({
      dds <- DESeqDataSetFromMatrix(countData=counts_merged_all[take_row,take_col],
                                    colData=colData_all[take_col,,drop=FALSE],
                                    design=~timepoint)
      dds <- DESeq(dds)
      for (the_comparison in grep("Intercept", resultsNames(dds), value=TRUE, invert = TRUE)) {
        #Do comparison only when there are three entries in both regions
        nos=table(take_col)[[2]]
        if (nos==6) {
          res[[paste(the_condition,the_distgroup,the_comparison, sep='_')]] <- lfcShrink(dds,
                                                                                         coef=the_comparison,
                                                                                         type='apeglm',
                                                                                         format="DataFrame") %>%
            as.data.frame() %>%
            tibble::rownames_to_column('gene_name') %>%
            dplyr::mutate(contrast=the_comparison,condition=the_condition,distgroup=the_distgroup)
        }
      }
    })
  }
}

pseudobulk_d0_d15_distgroupcomparison_allcelltypesmerged <- do.call(rbind,res)
saveRDS(pseudobulk_d0_d15_distgroupcomparison_allcelltypesmerged, "pseudobulk_d0_d15_distgroupcomparison_allcelltypesmerged.rds")
contrasts_to_use <- unique(pseudobulk_d0_d15_distgroupcomparison_allcelltypesmerged$contrast)
for (genes in c("Tnf", "Il6")) {
  the_name = paste("allcelltypes", genes, sep="_")
  
  df <- pseudobulk_d0_d15_distgroupcomparison_allcelltypesmerged %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::mutate(group=paste(condition,distgroup,contrast,sep='_'))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  df[is.na(df)] <- 0
  if (length(genes) > 1) {
    hc <- hclust(dist(df))
    gene.order <- row.names(df)[order.hclust(hc)]
  } else {
    gene.order <- genes
  }
  #plot
  minmax <- as.integer(max(max(df), abs(min(df))))+1
  ggplot(pseudobulk_d0_d15_distgroupcomparison_allcelltypesmerged %>%
           dplyr::filter(gene_name %in% genes) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::mutate(gene=factor(gene_name,levels=gene.order)),
         aes(y=condition,x=factor(distgroup, levels=c("baseandlesion", "proximal", "outside")),size=-log10(padj),color=log2FoldChange)) +
    geom_point(shape=16) +
    #scale_size_continuous(name='adj. p-value',
    #                      breaks=c(2,4,6),
    #                      labels=c("1e-2","1e-4", "<1e-6"),
    #                      limits=c(0,6)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-minmax,minmax),oob=scales::squish) +
    facet_wrap(~gene, nrow=1) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  ggsave(paste0("d15vsd0_pseudobulk_distgroups_", the_name ,".pdf"), width=3, height = 2, units="in", useDingbats=FALSE)
}


