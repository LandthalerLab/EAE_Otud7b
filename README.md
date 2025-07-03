# Code for: The deubiquitinase OTUD7B ameliorates central nervous system autoimmunity by inhibiting degradation of glial fibrillary acidic protein and astrocyte hyperinflammation

This repository contains the code used for analyzing the newly generated spatial transcriptomics data, as well as reanalysis of published bulk/single-nucleus RNA-sequencing. It starts with the raw output from the Xenium device, available on GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286422). Under the individual entries (e.g., https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8727431), download the Xenium output folder as tar archive (e.g., GSM8727431_d0_KO1_slice1.tar.gz). Note that, if you untar the archive, the folder name is the original output folder name (e.g., output-XETG00046__00181777__Region_6__20240207__143847).
From the GEO repository, the fully annotated VoltRon objects used in the paper can be downloaded (tar files GSE286422_all_d0_mols_reg_ann_2.rds.gz for day 0 and GSE286422_all_d15_mols_reg_ann_lesions_v4.rds.gz at day 15 after administration of the MOG peptide). 

# Step 1: resegmentation
Since the default nuclear expansion of 15 µm turned out to lead to fuzzy clustering of the cells, all samples were resegmented with 10 µm nuclear expansion radius using XeniumRanger 1.7 (Software from 10x Genomics), as following:
xeniumranger resegment --expansion-distance 10 --xenium-bundle $nameOfInputFolder --id dist10_"$sampleName" --localmem 90 --localcores 1
Examples for varialble names:
nameOfInputFolder="output-XETG00046__00181777__Region_6__20240207__143847"
sampleName="d0_KO1_slice1"

# Step 1: collect and merge data
The samples were defined as individual field of view when initiating the Xenium run. In the first step, we read in all the output folders as they came from the Xenium device, and merge the individual samples to a day 0 and a day 15 object. This step is rather memory-intense.

# Step 2: annotating cell types
In this step, we annotate cell types. This was done manually and purely based on the transcriptome, i.e. based on transcriptome-derived clusters. We have also tried the BANKSY package (https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md) and the BuildNicheAssay from Seurat, however the annotation was less reliable than transcriptome-only. Note that depending on the clustering and manual annotation, results may vary.

# Step 3: registration of H&E images
In this step, manual image registration for the H&E images was performed. The H&E images come from the staining of the Xenium slides following the Xenium run, i.e. were performed on the same section as the spatial transcriptomics.

# Step 4: annotating regions and lesion
For the sections from day 0, we annotated only dorsal and ventral regions, which can be later used for differential expression analysis. For day 15, we annotated lesion. This was done as following – first, groups of T cells are defined, however not based on the celltype definition from above but based on whether they expression T cell marker genes (Cd3e, Cd4, Cd8a). In order to form a group, at least n cells have to be within k µm (variables intra_lesion_distance and min_group_size). T cells outside of such groups are in "small clusters" and not taken into account. Second, the T cell groups are spatially expanded, i.e. proximal and distal regions are defined. At the end of this script, plots are generated with cells colored by lesion/region. However, these are not the lesions/regions that were used in the paper, as we subsequently visually examined these lesions in H&E staining. Some of the computationally defined lesions were excluded from analysis, and others merged into larger lesions. Our annotation can be transfered from the objects available on GEO. 
