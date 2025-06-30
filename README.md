# Code for: The deubiquitinase OTUD7B ameliorates central nervous system autoimmunity by inhibiting degradation of glial fibrillary acidic protein and astrocyte hyperinflammation

This repository contains the code used for analyzing the newly generated spatial transcriptomics data, as well as reanalysis of published bulk/single-nucleus RNA-sequencing. It starts with the raw output from the Xenium device, available on GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286422). Under the individual entries (e.g., https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8727431), download the Xenium output folder as tar archive (e.g., GSM8727431_d0_KO1_slice1.tar.gz). Note that, if you untar the archive, the folder name is the original output folder name (e.g., output-XETG00046__00181777__Region_6__20240207__143847).

# Step 1: Resegmentation
Since the default nuclear expansion of 15 µm turned out to lead to fuzzy clustering of the cells, all samples were resegmented with 10 µm nuclear expansion radius using XeniumRanger 1.7 (Software from 10x Genomics), as following:
xeniumranger resegment --expansion-distance 10 --xenium-bundle $nameOfInputFolder --id dist10_"$sampleName" --localmem 90 --localcores 1
Examples for varialble names:
nameOfInputFolder="output-XETG00046__00181777__Region_6__20240207__143847"
sampleName="d0_KO1_slice1"

# Step 2: Collect and merge data
The samples were defined as individual field of view when initiating the Xenium run. In the first step, we read in all the output folders as they came from the Xenium device, and merge the individual samples to a day 0 and a day 15 object. This step is rather memory-intense.

# Step 3: Annotating cell types
In this step, we annotate cell types. This was done purely based on the transcriptome. We have also tried the BANKSY package (https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md) and the BuildNicheAssay from Seurat, however the annotation was less reliable than transcriptome-only.

