# Code for: The deubiquitinase OTUD7B ameliorates central nervous system autoimmunity by inhibiting degradation of glial fibrillary acidic protein and astrocyte hyperinflammation

This repository contains the code used for analyzing the newly generated spatial transcriptomics data, as well as reanalysis of published bulk/single-nucleus RNA-sequencing.

# Step 1: Resegmentation
Since the default nuclear expansion of 15 µm turned out to lead to fuzzy clustering of the cells, all samples were resegmented with 10 µm nuclear expansion radius using XeniumRanger 1.7 (Software from 10x Genomics), as following:
xeniumranger resegment --expansion-distance 10 --xenium-bundle $nameOfInputFolder --id dist10_"$outputName" --localmem 90 --localcores 1

# Step 2: Merging individual samples
The samples were defined as individual field of view when initiating the Xenium run. 

