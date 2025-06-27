# Code for: The deubiquitinase OTUD7B ameliorates central nervous system autoimmunity by inhibiting degradation of glial fibrillary acidic protein and astrocyte hyperinflammation

This repository contains the code used for analyzing the newly generated spatial transcriptomics data, as well as reanalysis of published bulk/single-nucleus RNA-sequencing.

# Step 1: Collect and merge data
The samples were defined as individual field of view when initiating the Xenium run. In the first step, we read in all the output folders as they came from the Xenium device, and merge the individual samples to a day 0 and a day 15 object. This step is rather memory-intense.

# Step 2: Annotating cell types
In this step, we annotate cell types. This was done purely based on the transcriptome. We have also tried the BANKSY package (https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md) and the BuildNicheAssay from Seurat, however the annotation was less reliable than transcriptome-only.

