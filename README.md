# Computational approach to senescence processes in lung cancer 
This repository contains two scripts with the codes used for determining how senescence is affecting lung cancer in both, human and murine datasets. 

## LUAD-project
This script contains the analysis of a human dataset extracted from the TCGA repository (The Cancer Genome Atlas), in order to determine how is senescence affecting those patients that were under chemotherapy treatment conditions compared to the ones that were not under any treatment at all. 

The analysis were performed based on some established goals that we mention bellow: 
*1.* Determining the survival rate of those patients that were under a pharmacological treatment compared to those patients that were not under any treatment condition. 
*2.* Analysing the gene expression profiles of both conditions (treated patients vs. non-treated patients). 
*3.* Dividing the patients that were treated between the ones that have high expression of senescence biomarkers and the ones that have low expression of senescence biomarkers, we want to determine the survival rate of both conditions. 



## KRAS-project 
This script contains the analysis of a murine dataset extracted from the GEO repository (Gene Expression Omnibus), in order to determine how is senescence affecting different type of lesions: hyperplastic lesions (early injuries) and full-blown adenocarcinomas (advanced injuries). 
Hence, different approaches were considered: 
1. Dividing the dataset between both conditions mentioned above, we want to determine if there is presence of a senescence signature in early stages of cancer and therefore, in hyperplastic lesions, and if this signature is significally increased comparing with the signature of advanced stages or full-blown adenocarcinomas. 



## Observations
- Note that both scripts are supplemented with 'hashtags' explaining each code line, just to make sure every line is fully understanded and explained.
- Nevertheless, those hashtags were writen in Catalan. 
