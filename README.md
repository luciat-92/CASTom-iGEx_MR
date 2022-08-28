# CASTom-iGEx_MR

Mendelian randomization applied to castom-igex TWAS and PALAS output

## About the project

This repository contains R scripts to perform systematic correlation and Mendelian randomization analyses based on TWAS (transcriptome-wide association study) and PALAS (pathway level association study) from [CASTom-iGEx](https://gitlab.mpcdf.mpg.de/luciat/castom-igex.git) pipeline.

The aim is to test correlation mediated via tissue-specific genetically derived genes and biological pathways between a trait of interest (e.g. CAD) and a plausible endophenotype (e.g. LDL). After an initial pruning of genes and pathways based on transcription starting site distance and jaccard index similarity, correlation of a trait-endophenotype pair is computed via Spearman correlation. Those passing nominal P-value threshold of 0.05 are then investigated for causal effect via two-sample Mendelian Randomization using random-effect inverse-weighted variance method [[1]](#1), accounting for correlation among genetic instruments. Genetic instruments are the derived imputed genes and pathways from dosages information.
Endophenotypes are obtained from UK Biobank large phenotype collection but scripts could be modified for a generla input. 

## Built with 
* R (3.5.3)
### Required R packages
- 

## References
<a id="1">[1]</a> Bowden, Jack et al. “A framework for the investigation of pleiotropy in two-sample summary data Mendelian randomization.” Statistics in medicine vol. 36,11 (2017): 1783-1802. [doi:10.1002/sim.7221](doi:10.1002/sim.7221)
