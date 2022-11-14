# Correlation and causality assessment of a trait-endophenotype pair
From [CASTom-iGEx Module 2](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_prediction) output, the following R based pipeline performs 
1. correlation of Z-statistic association for genes and pathways between a trait and considered endophenotypes, 
2. mendelian randomization analysis to assess causality endophenotype -> trait or trait -> endophenotype using as genetic instruments imputed genes and pathways.

## Input Files
- **Imputed gene expression** (*--geneExpPred_file*): output from CASTom-iGEx (Module 2), tissue-specific imputed gene expression (predictedExpression.tar.gz)
- **Pathway-scores** (*--pathScore_file*): output from CASTom-iGEx (Module 2), tissue-specific pathway-scores (e.g. Pathway_Reactome_scores.txt)
- **Summary TWAS** (*--tscore_pheno_file*): summary statistic from TWAS in a txt file format, all tissues combined
- **Summary PALAS (Reactome)** (*--pathR_pheno_file*): summary statistic from PALAS for Reactome in a txt file format, all tissues combined
- **Summary PALAS (GO)** (*--pathGO_pheno_file*): summary statistic from PALAS for GO in a txt file format, all tissues combined
- **Path to UKBB .RData associations** (*--inputFold_rel*): complete path to TWAS and PALAS output in .RData format for UKBB pehnotypes. .RData files are obtained from *pheno_association_<>.R* scripts of CASTom-iGEx (Module 2)
- 

## Workflow
### Initial filtering if datasets are not harmonized 

- compare_geneExp_matchedDataset_run.R
- compare_pathScore_matchedDataset_run.R

### Correlation based on gene and pathway association of a trait of interest with multiple endophenotypes 

- correlation_pheno_relatedPheno_run.R

### Mendelian randomization based on genes and pathways association 

correlation_features_run.R

#### (direct)
mendelianRand_pheno_relatedPheno_run.R

#### (reverse)

mendelianRand_reverse_pheno_relatedPheno_run.R


