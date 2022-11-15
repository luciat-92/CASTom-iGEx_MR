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

## Workflow
### Initial filtering if datasets are not harmonized 

Compute genes correlation imputed from 2 different models. Genes are imputed on the reference panel from which the gene expression models are estimated, see [CASTom-iGEx Module 1](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_training)

```sh
Rscript ./compare_geneExp_matchedDataset_run.R \
    --geneExpPred_file (2 files necessary) \
    --corr_thr (default 0.8) \
    --tissue_name \
    --outFold 
```
The output includes (saved in *--outFold*):
- <tissue_name>_filter_genes_matched_datasets.txt 


Compute pathways correlation imputed from 2 different models. Pathways-scores are computed on the reference panel from which the gene expression models are estimated (see [CASTom-iGEx Module 1](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_training)), after the computation of gene T-scores (see [CASTom-iGEx Module 2](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_prediction))

```sh
Rscript ./compare_pathScore_matchedDataset_run.R \
    --pathScore_file (2 files necessary) \
    --corr_thr (default 0.8) \
    --type_path \
    --tissue_name \
    --outFold 
```
The output includes (saved in *--outFold*):
- <tissue_name>_filter_path_<type_path>_matched_datasets.txt

### Correlation based on gene and pathway association of a trait of interest with multiple endophenotypes 

- *--phenoFold*: includes phenotype info, input files also available in Usage/Data/
- *--refFold*: path to folder including gene TSS annotation (use [CASTom-iGEx refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData))
- *--feat_filt*: if not NULL, path to genes, Reactome and GO files for feature filtering based on correlation among models 

```sh
./correlation_pheno_relatedPheno_run.R \
	--phenoFold  \
	--tissue_name \
	--pheno_name_comp (e.g. SCZ) \
	--inputFold_rel \
	--tscore_pheno_file \
	--pathR_pheno_file \
	--pathGO_pheno_file \
	--outFold \
	--pval_FDR_pheno (default 0.05) \
	--pval_FDR_rel (default 0.05) \
	--perc_par (default 0.1) \
	--thr_dist_par (default 250000) \
	--refFold \
	--feat_filt (default NULL)
```
The output includes (saved in *--outFold*):
- correlation_enrich_<pheno_name_comp>_relatedPheno.RData
list including correlation and pvalue between the trait (pheno_name_comp) and the UKBB endophenotypes for genes and pathways separately


### Mendelian randomization based on genes and pathways association 

- correlation_features_run.R

#### (direct)
- mendelianRand_pheno_relatedPheno_run.R

#### (reverse)

- mendelianRand_reverse_pheno_relatedPheno_run.R


