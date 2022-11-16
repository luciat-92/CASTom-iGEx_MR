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

#### Combine results across tissues and filter significant correlations

- *--corrRes_tissue_file* vector including .RData object from the previous command, one per tissue
- *--tissue_name* tissues name, must be same length as *--corrRes_tissue_file*

```sh
./combine_correlation_allTissues_run.R \
	--corrRes_tissue_file \
	--tissue_name \
	--type_data \
	--thr_pval (default 0.05) \
	--outFold 
```
The output includes (saved in *--outFold*):
- <type_data>_correlation_allTissues.txt, <type_data>_correlation_pvalue_allTissues.txt, <type_data>_correlation_sign_pval%s_allTissues.txt
phenotypes x tissues matrices with Spearman correlation, p-value, and p-value only for nominal significant results

### Mendelian randomization (MR) based on genes and pathways association 

#### Compute features correlation

Compute features (genes or pathways) correlation to be used in MR
- *--inputFile* pathway-score or gene T-scores (.RData) from which the correlation is computed. If pathway-scores, it must be a vector of 2 .RData, one for Reactome and one for GO
- *--sampleAnnFile* samples matching *--inputFile* or a subset used to compute features correlation
- *--split_tot* depends on the dimensionality of input, if 0 a single matrix is loaded

```sh
./correlation_features_run.R \
	--inputFile \
	--sampleAnnFile \
	--tissue_name \
	--split_tot (default 0) \
	--type_data \
	--outFold
```
The output includes (saved in *--outFold*):
- correlation_estimate_<type_data>.RData, R object with correlation matrix and samples info 

#### Direct MR
Mendelian randomization, endophenotypes considered as exposure and trait as outcome
- *--MR_pheno_file* file including UKBB phenotype categories to test
- *--pheno_file* .txt file with TWAS summary statistic or 2 .txt files (Reactome and GO) from PALAS for the trait
- *--cor_file* output from the previous command
- *--enrich_file* .RData output of correlation_pheno_relatedPheno_run.R

```sh
./mendelianRand_pheno_relatedPheno_run.R \
	--MR_pheno_file (default NULL) \
	--tissue_name \
	--pheno_name_comp \
	--inputFold_rel \
	--pheno_file \
	--outFold \
	--pval_corr (default 0.05) \
	--pval_FDR_rel (default 0.05) \
	--type_data \
	--cor_file \
	--enrich_file
```
The output includes (saved in *--outFold*):
- Mendelian_randomization_Egger_<type_data>_pvalFDRrel<pval_FDR_rel>.txt output of Egger method
- Mendelian_randomization_IVW_<type_data>_pvalFDRrel<pval_FDR_rel>.txt output of inverse-variance weighted method


#### Reverse MR
Mendelian randomization, endophenotypes considered as outcome and trait as exposure

```sh
./mendelianRand_reverse_pheno_relatedPheno_run.R \
	--MR_pheno_file (default NULL) \
	--tissue_name \
	--pheno_name_comp SCZ \
	--inputFold_rel \
	--pheno_file \
	--outFold \
	--pval_corr (default 0.05) \
	--pval_FDR_rel (default 0.05) \
	--type_data \
	--cor_file \
	--enrich_file
```
The output includes (saved in *--outFold*):
- Mendelian_randomization_reverse_Egger_<type_data>_pvalFDRrel<pval_FDR_rel>.txt output of Egger method
- Mendelian_randomization_reverse_IVW_<type_data>_pvalFDRrel<pval_FDR_rel>.txt output of inverse-variance weighted method

#### Combine results across tissues
Combine results (MR_egger or MR_IVW) across tissues
- *--mrRes_tissue_file* output from the previous command
- *--comb_corr_file* <type_data>_correlation_sign_pval%s_allTissues.txt file from combine correlation script

```sh
./combine_MR_allTissues_run.R \
	--MR_pheno_file \
	--mrRes_tissue_file \
	--tissue_name \
	--type_data \
	--thr_pval (default 0.05) \
	--outFold \
	--comb_corr_file \
	--type_analysis \
	--mr_type 
```
The output includes the following phenotypes x tissues matrices (saved in *--outFold*):
- <type_data>_Mendelian_randomization_<>_estimates_allTissues.txt MR estimates
- <type_data>_Mendelian_randomization_<>_pvalue_allTissues.txt MR p-value
- <type_data>_Mendelian_randomization_<>_FDRpvalue_allTissues.txt MR FDR p-value (tissue-spec corrections)
- <type_data>_Mendelian_randomization_<>_sign_pval<thr_pval>_allTissues.txt MR sign(estimate)*-log10(p-value)
If mr_type == "Egger", it also includes same files with intercept estimate, p-value and FDR pvalue 
