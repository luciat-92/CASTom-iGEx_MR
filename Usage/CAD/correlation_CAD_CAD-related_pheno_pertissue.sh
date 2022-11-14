#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/correlation_CAD_CADrel_t%a_perc0.3_dist200kb.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/correlation_CAD_CADrel_t%a_perc0.3_dist200kb.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=20G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < OUTPUT_GTEx/Tissue_CADgwas
t=$(eval echo "\${tissues[${id_t}-1]}")

fold=OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/

Rscript ${git_fold}correlation_pheno_relatedPheno_run.R \
	--phenoFold INPUT_DATA_GTEx/CAD/Covariates/UKBB/ \
	--tissue_name ${t} \
	--pheno_name_comp CADCardioG \
	--inputFold_rel OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/ \
	--tscore_pheno_file OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/tscore_pval_Dx_pheno_covCorr.txt \
	--pathR_pheno_file OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_Reactome_pval_Dx_covCorr_filt.txt \
	--pathGO_pheno_file OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/path_GO_pval_Dx_covCorr_filt.txt \
	--outFold ${fold}/enrichment_CADHARD_res/CADCardioG_perc0.3_dist200000b_ \
	--pval_FDR_pheno 0.05 \
	--pval_FDR_rel 0.05 \
	--perc_par 0.3 \
	--thr_dist_par 200000 \
	--refFold /psycl/g/mpsziller/lucia/priler_project/refData/ \

