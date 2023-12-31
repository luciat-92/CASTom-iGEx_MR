#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/correlation_totPath_t%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/correlation_totPath_t%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=70G

module load R/3.5.3

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < OUTPUT_GTEx/Tissue_CADgwas
t=$(eval echo "\${tissues[${id_t}-1]}")

fold=OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/
cov_fold=INPUT_DATA_GTEx/CAD/Covariates/UKBB/
git_fold=/psycl/g/mpsziller/lucia/castom-igex_mr/R/

./${git_fold}correlation_features_run.R \
	--inputFile ${fold}Pathway_Reactome_scores.RData \
	${fold}Pathway_GO_scores.RData \
	--sampleAnnFile ${cov_fold}covariateMatrix_forCorrelation.txt \
	--type_data Reactome GO \
	--outFold ${fold} \
	--tissues_name ${t}

