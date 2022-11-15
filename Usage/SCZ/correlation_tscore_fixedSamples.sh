#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/correlation_tscore_t%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/correlation_tscore_t%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=30G

module load R/3.5.3

cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT
id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Tissues_PGC_red2
t=$(eval echo "\${tissues[${id_t}-1]}")

if [[ "${t}" == "DLPC_CMC" ]]
then
fold=/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/
else
fold=/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/
fi

cov_fold=INPUT_DATA_GTEx/CAD/Covariates/UKBB/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

Rscript ${git_fold}correlation_features_run.R \
	--inputFile ${fold}predictedTscores_splitGenes \
	--sampleAnnFile /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/covariateMatrix_forCorrelation.txt \
	--split_tot 100	\
	--type_data tscore \
	--outFold ${fold} \
	--tissues_name ${t}

