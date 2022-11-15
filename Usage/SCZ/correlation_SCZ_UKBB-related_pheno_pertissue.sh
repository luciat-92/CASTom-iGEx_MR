#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/correlation_SCZ_UKBBrel_t%a_perc0.3_dist200.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/correlation_SCZ_UKBBrel_t%a_perc0.3_dist200.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=100G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/castom-igex_mr/R/


id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < Meta_Analysis_SCZ/Tissues_PGC_red2
t=$(eval echo "\${tissues[${id_t}-1]}")

if [[ "${t}" == "DLPC_CMC" ]]
then
fold_UKBB=/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/
else
fold_UKBB=/psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/
fi

fold=Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/
fold_out=Meta_Analysis_SCZ/OUTPUT_all/
fold_filt=/psycl/g/mpsziller/lucia/compare_prediction_UKBB_SCZ-PGC/

./${git_fold}correlation_pheno_relatedPheno_run.R \
	--phenoFold /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/ \
	--tissue_name ${t} \
	--pheno_name_comp SCZ \
	--inputFold_rel ${fold_UKBB} \
	--tscore_pheno_file ${fold_out}tscore_pval_SCZ_covCorr.txt \
	--pathR_pheno_file ${fold_out}path_Reactome_pval_SCZ_covCorr_filt.txt \
	--pathGO_pheno_file ${fold_out}path_GO_pval_SCZ_covCorr_filt.txt \
	--outFold ${fold}/matchUKBB_perc0.3_dist200000b_ \
	--pval_FDR_pheno 0.05 \
	--pval_FDR_rel 0.05 \
	--perc_par 0.3 \
	--thr_dist_par 250000 \
	--refFold /psycl/g/mpsziller/lucia/castom-igex/refData/ \
	--feat_filt ${fold_filt}${t}_filter_genes_matched_datasets.txt ${fold_filt}${t}_filter_path_Reactome_matched_datasets.txt ${fold_filt}${t}_filter_path_GO_matched_datasets.txt
