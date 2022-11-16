#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/mendRandom_rev_path_SCZ_UKBBrel_t%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/mendRandom_rev_path_SCZ_UKBBrel_t%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=30G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

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

Rscript ${git_fold}mendelianRand_reverse_pheno_relatedPheno_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
	--tissue_name ${t} \
	--pheno_name_comp SCZ \
	--inputFold_rel ${fold_UKBB} \
	--pheno_file ${fold_out}path_Reactome_pval_SCZ_covCorr_filt.txt ${fold_out}path_GO_pval_SCZ_covCorr_filt.txt \
	--outFold ${fold}matchUKBB_perc0.3_ \
	--pval_corr 0.05 \
	--pval_FDR_pheno 0.05 \
	--type_data tot_path \
	--cor_file ${fold_UKBB}correlation_estimate_tot_path.RData \
	--enrich_file ${fold}matchUKBB_perc0.3_dist200000b_correlation_enrich_SCZ_relatedPheno.RData
