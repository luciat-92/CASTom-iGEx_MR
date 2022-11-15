#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/plot_corr_mendRandom_tscore_UKBB.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/plot_corr_mendRandom_tscore_UKBB.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=30G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/castom-igex_mr/Figures/

tissues=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Colon_Sigmoid Colon_Transverse Heart_Atrial_Appendage Heart_Left_Ventricle Liver Whole_Blood)
fold=OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/

Rscript ${git_fold}plot_correlation_MR_pheno_related_run.R \
	--tissue_name ${tissues[@]} \
	--corrRes_file ${fold}dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt \
	--mrRes_Egg_file ${fold}dist200000b_tscore_Mendelian_randomization_Egger_estimates_allTissues.txt \
	--mrRes_Egg_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_allTissues.txt \
	--mrRes_Egg_pval_file ${fold}dist200000b_tscore_Mendelian_randomization_Egger_estimates_pvalue_allTissues.txt \
	--mrRes_Egg_pval_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_pvalue_allTissues.txt \
	--mrRes_Egg_intpval_file ${fold}dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_pvalue_allTissues.txt \
	--mrRes_Egg_intpval_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_pvalue_allTissues.txt \
	--mrRes_IVW_file ${fold}dist200000b_tscore_Mendelian_randomization_IVW_estimates_allTissues.txt \
	--mrRes_IVW_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_allTissues.txt \
	--mrRes_IVW_pval_file ${fold}dist200000b_tscore_Mendelian_randomization_IVW_estimates_pvalue_allTissues.txt \
	--mrRes_IVW_pval_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_pvalue_allTissues.txt \
	--mrRes_IVW_FDRpval_file ${fold}dist200000b_tscore_Mendelian_randomization_IVW_estimates_FDRpvalue_allTissues.txt \
	--mrRes_IVW_FDRpval_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_intFDRpval_file ${fold}dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_intFDRpval_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_FDRpval_file ${fold}dist200000b_tscore_Mendelian_randomization_Egger_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_FDRpval_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_int_file ${fold}dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_allTissues.txt \
	--mrRes_Egg_int_rev_file ${fold}dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_allTissues.txt \
	--color_pheno_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/color_pheno_type_UKBB.txt \
	--color_tissues_file /psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt \
	--pheno_name CAD \
	--outFold ${fold}dist200000b_ \
	--pheno_list_MR_file ${git_fold}../Usage/Data/plot_pheno_MRext_CAD.txt \
	--data_type tscore \



