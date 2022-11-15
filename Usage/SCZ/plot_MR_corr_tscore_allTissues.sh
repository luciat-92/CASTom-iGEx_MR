#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/plot_corr_mendRandom_tscore.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/plot_corr_mendRandom_tscore.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=30G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/castom-igex_mr/Figures/

readarray -t tissues < /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Tissues_PGC_red2
fold=Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/

Rscript ${git_fold}plot_correlation_MR_pheno_related_run.R \
	--tissue_name ${tissues[@]} \
	--corrRes_file ${fold}matchUKBB_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt \
	--mrRes_Egg_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_Egger_estimates_allTissues.txt \
	--mrRes_Egg_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_allTissues.txt \
	--mrRes_Egg_pval_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_Egger_estimates_pvalue_allTissues.txt \
	--mrRes_Egg_pval_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_pvalue_allTissues.txt \
	--mrRes_Egg_intpval_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_pvalue_allTissues.txt \
	--mrRes_Egg_intpval_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_pvalue_allTissues.txt \
	--mrRes_IVW_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_IVW_estimates_allTissues.txt \
	--mrRes_IVW_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_allTissues.txt \
	--mrRes_IVW_pval_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_IVW_estimates_pvalue_allTissues.txt \
	--mrRes_IVW_pval_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_pvalue_allTissues.txt \
	--mrRes_IVW_FDRpval_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_IVW_estimates_FDRpvalue_allTissues.txt \
	--mrRes_IVW_FDRpval_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_IVW_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_intFDRpval_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_intFDRpval_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_FDRpval_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_Egger_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_FDRpval_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_Egger_estimates_FDRpvalue_allTissues.txt \
	--mrRes_Egg_int_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_Egger_inter_estimates_allTissues.txt \
	--mrRes_Egg_int_rev_file ${fold}matchUKBB_dist200000b_tscore_Mendelian_randomization_reverse_Egger_inter_estimates_allTissues.txt \
	--color_pheno_file /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/INPUT_DATA/Covariates/color_pheno_type_UKBB.txt \
	--color_tissues_file /psycl/g/mpsziller/lucia/priler_project/Figures/color_tissues.txt \
	--pheno_name SCZ \
	--outFold ${fold}matchUKBB_dist200000b_ \
	--pheno_list_MR_file ${git_fold}../Usage/Data/plot_pheno_MRext_SCZ.txt \
	--data_type tscore \












