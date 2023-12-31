#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/combine_MR_path_SCZ_SCZrel_perc0.3.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/combine_MR_path_SCZ_SCZrel_perc0.3.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=20G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/castom-igex_mr/R/

readarray -t tissues < /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Tissues_PGC_red2

fold=Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/

mrRes_Egg=()
mrRes_rev_Egg=()
mrRes_IVW=()
mrRes_rev_IVW=()

for t in ${tissues[@]}
do
	mrRes_Egg+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_perc0.3_Mendelian_randomization_Egger_tot_path_pvalFDRrel0.05.txt)
	mrRes_rev_Egg+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_perc0.3_Mendelian_randomization_reverse_Egger_tot_path_pvalFDRrel0.05.txt)
	mrRes_IVW+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_perc0.3_Mendelian_randomization_IVW_tot_path_pvalFDRrel0.05.txt)
	mrRes_rev_IVW+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_perc0.3_Mendelian_randomization_reverse_IVW_tot_path_pvalFDRrel0.05.txt)
done

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
	--mrRes_tissue_file ${mrRes_Egg[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tot_path \
	--thr_pval 0.05 \
	--outFold ${fold}matchUKBB_perc0.3_ \
	--comb_corr_file ${fold}matchUKBB_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
	--type_analysis not_reverse \
	--mr_type Egger 

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
        --mrRes_tissue_file ${mrRes_rev_Egg[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tot_path \
        --thr_pval 0.05 \
        --outFold ${fold}matchUKBB_perc0.3_ \
        --comb_corr_file ${fold}matchUKBB_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis	reverse \
        --mr_type Egger	

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
        --mrRes_tissue_file ${mrRes_IVW[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tot_path \
        --thr_pval 0.05 \
        --outFold ${fold}matchUKBB_perc0.3_ \
        --comb_corr_file ${fold}matchUKBB_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis	not_reverse \
        --mr_type IVW	

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
        --mrRes_tissue_file ${mrRes_rev_IVW[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tot_path \
        --thr_pval 0.05 \
        --outFold ${fold}matchUKBB_perc0.3_ \
        --comb_corr_file ${fold}matchUKBB_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis reverse \
        --mr_type IVW


