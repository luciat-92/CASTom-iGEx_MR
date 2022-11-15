#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/combine_MR_tscore_SCZ_SCZrel_dist200kb.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/combine_MR_tscore_SCZ_SCZrel_dist200kb.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=20G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Figures/model_prediction/
readarray -t tissues < /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Tissues_PGC_red2

fold=Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/

mrRes_Egg=()
mrRes_rev_Egg=()
mrRes_IVW=()
mrRes_rev_IVW=()

for t in ${tissues[@]}
do
	mrRes_Egg+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_dist200000b_Mendelian_randomization_Egger_tscore_pvalFDRrel0.05.txt)
	mrRes_rev_Egg+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_dist200000b_Mendelian_randomization_reverse_Egger_tscore_pvalFDRrel0.05.txt)
	mrRes_IVW+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_dist200000b_Mendelian_randomization_IVW_tscore_pvalFDRrel0.05.txt)
	mrRes_rev_IVW+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_dist200000b_Mendelian_randomization_reverse_IVW_tscore_pvalFDRrel0.05.txt)
done

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
	--mrRes_tissue_file ${mrRes_Egg[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tscore \
	--thr_pval 0.05 \
	--outFold ${fold}matchUKBB_dist200000b_ \
	--comb_corr_file ${fold}matchUKBB_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt \
	--type_analysis not_reverse \
	--mr_type Egger 

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
        --mrRes_tissue_file ${mrRes_rev_Egg[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tscore \
        --thr_pval 0.05 \
        --outFold ${fold}matchUKBB_dist200000b_ \
        --comb_corr_file ${fold}matchUKBB_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis	reverse \
        --mr_type Egger	

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
        --mrRes_tissue_file ${mrRes_IVW[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tscore \
        --thr_pval 0.05 \
        --outFold ${fold}matchUKBB_dist200000b_ \
        --comb_corr_file ${fold}matchUKBB_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis	not_reverse \
        --mr_type IVW	

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file subset_pheno_MR.txt \
        --mrRes_tissue_file ${mrRes_rev_IVW[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tscore \
        --thr_pval 0.05 \
        --outFold ${fold}matchUKBB_dist200000b_ \
        --comb_corr_file ${fold}matchUKBB_dist200000b_tscore_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis reverse \
        --mr_type IVW


