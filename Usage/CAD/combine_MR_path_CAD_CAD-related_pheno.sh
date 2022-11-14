#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/combine_MR_path_CAD_CADrel_perc0.3.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/combine_MR_path_CAD_CADrel_perc0.3.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=20G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Figures/model_prediction/
tissues=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Colon_Sigmoid Colon_Transverse Heart_Atrial_Appendage Heart_Left_Ventricle Liver Whole_Blood)
fold=OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/

mrRes_Egg=()
mrRes_rev_Egg=()
mrRes_IVW=()
mrRes_rev_IVW=()

for t in ${tissues[@]}
do
	mrRes_Egg+=(OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_perc0.3_Mendelian_randomization_Egger_tot_path_pvalFDRrel0.05.txt)
	mrRes_rev_Egg+=(OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_perc0.3_Mendelian_randomization_reverse_Egger_tot_path_pvalFDRrel0.05.txt)
	mrRes_IVW+=(OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_perc0.3_Mendelian_randomization_IVW_tot_path_pvalFDRrel0.05.txt)
	mrRes_rev_IVW+=(OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/CADCardioG_perc0.3_Mendelian_randomization_reverse_IVW_tot_path_pvalFDRrel0.05.txt)
done

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/MR_subset_pheno.txt \
	--mrRes_tissue_file ${mrRes_Egg[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tot_path \
	--thr_pval 0.05 \
	--outFold ${fold}CADCardioG_perc0.3_ \
	--comb_corr_file ${fold}CADCardioG_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
	--type_analysis not_reverse \
	--mr_type Egger 

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/MR_subset_pheno.txt \
        --mrRes_tissue_file ${mrRes_rev_Egg[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tot_path \
        --thr_pval 0.05 \
        --outFold ${fold}CADCardioG_perc0.3_ \
        --comb_corr_file ${fold}CADCardioG_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis	reverse \
        --mr_type Egger	

Rscript ${git_fold}combine_MR_allTissues_run.R \
	--MR_pheno_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/MR_subset_pheno.txt \
        --mrRes_tissue_file ${mrRes_IVW[@]} \
        --tissue_name ${tissues[@]} \
        --type_data tot_path \
        --thr_pval 0.05 \
        --outFold ${fold}CADCardioG_perc0.3_ \
        --comb_corr_file ${fold}CADCardioG_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis	not_reverse \
        --mr_type IVW	

Rscript ${git_fold}combine_MR_allTissues_run.R \
        --mrRes_tissue_file ${mrRes_rev_IVW[@]} \
	--MR_pheno_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/MR_subset_pheno.txt \
        --tissue_name ${tissues[@]} \
        --type_data tot_path \
        --thr_pval 0.05 \
        --outFold ${fold}CADCardioG_perc0.3_ \
        --comb_corr_file ${fold}CADCardioG_perc0.3_tot_path_correlation_sign_pval0.05_allTissues.txt \
        --type_analysis reverse \
        --mr_type IVW


