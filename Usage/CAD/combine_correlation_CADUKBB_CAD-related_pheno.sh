#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/combine_correlation_CADUKBB_CADrel_perc0.3_dist200kb.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/combine_correlation_CADUKBB_CADrel_perc0.3_dist200kb.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=20G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Figures/model_prediction/
tissues=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Colon_Sigmoid Colon_Transverse Heart_Atrial_Appendage Heart_Left_Ventricle Liver Whole_Blood)
fold=OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/enrichment_CADHARD_res/

corrRes_tissue_file=()
for t in ${tissues[@]}
do
	corrRes_tissue_file+=(OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/enrichment_CADHARD_res/perc0.3_dist200000b_correlation_enrich_CAD_HARD_relatedPheno.RData)
done

Rscript ${git_fold}combine_correlation_allTissues_run.R \
	--corrRes_tissue_file ${corrRes_tissue_file[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tscore \
	--thr_pval 0.05 \
	--outFold ${fold}dist200000b_

Rscript ${git_fold}combine_correlation_allTissues_run.R \
	--corrRes_tissue_file ${corrRes_tissue_file[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tot_path \
	--thr_pval 0.05 \
	--outFold ${fold}perc0.3_


