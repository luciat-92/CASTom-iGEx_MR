#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/combine_correlation_SCZ_SCZrel_perc0.3_dist200kb.out
#SBATCH -e /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/err_out_fold/combine_correlation_SCZ_SCZrel_perc0.3_dist200kb.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=20G


module load R/3.5.3

cd /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/castom-igex_mr/R/

readarray -t tissues < /psycl/g/mpsziller/lucia/SCZ_PGC/eQTL_PROJECT/Meta_Analysis_SCZ/Tissues_PGC_red2

fold=Meta_Analysis_SCZ/OUTPUT_all/enrichment_SCZ-UKBB_res/

corrRes_tissue_file=()
for t in ${tissues[@]}
do
	corrRes_tissue_file+=(Meta_Analysis_SCZ/${t}/enrichment_SCZ-UKBB_res/matchUKBB_perc0.3_dist200000b_correlation_enrich_SCZ_relatedPheno.RData)
done

Rscript ${git_fold}combine_correlation_allTissues_run.R \
	--corrRes_tissue_file ${corrRes_tissue_file[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tscore \
	--thr_pval 0.05 \
	--outFold ${fold}matchUKBB_dist200000b_

Rscript ${git_fold}combine_correlation_allTissues_run.R \
	--corrRes_tissue_file ${corrRes_tissue_file[@]} \
	--tissue_name ${tissues[@]} \
	--type_data tot_path \
	--thr_pval 0.05 \
	--outFold ${fold}matchUKBB_perc0.3_

