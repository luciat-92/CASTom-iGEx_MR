##############################################################################################
# enrichment CAD realted phenotypes
sbatch --job-name=enr --array=1-5,7-11,15 Enrichment_CAD_CAD-related_pheno_pertissue.sh
sbatch --job-name=enr_tot Enrichment_CAD_CAD-related_pheno_alltissues.sh

# correlation CAD related phenotypes (from CardioGram and UKBB)
sbatch --job-name=corr --array=1-5,7-11,15 correlation_CAD_CAD-related_pheno_pertissue.sh
sbatch --job-name=corr --array=1-5,7-11,15 correlation_CADUKBB_CAD-related_pheno_pertissue.sh
# combine results
sbatch --job-name=combine combine_correlation_CAD_CAD-related_pheno.sh
sbatch --job-name=combine combine_correlation_CADUKBB_CAD-related_pheno.sh

#sbatch --job-name=corr_tot correlation_CAD_CAD-related_pheno_alltissues.sh
#sbatch --job-name=corr_plot plot_correlation_CAD_CAD-related_pheno.sh

##############################################################################################
# mendelian randomization
sbatch --job-name=subset_corr create_samples_forcorr.sh
sbatch --job-name=cor_feat --array=1-5,7-11,15 correlation_tscore_fixedSamples.sh
sbatch --job-name=cor_feat --array=1-5,7-11,15 correlation_path_fixedSamples.sh

sbatch --job-name=mr_tscore --array=1-5,7-11,15 mendelianRandom_tscore_CAD_CAD-related_pertissue.sh
sbatch --job-name=perc0.3 --array=1-5,7-11,15 mendelianRandom_path_CAD_CAD-related_pertissue.sh 0.3
sbatch --job-name=mr_tscore --array=1-5,7-11,15 mendelianRandom_tscore_CADUKBB_CAD-related_pertissue.sh
sbatch --job-name=perc0.3 --array=1-5,7-11,15 mendelianRandom_path_CADUKBB_CAD-related_pertissue.sh 0.3
# sbatch --job-name=perc0.7 --array=1-5,7-11,15 mendelianRandom_path_CAD_CAD-related_pertissue.sh	0.7

sbatch --job-name=mr_tscore --array=1-5,7-11,15 mendelianRandom_reverse_tscore_CAD_CAD-related_pertissue.sh
sbatch --job-name=perc0.3 --array=1-5,7-11,15 mendelianRandom_reverse_path_CAD_CAD-related_pertissue.sh 0.3
sbatch --job-name=mr_tscore --array=1-5,7-11,15 mendelianRandom_reverse_tscore_CADUKBB_CAD-related_pertissue.sh
sbatch --job-name=perc0.3 --array=1-5,7-11,15 mendelianRandom_reverse_path_CADUKBB_CAD-related_pertissue.sh 0.3

# combine results
sbatch --job-name=combine combine_MR_tscore_CAD_CAD-related_pheno.sh
sbatch --job-name=combine combine_MR_path_CAD_CAD-related_pheno.sh
sbatch --job-name=combine combine_MR_tscore_CADUKBB_CAD-related_pheno.sh
sbatch --job-name=combine combine_MR_path_CADUKBB_CAD-related_pheno.sh

#sbatch --job-name=plot_CAD plot_MR_corr_OR_allTissues.sh
sbatch --job-name=pl_tscore plot_MR_corr_tscore_allTissues.sh
sbatch --job-name=pl_path plot_MR_corr_path_allTissues.sh
sbatch --job-name=pl_tscore plot_MR_corr_tscore_CADUKBB_allTissues.sh
sbatch --job-name=pl_path plot_MR_corr_path_CADUKBB_allTissues.sh

# specific MR plot
sbatch --job-name=pl plot_MRreg_tscore_spec.sh
sbatch --job-name=pl plot_MRreg_path_spec.sh # use this on paper figure
sbatch --job-name=pl plot_MRreg_tscore_BFP.sh
sbatch --job-name=pl plot_MRreg_path_BFP.sh

