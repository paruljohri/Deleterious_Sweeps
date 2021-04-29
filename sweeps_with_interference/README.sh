#Here is the pipeline and the set of commandlines used to simulate fixations of deleterious mutations in the presence of potential interference from other deleterious mutations and get the effect of such fixations on neutral  diversity at linked sites.
#All required scripts are provided in this folder.

#Step 1: run the simulations:
slim -s ${repID} -d "d_folder='${output_folder}'" eqm_10kb_del_fixation_Droso_N10k_mean_rec.slim
slim -s ${repID} -d "d_folder='${output_folder}'" eqm_10kb_del_fixation_Droso_N10k_half_rec.slim
slim -s ${repID} -d "d_folder='${output_folder}'" eqm_10kb_del_fixation_Droso_N10k_tenth_rec.slim
slim -s ${repID} -d "d_folder='${output_folder}'" eqm_10kb_del_fixation_Human_N10k_mean_rec.slim
slim -s ${repID} -d "d_folder='${output_folder}'" eqm_10kb_del_fixation_Human_N10k_half_rec.slim
slim -s ${repID} -d "d_folder='${output_folder}'" eqm_10kb_del_fixation_Human_N10k_tenth_rec.slim

#Step 2: Get neutral thetapi at linked sites:
python get_neutral_pi_from_outputfull_N10k.py -folder Human_N10k_mean_rec -NumReps 600
python get_neutral_pi_from_outputfull_N10k.py -folder Human_N10k_half_rec -NumReps 600
python get_neutral_pi_from_outputfull_N10k.py -folder Human_N10k_tenth_rec -NumReps 600
python get_neutral_pi_from_outputfull_N10k.py -folder Droso_N10k_mean_rec -NumReps 100
python get_neutral_pi_from_outputfull_N10k.py -folder Droso_N10k_half_rec -NumReps 100
python get_neutral_pi_from_outputfull_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100

#Step 3:
python get_prob_time_fixation_N10k.py -folder Human_N10k_mean_rec -NumReps 600 -mutnType m1
python get_prob_time_fixation_N10k.py -folder Human_N10k_mean_rec -NumReps 600 -mutnType m2
python get_prob_time_fixation_N10k.py -folder Human_N10k_mean_rec -NumReps 600 -mutnType m5

python get_prob_time_fixation_N10k.py -folder Human_N10k_half_rec -NumReps 600 -mutnType m1
python get_prob_time_fixation_N10k.py -folder Human_N10k_half_rec -NumReps 600 -mutnType m2
python get_prob_time_fixation_N10k.py -folder Human_N10k_half_rec -NumReps 600 -mutnType m5

python get_prob_time_fixation_N10k.py -folder Human_N10k_tenth_rec -NumReps 600 -mutnType m1
python get_prob_time_fixation_N10k.py -folder Human_N10k_tenth_rec -NumReps 600 -mutnType m2
python get_prob_time_fixation_N10k.py -folder Human_N10k_tenth_rec -NumReps 600 -mutnType m5

python get_prob_time_fixation_N10k.py -folder Droso_N10k_mean_rec -NumReps 100 -mutnType m1
python get_prob_time_fixation_N10k.py -folder Droso_N10k_mean_rec -NumReps 100 -mutnType m2
python get_prob_time_fixation_N10k.py -folder Droso_N10k_mean_rec -NumReps 100 -mutnType m5

python get_prob_time_fixation_N10k.py -folder Droso_N10k_half_rec -NumReps 100 -mutnType m1
python get_prob_time_fixation_N10k.py -folder Droso_N10k_half_rec -NumReps 100 -mutnType m2
python get_prob_time_fixation_N10k.py -folder Droso_N10k_half_rec -NumReps 100 -mutnType m5

python get_prob_time_fixation_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100 -mutnType m1
python get_prob_time_fixation_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100 -mutnType m2
python get_prob_time_fixation_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100 -mutnType m5

#Step 4: Summarize prob and time of fixation:
Rscript ./get_final_stats.r Droso_N10k_mean_rec
Rscript ./get_final_stats_BGS.r Droso_N10k_mean_rec

Rscript ./get_final_stats.r Droso_N10k_half_rec
Rscript ./get_final_stats_BGS.r Droso_N10k_half_rec

Rscript ./get_final_stats.r Droso_N10k_tenth_rec
Rscript ./get_final_stats_BGS.r Droso_N10k_tenth_rec

Rscript ./get_final_stats.r Human_N10k_mean_rec
Rscript ./get_final_stats_BGS.r Human_N10k_mean_rec

Rscript ./get_final_stats.r Human_N10k_half_rec
Rscript ./get_final_stats_BGS.r Human_N10k_half_rec

Rscript ./get_final_stats.r Human_N10k_tenth_rec
Rscript ./get_final_stats_BGS.r Human_N10k_tenth_rec


#Step 5: Get pi around sweeps:
python get_pi_from_outputfull_N10k.py -folder Droso_N10k_mean_rec -NumReps 100 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Droso_N10k_half_rec -NumReps 100 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100 -MutnType del

python get_pi_from_outputfull_N10k.py -folder Human_N10k_mean_rec -NumReps 600 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Human_N10k_half_rec -NumReps 600 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Human_N10k_tenth_rec -NumReps 600 -MutnType del

#Step 6: Plot the relative diversity around a fixed site:
#Use this script to plot - plot_single_sweeps_del.r


