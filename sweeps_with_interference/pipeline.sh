#Step 1: run the sims:
bash queue_sims.sh #modify the folder name in the script

#Step 2: Get neutral pi:
python get_neutral_pi_from_outputfull_N10k.py -folder Human_N10k_mean_rec -NumReps 600
python get_neutral_pi_from_outputfull_N10k.py -folder Human_N10k_half_rec -NumReps 600
python get_neutral_pi_from_outputfull_N10k.py -folder Human_N10k_tenth_rec -NumReps 600
python get_neutral_pi_from_outputfull_N10k.py -folder Droso_N10k_mean_rec -NumReps 100
python get_neutral_pi_from_outputfull_N10k.py -folder Droso_N10k_half_rec -NumReps 100
python get_neutral_pi_from_outputfull_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100
#Step 3: Add values of neutral pi into the scripts below:

#Step 4:
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

#Summarize prob and time of fixation:
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


#Get pi around sweeps:
python get_pi_from_outputfull_N10k.py -folder Droso_N10k_mean_rec -NumReps 100 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Droso_N10k_half_rec -NumReps 100 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Droso_N10k_tenth_rec -NumReps 100 -MutnType del

python get_pi_from_outputfull_N10k.py -folder Human_N10k_mean_rec -NumReps 600 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Human_N10k_half_rec -NumReps 600 -MutnType del
python get_pi_from_outputfull_N10k.py -folder Human_N10k_tenth_rec -NumReps 600 -MutnType del




