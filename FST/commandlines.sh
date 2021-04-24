#FST simulations:
slim -d "d_ben_s='small'" -d d_ben_p=0.01 -d "replicate='${repID}'" -d "d_folder='/mnt/storage/pjohri/FST/${folder}'" arguello_beneficials_mean.slim
slim -d "d_ben_s='large'" -d d_ben_p=0.01 -d "replicate='${repID}'" -d "d_folder='/mnt/storage/pjohri/FST/${folder}'" arguello_beneficials_mean.slim

#FST simulations performed for extreme values of the confidence intervals:
#This was run for both _maxtime and _mintime scripts.
pop_afr_min=3.02e6
pop_afr_max=4.69e6
pop_eur_min=2.03e5
pop_eur_max=3.89e6

#run slim on combination 1
/mnt/storage/software/slim.3.1/build/slim -d "d_ben_s='large'" -d d_ben_p=0.01 -d d_N_Afr=${pop_afr_min} -d d_N_Eur=${pop_eur_min} -d "replicate='${repID}'" -d "d_folder='/mnt/storage/pjohri/FST/${folder}1'" arguello_beneficials_posterior_maxtime.slim

#run slim on combination 2
/mnt/storage/software/slim.3.1/build/slim -d "d_ben_s='large'" -d d_ben_p=0.01 -d d_N_Afr=${pop_afr_max} -d d_N_Eur=${pop_eur_max} -d "replicate='${repID}'" -d "d_folder='/mnt/storage/pjohri/FST/${folder}2'" arguello_beneficials_posterior_maxtime.slim

#run slim on combination 3
/mnt/storage/software/slim.3.1/build/slim -d "d_ben_s='large'" -d d_ben_p=0.01 -d d_N_Afr=${pop_afr_min} -d d_N_Eur=${pop_eur_max} -d "replicate='${repID}'" -d "d_folder='/mnt/storage/pjohri/FST/${folder}3'" arguello_beneficials_posterior_maxtime.slim

#run slim on combination 4
/mnt/storage/software/slim.3.1/build/slim -d "d_ben_s='large'" -d d_ben_p=0.01 -d d_N_Afr=${pop_afr_max} -d d_N_Eur=${pop_eur_min} -d "replicate='${repID}'" -d "d_folder='/mnt/storage/pjohri/FST/${folder}4'" arguello_beneficials_posterior_maxtime.slim

##########################Statistics###########################################
#To get FST for slinding windows:
python statistics_slidingwindow_pylibseq_general_reps_fst.py -winSize 500 -stepSize 500 -regionLen 10000 -numIndvPop1 100 -numIndvPop2 100 -input_folder /path/to/${foldername} -output_folder /path/to/output/directory -output_prefix ${foldername}

#To make separate .ms files by mutation type:
python separate_mutn_types_fst.py -chromLen 10000 -mutnType m1 -input_folder /path/to/${foldername} -output_folder /path/to/${foldername}_m1
python separate_mutn_types_fst.py -chromLen 10000 -mutnType m2 -input_folder /path/to/${foldername} -output_folder /path/to/${foldername}_m2
python separate_mutn_types_fst.py -chromLen 10000 -mutnType m5 -input_folder /path/to/${foldername} -output_folder /path/to/${foldername}_m5

#where m1: effectively neutral mutations
#      m2: weakly deleterious mutations
#      m5: beneficial mutations
#Note that chromLen here is 1 base pair less than the full length. That is, our chromosome is 10001 bp long.


#To get statistics for each mutation type separately:
python statistics_slidingwindow_pylibseq_general_reps_fst.py -numIndvPop1 100 -numIndvPop2 100 -winSize 1000 -stepSize 1000 -regionLen 10000 -input_folder /scratch/pjohri1/DelSweeps/FST/arguello_beneficials_strong_m1 -output_folder /scratch/pjohri1/DelSweeps/FST -output_prefix arguello_beneficials_strong_m1
python statistics_slidingwindow_pylibseq_general_reps_fst.py -numIndvPop1 100 -numIndvPop2 100 -winSize 1000 -stepSize 1000 -regionLen 10000 -input_folder /scratch/pjohri1/DelSweeps/FST/arguello_beneficials_strong_m2 -output_folder /scratch/pjohri1/DelSweeps/FST -output_prefix arguello_beneficials_strong_m2
python statistics_slidingwindow_pylibseq_general_reps_fst.py -numIndvPop1 100 -numIndvPop2 100 -winSize 1000 -stepSize 1000 -regionLen 10000 -input_folder /scratch/pjohri1/DelSweeps/FST/arguello_beneficials_strong_m5 -output_folder /scratch/pjohri1/DelSweeps/FST -output_prefix arguello_beneficials_strong_m5

