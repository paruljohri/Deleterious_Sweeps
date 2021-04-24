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

