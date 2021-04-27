#folders=("arguello_neutral" "arguello_deleterious" "arguello_fitrescaled_deleterious" "arguello_beneficials_weak" "arguello_beneficials_strong" "li_and_stephan_neutral_bneck")

python convert_ms_to_fasta.py ${folder} outputrepID.ms num_replicates window_size

#where num_replicates = 100
#and window size = 500 / 1000 / 2000
