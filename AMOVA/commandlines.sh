#folders=("arguello_neutral" "arguello_deleterious" "arguello_fitrescaled_deleterious" "arguello_beneficials_weak" "arguello_beneficials_strong" "li_and_stephan_neutral_bneck")

#Step 1: Convert .ms files into fasta files
python convert_ms_to_fasta.py ${folder} outputrepID.ms num_replicates window_size

#where num_replicates = 100
#and window size = 500 / 1000 / 2000

#Step 2: Perform Amova:

#Step 3:Get the proportion of mutations in those outlier regions:
python get_proportion_mutn_types_amova_outliers.py 500
python get_proportion_mutn_types_amova_outliers.py 1000
python get_proportion_mutn_types_amova_outliers.py 2000

#Step 4: Get the proportion of mutation types in all segregating mutations:
python get_proportion_mutn_types_amova_segregating.py
