#We used simulations performed for the FST analyses provided here- https://github.com/paruljohri/Deleterious_Sweeps/tree/main/FST to conduct AMOVA and identify outliers.

#Step 1: Convert .ms files into fasta files
python convert_ms_to_fasta.py ${folder} outputrepID.ms num_replicates window_size

#where folders=("arguello_neutral" "arguello_deleterious" "arguello_fitrescaled_deleterious" "arguello_beneficials_weak" "arguello_beneficials_strong" "li_and_stephan_neutral_bneck")
#num_replicates = 100
#and window size = 500 / 1000 / 2000
#The exact commanndlines used are provided below:
python convert_ms_to_fasta.py neutral_arguello recurrent_sweeps_split_strictly_neutral_repID.ms 100 500
python convert_ms_to_fasta.py del_dfe_fitrescaled recurrent_sweeps_split_arguello_bnz_repID_fitrescaled.ms 100 500
python convert_ms_to_fasta.py del_dfe_fitnotrescaled recurrent_sweeps_split_arguello_bnz_norescale_repID_full_muts.ms 100 500
python convert_ms_to_fasta.py arguello_beneficials_weak outputrepID.ms 100 500
python convert_ms_to_fasta.py arguello_beneficials_strong outputrepID.ms 100 500
python convert_ms_to_fasta.py neutral_bneck_li_stephan recurrent_sweeps_split_li_stephan_repID_neutral_bneck.ms 10 500

#Step 2: Perform Amova:
#Use amova.r to perform that analysis for all folders mentioned above

#Step 3:Get the proportion of mutations in those outlier regions:
python get_proportion_mutn_types_amova_outliers.py 500
python get_proportion_mutn_types_amova_outliers.py 1000
python get_proportion_mutn_types_amova_outliers.py 2000

#Step 4: Get the proportion of mutation types in all segregating mutations:
python get_proportion_mutn_types_amova_segregating.py
