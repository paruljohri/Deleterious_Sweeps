#Steps to simulate deleterious sweeps and get the effects on neutral diversity at linked sites:
#All relevant scripts are provided in this folder.

#Step 1: Simulate sweeps using SLiM:
slim -d gamma=5 -d replicate=<enter replicate number here> -d dominance=<dominance coefficient> single_sweeps.slim
slim -d gamma=-5 -d replicate=<enter replicate number here> -d dominance=<dominance coefficient> single_sweeps.slim
slim -d gamma=0 -d replicate=<enter replicate number here> -d dominance=<dominance coefficient> single_sweeps.slim

#Step 2: Get diversity and other stats in sliding windows near the target site:
python statistics_slidingwindow_pylibseq_general_reps.py -winSize 10 -stepSize 10 -regionLen 10000 -folder /path/to/output/folder/gamma_pos_five -output_prefix gamma_pos_five
python statistics_slidingwindow_pylibseq_general_reps.py -winSize 10 -stepSize 10 -regionLen 10000 -folder /path/to/output/folder/gamma_neg_five -output_prefix gamma_neg_five
python statistics_slidingwindow_pylibseq_general_reps.py -winSize 10 -stepSize 10 -regionLen 10000 -folder /path/to/output/folder/gamma_zero -output_prefix gamma_zero

#Step 3: Summarize diversity at linked sites in a tabular form good for plotting:
Rscript ./get_winsummary_general_SEs.R /path/to/working/directory gamma_pos_five 10000 10
python get_final_stats_single_sweeps.py gamma_pos_five 10

Rscript ./get_winsummary_general_SEs.R /path/to/working/directory gamma_neg_five 10000 10
python get_final_stats_single_sweeps.py gamma_neg_five 10

Rscript ./get_winsummary_general_SEs.R /path/to/working/directory gamma_zero 10000 10
python get_final_stats_single_sweeps.py gamma_zero 10

#where 10000 is the size of the simulated region and 10 is the sliding window size

#Step 4: Plot it
#Use the script - plot_single_sweeps.r to plot.
