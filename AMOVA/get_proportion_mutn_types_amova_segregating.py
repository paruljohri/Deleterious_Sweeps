#This is to get the proportion of mutations of each type in the outliers of AMOVA analysis:

import sys

l_folders = ["neutral_arguello", "del_dfe_fitnotrescaled", "del_dfe_fitrescaled", "arguello_beneficials_weak", "arguello_beneficials_strong", "neutral_bneck_li_stephan"]
l_prefix = ["recurrent_sweeps_split_strictly_neutral_repID.full", "recurrent_sweeps_split_arguello_bnz_norescale_repID_full_muts.full", "recurrent_sweeps_split_arguello_bnz_repID_fitrescaled.full", "outputrepID.full", "outputrepID.full", "recurrent_sweeps_split_li_stephan_repID_neutral_bneck.full"]
l_num_reps = [100, 100, 100, 20, 20, 10]

result = open("/home/pjohri1/DelSweeps/AMOVA/segregating.mutntypes", 'w+')
result.write("folder" + '\t' + "m1" + '\t' + "m2" + '\t' + "m3" + '\t' + "m4" + '\t' + "m5" + '\n')

#Look through the outliers:
m = 0
for s_folder in l_folders:
    l_mutn_types = []
    
    #go through replicates:
    repID = 1
    while repID <= l_num_reps[m]:
        print ("repID: " + str(repID))
        f_txt = open("/scratch/pjohri1/DelSweeps/FST/" + s_folder + "/" + l_prefix[m].replace("repID", str(repID)), 'r')
        for line in f_txt:
            line1 = line.strip('\n')
            if "#" not in line1:
                if "Mutations" not in line1:
                    if "m" in line1:
                        line2 = line1.split()
                        if len(line2) == 9:
                            mutn_type = line2[2]
                            l_mutn_types.append(mutn_type)
        f_txt.close()
        repID = repID + 1
    #write these out
    result.write(s_folder + '\t' + str(l_mutn_types.count("m1")) + '\t' + str(l_mutn_types.count("m2")) + '\t' + str(l_mutn_types.count("m3")) + '\t' + str(l_mutn_types.count("m4")) + '\t' + str(l_mutn_types.count("m5")) + '\n')

    m = m + 1

result.close()
print ("done")

