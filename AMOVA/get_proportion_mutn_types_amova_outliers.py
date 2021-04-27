#This is to get the proportion of mutations of each type in the outliers of AMOVA analysis:

import sys
win_size = int(sys.argv[1])
region_size = 10000

l_folders = ["neutral_arguello", "del_dfe_fitnotrescaled", "del_dfe_fitrescaled", "arguello_beneficials_weak", "arguello_beneficials_strong", "neutral_bneck_li_stephan"]
l_prefix = ["recurrent_sweeps_split_strictly_neutral_repID.full", "recurrent_sweeps_split_arguello_bnz_norescale_repID_full_muts.full", "recurrent_sweeps_split_arguello_bnz_repID_fitrescaled.full", "outputrepID.full", "outputrepID.full", "recurrent_sweeps_split_li_stephan_repID_neutral_bneck.full"]

result = open("/home/pjohri1/DelSweeps/AMOVA/outliers_1percent_" + str(win_size) + ".mutntypes", 'w+')
result.write("folder" + '\t' + "m1" + '\t' + "m2" + '\t' + "m3" + '\t' + "m4" + '\t' + "m5" + '\n')

#define the windows
d_window = {}
win = 1
posn = 0
while posn < region_size-1:
    #l_windows.append(str(win))
    s_start = posn
    s_end = posn + win_size - 1
    i = s_start
    while i <= s_end:
        d_window[str(i)] = str(win)
        i += 1
    posn = posn + win_size
    win = win + 1
d_window[str(10000)] = str(win-1)

#Look through the outliers:
m = 0
for s_folder in l_folders:
    l_mutn_types = []
    l_regions = []
    f_out = open("/home/pjohri1/DelSweeps/AMOVA/" + s_folder + "_outliers_1percent_" + str(win_size) + ".phiST", 'r') 
    for line in f_out:
        line1 = line.strip('\n')
        line2 = line1.split('\t')
        if "replicate" not in line2[0]:
            l_regions.append(line2[0] + '\t' + line2[1])
    f_out.close()
    
    #read each of those files and those regions and store the number of mutations
    for region in l_regions:
        print ("outlier region: " + region)
        repID = region.split('\t')[0]
        windowID = region.split('\t')[1]
        f_txt = open("/scratch/pjohri1/DelSweeps/FST/" + s_folder + "/" + l_prefix[m].replace("repID", str(repID)), 'r')
        for line in f_txt:
            line1 = line.strip('\n')
            if "#" not in line1:
                if "Mutations" not in line1:
                    if "m" in line1:
                        line2 = line1.split()
                        if len(line2) == 9:
                            posn = line2[3]
                            if d_window[posn] == windowID:
                                mutn_type = line2[2]
                                l_mutn_types.append(mutn_type)
                                print("position, window, mutntype: " + str(posn) + ", " + d_window[posn] + ", " + mutn_type)
        f_txt.close()

    #write these out
    result.write(s_folder + '\t' + str(l_mutn_types.count("m1")) + '\t' + str(l_mutn_types.count("m2")) + '\t' + str(l_mutn_types.count("m3")) + '\t' + str(l_mutn_types.count("m4")) + '\t' + str(l_mutn_types.count("m5")) + '\n')

    m = m + 1

result.close()
print ("done")

