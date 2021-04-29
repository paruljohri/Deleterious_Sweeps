#This is to summarize the final values of single sweeps.
#Converting base pairs to rho:
#python get_final_stats_single_sweeps.py gamma_pos_five 5

import sys

folder = sys.argv[1] #gamma_pos_five
bases = int(sys.argv[2])
if "pos" in folder:
    s_col = "blue"
elif "neg" in folder:
    s_col = "red"
else:
    s_col = "black"

pi_neutral = 0.012
Ne = 1000000.0
rec = 0.00000001
region_size = 10000
mid_size = region_size/2
start = -1.0*round((float(mid_size) - (float(bases)/2.0)),1)

#reading values:
d_pi = {}
d_td = {}
f_win = open("/scratch/pjohri1/DelSweeps/single_sweeps/" + folder + "_" + str(bases) + "_SE.winsummary", 'r')
for line in f_win:
    line1 = line.strip('\n')
    line2 = line1.split('\t')
    if "thetapi_m" in line2[0]:
        d_pi["m"] = line2[1:]
    elif "thetapi_sd" in line2[0]:
        d_pi["sd"] = line2[1:]
    elif "thetapi_se" in line2[0]:
        d_pi["se"] = line2[1:]
    elif "tajimasd_m" in line2[0]:
        d_td["m"] = line2[1:]
    elif "tajimasd_sd" in line2[0]:
        d_td["sd"] = line2[1:]
    elif "tajimasd_se" in line2[0]:
        d_td["se"] = line2[1:]
f_win.close()

#writing it out:
result = open("/home/pjohri1/DelSweeps/single_sweeps/" + folder + "_" + str(bases) + ".final", 'w+')
result.write("bases" + '\t' + "rho" + '\t' + "thetapi_m" + '\t' + "thetapi_sd" + '\t' + "thetapi_se" + '\t' + "tajimasd_m" + '\t' + "tajimasd_sd" + '\t' + "tajimasd_se" + '\t' + "color" + '\n')
posn = start
i = 0
while posn <= mid_size:
    result.write(str(posn) + '\t' + str(2.0*Ne*rec*posn) + '\t' + str(float(d_pi["m"][i])/pi_neutral) + '\t' + str(float(d_pi["sd"][i])/pi_neutral) + '\t' + str(float(d_pi["se"][i])/pi_neutral) + '\t' + str(d_td["m"][i]) + '\t' + str(d_td["sd"][i]) + '\t' + str(d_td["se"][i]) + '\t' + s_col + '\n')
    i = i + 1
    posn = posn + bases
result.close()

print ("Finished")


