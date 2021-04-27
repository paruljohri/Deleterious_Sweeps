#This is to covert a .ms file into a fasta file:
#How to run:
#python convert_ms_to_fasta.py neutral_arguello recurrent_sweeps_split_strictly_neutral

import sys
import os

s_folder = sys.argv[1] #neutral_arguello
s_prefix = sys.argv[2] # e.g. recurrent_sweeps_split_strictly_neutral
num_reps = int(sys.argv[3])
win_size = int(sys.argv[4]) #100/500/1000/2000

region_size = 10000

os.system("mkdir /scratch/pjohri1/DelSweeps/AMOVA/" + s_folder + "_FASTA_" + str(win_size) + "bp")

#read ms file:
def read_subset_ms(f_ms, start, end, region_len):
    l_Pos = [] #list of positions of SNPs
    d_Gt = {}
    l_int_posn = []
    for line in f_ms:
        line1 = line.strip('\n')
        if "positions" in line1:
            line2 = line1.split()
            i = 0
            for x in line2:
                if "position" not in x:
                    x_int = round(float(x) * (region_len-1))
                    if (x_int >= float(start)) and (x_int <= float(end)):
                        l_Pos.append(float(x))
                        d_Gt[str(i)] = ""
                        l_int_posn.append(i)
                    i = i + 1
        elif "//" not in line and "segsites" not in line:
            for j in l_int_posn:
                d_Gt[str(j)] = d_Gt[str(j)] + line1[j]
    return (d_Gt, l_int_posn)

#define 500 bp fixed windows:
d_windows, d_start, d_end = {}, {}, {}
win = 1
l_windows = []
posn = 0
while posn < region_size-1:
    l_windows.append(str(win))
    d_start[str(win)] = posn
    d_end[str(win)] = posn + win_size - 1
    posn = posn + win_size
    win = win + 1
print (l_windows)
print (d_start)
print (d_end)
#Read each replicate file and write out each window as a separate fasta file
repID = 1
while repID <= num_reps:
    print ("rep: " + str(repID))
    f_ms = open("/scratch/pjohri1/DelSweeps/FST/" + s_folder + "/" + s_prefix.replace("repID", str(repID)), 'r')
    for win in l_windows:
        t_msi = read_subset_ms(f_ms, d_start[win], d_end[win], region_size)
        d_gti = t_msi[0]
        l_posni = t_msi[1]
        #print(l_posni)
        #print(d_gti)
        f_ms.seek(0)
        result = open("/scratch/pjohri1/DelSweeps/AMOVA/" + s_folder + "_FASTA_" + str(win_size) + "bp/rep" + str(repID) + "_win" + str(win) + ".fasta", 'w+')
        num_indv = 1
        while num_indv <= len(d_gti[str(l_posni[0])]):
            result.write(">" + str(num_indv) + '\n')
            for posn in l_posni:
                allele = d_gti[str(posn)][num_indv - 1]
                if allele == "0":
                    result.write("A")
                elif allele == "1":
                    result.write("T")
                else:
                    print ("Something wrong with the ms file 0 and 1s")
            result.write('\n')
            num_indv += 1
        result.close()
    f_ms.close()
    repID += 1

print ("done")
