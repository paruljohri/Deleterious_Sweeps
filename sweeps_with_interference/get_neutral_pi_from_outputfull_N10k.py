#This is to get pi for each site from outputfull in SLiM, immediately post-fixation:

import sys
import os
import argparse
import numpy as np

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-folder', dest = 'folder', action='store', nargs = 1, type = str, help = 'only the folder name for the evolutionary scenario')
parser.add_argument('-NumReps', dest = 'NumReps', action='store', nargs = 1, type = int, help = 'total number of replicates')
args = parser.parse_args()
s_folder = args.folder[0]
num_reps = args.NumReps[0]

N = 10000
last_gen = 100*N
if "Human" in s_folder:
	f0 = 0.51
elif "Droso" in s_folder:
	f0 = 0.25

def get_pi_from_outputfull(f_txt, mutn_type):
	d_PI = {}
	for line in f_txt:
		line1 = line.strip('\n')
		if "#" not in line1:
			if "Mutations" not in line1:
				if mutn_type in line1:
					line2 = line1.split()
					if line2[2] == mutn_type:
						posn = line2[3]
						AF = int(line2[8])
						if line2[2] == "m1":
							d_PI[posn] = 2.0*(1.0/f0)*(float(AF)/100.0)*float(float(100-AF)/100.0)
						else:
							d_PI[posn] = 2.0*(float(AF)/100.0)*float(float(100-AF)/100.0)
	return(d_PI)



def get_pi_window(d_PI, start, end):
	s_pi = 0.0
	s_pos = int(start)
	while s_pos <= end:
		s_pi = s_pi + d_PI.get(str(s_pos), 0.0)
		s_pos += 1
	return(s_pi/float(win_size))



#result file:
result = open("/home/pjohri/DelSweeps/" + s_folder + "/pi_neutral.txt", 'w+')
result.write("repID" + '\t' + "pi_intergenic" + '\t' + "pi_coding_neutral" + '\n')


#read each file and calculate pi
if "Droso" in s_folder:
	exon_len = 300
	num_genes = 2
	num_exons_per_gene = 5
	inter_len = 4000
	intron_len = 100
	num_introns_per_gene = 4
elif "Human" in s_folder:
	exon_len = 300
	num_genes = 2
	num_exons_per_gene = 5
	inter_len = 15000
	intron_len = 2000
	num_introns_per_gene = 4
coding_len = exon_len*num_genes*num_exons_per_gene
inter_len = (inter_len*(num_genes+1)) + (intron_len*num_genes*num_introns_per_gene)
repID = 1
while repID <= num_reps:
	print ("rep number: " + str(repID))
	f_txt = open("/mnt/storage/pjohri/" + s_folder + "/rep" + str(repID) + "/output_gen" + str(last_gen) + ".txt", 'r')
	d_pi_inter = get_pi_from_outputfull(f_txt, "m5")
	f_txt.seek(0)
	d_pi_syn = get_pi_from_outputfull(f_txt, "m1")
	f_txt.close()
	#get pi values
	s_pi_inter, s_pi_syn = 0.0 , 0.0
	for x in d_pi_inter.keys():
		s_pi_inter = s_pi_inter + d_pi_inter[x]
	for y in d_pi_syn.keys():
		s_pi_syn = s_pi_syn + d_pi_syn[y]
	#write out the pi values:
	result.write(str(repID) + '\t' + str(s_pi_inter/inter_len) + '\t' + str(s_pi_syn/coding_len) + '\n')
	repID += 1

result.close()
print ("done")

