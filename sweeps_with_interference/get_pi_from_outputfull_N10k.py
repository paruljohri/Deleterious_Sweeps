#This is to get pi for each site from outputfull in SLiM, immediately post-fixation:

import sys
import os
import argparse

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-folder', dest = 'folder', action='store', nargs = 1, type = str, help = 'only the folder name for the evolutionary scenario')
parser.add_argument('-NumReps', dest = 'NumReps', action='store', nargs = 1, type = int, help = 'total number of replicates')
parser.add_argument('-MutnType', dest = 'MutnType', action='store', nargs = 1, type = str, help = 'del/neutral')
args = parser.parse_args()
s_folder = args.folder[0]
num_reps = args.NumReps[0]
mutn_type = args.MutnType[0]

if "Droso" in s_folder:
	s_range = 500 #get  diversity upto 500 bp on both sides of the fixation:
	f0 = 0.25
	win_size = 10
elif "Human" in s_folder:
	s_range = 40000 #10 kb
	f0 = 0.51
	win_size = 200
def get_pi_from_outputfull(f_txt):
	d_PI = {}
	for line in f_txt:
		line1 = line.strip('\n')
		if "#" not in line1:
			if "Mutations" not in line1:
				if "m5" in line1 or "m1" in line1:
					line2 = line1.split()
					if line2[2] == "m5" or line2[2] == "m1":
						posn = line2[3]
						AF = int(line2[8])
						if line2[2] == "m1":
							d_PI[posn] = 2.0*(1.0/f0)*(float(AF)/100.0)*float(float(100-AF)/100.0)
						else:
							d_PI[posn] = 2.0*(float(AF)/100.0)*float(float(100-AF)/100.0)
	return(d_PI)

def get_pi_AF(d_af):
	d_PI = {}
	for x in d_af.keys():
		d_PI[x] = 2.0*(float(d_af[x])/100.0)*float(float(100-d_af[x])/100.0)
	return(d_PI)

def get_pi_window(d_PI, start, end):
	s_pi = 0.0
	s_num = 0
	s_pos = int(start)
	while s_pos <= end:
		if s_pos >= 0 and s_pos <= 63999:
			s_pi = s_pi + d_PI.get(str(s_pos), 0.0)
			s_num += 1
		s_pos += 1
	#print (s_num)
	if s_num > 0:
		return(s_pi/float(s_num))
	else:
		return("NA")


#get the list of fixations:
if mutn_type == "del":
	f_time = open("/home/pjohri/DelSweeps/" + s_folder + "/time_to_fixation_m2.txt", 'r')
elif mutn_type == "neutral":
	f_time = open("/home/pjohri/DelSweeps/" + s_folder + "/time_to_fixation_m5.txt", 'r')
d_subs = {}
d_col = {}
subs_num = 0
for line in f_time:
	line1 = line.strip('\n')
	line2 = line1.split('\t')
	if line2[0] == "repID":
		col = 0
		for x in line2:
			d_col[x] = col
			col += 1
	else:
		subs_num += 1
		d_subs[subs_num] = {}
		d_subs[subs_num]["repID"] = line2[d_col["repID"]]
		d_subs[subs_num]["posn"] = line2[d_col["base_posn"]]
		d_subs[subs_num]["gamma"] = line2[d_col["gamma"]]
		d_subs[subs_num]["gen_of_fixation"] = line2[d_col["gen_of_fixation"]]
f_time.close()

#result file:
result = open("/home/pjohri/DelSweeps/" + s_folder + "/pi_postfixation_" + mutn_type + "_" + str(win_size) + ".txt", 'w+')
result.write("num" + '\t' + "repID" + '\t' + "gamma")
s_pos = -1*(s_range-(win_size/2))
while s_pos < s_range:
	result.write('\t' + str(s_pos))
	s_pos += win_size
result.write('\n')

#read each file and calculate pi
i = 1
while i <= subs_num:
	print ("substitution number: " + str(i))
	f_txt = open("/mnt/storage/pjohri/" + s_folder + "/" + d_subs[i]["repID"] + "/output_postfixation_" + str(d_subs[i]["gen_of_fixation"]) + ".txt", 'r')
	d_pi = get_pi_from_outputfull(f_txt)
	f_txt.close()
	#get pi values for windows on size 10 
	d_pi_win = {}
	startL = int(d_subs[i]["posn"]) - 1
	startR = int(d_subs[i]["posn"]) + 1
	endL = startL - win_size + 1
	endR = startR + win_size - 1
	s_dist = win_size/2
	while s_dist < s_range:
		d_pi_win[s_dist] = get_pi_window(d_pi, startR, endR)
		d_pi_win[-1*s_dist] = get_pi_window(d_pi, endL, startL)
		startR = endR + 1
		endR = startR + win_size - 1
		startL = endL - 1
		endL = startL - win_size + 1
		s_dist = s_dist + win_size
	#write out the pi values:
	result.write(str(i) + '\t' + d_subs[i]["repID"] + '\t' + d_subs[i]["gamma"])
	s_pos = -1*(s_range-(win_size/2))
	while s_pos < s_range:
		result.write('\t' + str(d_pi_win[s_pos]))
		s_pos = s_pos + win_size
	result.write('\n')
	i += 1

print ("done")
