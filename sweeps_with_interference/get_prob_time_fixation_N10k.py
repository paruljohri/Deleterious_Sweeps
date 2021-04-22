#This is to get the probability and time to fixation of a mildly deleterious mutation:

import sys
import argparse
import numpy

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-folder', dest = 'folder', action='store', nargs = 1, type = str, help = 'only the folder name for the evolutionary scenario')
parser.add_argument('-NumReps', dest = 'NumReps', action='store', nargs = 1, type = int, help = 'total number of replicates')
parser.add_argument('-mutnType', dest = 'mutnType', action='store', nargs = 1, type = str, help = 'm1/m2/m5')
args = parser.parse_args()
s_folder = args.folder[0]
num_reps = args.NumReps[0]
mutn_type = args.mutnType[0]

Ne = 10000
gen_burnin = 10*Ne
last_gen = 100*Ne

if "Droso" in s_folder:
	mutn_rate = 5.85e-7
	#if "mean_rec" in s_folder:
	#	B = 0.01899446/float(4*Ne*mutn_rate)
	#elif "half_rec" in s_folder:
	#	B = 0.01851208/float(4*Ne*mutn_rate)
	#elif "tenth_rec" in s_folder:
	#	B = 0.01563094/float(4*Ne*mutn_rate)
elif "Human" in s_folder:
	mutn_rate = 1.2e-8
	#if "mean_rec" in s_folder:
	#	B = 0.0004082986/float(4*Ne*mutn_rate)
	#elif "half_rec" in s_folder:
	#	B = 0.000395343/float(4*Ne*mutn_rate)
	#elif "tenth_rec" in s_folder:
	#	B = 0.0003986758/float(4*Ne*mutn_rate)
#Get B:
f_neu = open("/home/pjohri/DelSweeps/" + s_folder + "/pi_neutral.txt", 'r')
l_pi = []
for line in f_neu:
	line1 = line.strip('\n')
	if "repID" not in line:
		line2 = line1.split('\t')
		l_pi.append(float(line2[1]))
f_neu.close()
B = numpy.mean(l_pi)/float(4*Ne*mutn_rate)

result = open("/home/pjohri/DelSweeps/" + s_folder + "/time_to_fixation_" + mutn_type + ".txt", 'w+')
result.write("repID" + '\t' + "gamma" + '\t' + "gamma_BGS" + '\t' + "time_gen" + '\t' + "time_Ngen" + '\t' + "gen_of_fixation" + '\t' + "base_posn" + '\n')
repID = 1
while repID <= num_reps:
	print ("repID: " + str(repID))
	t_fixed = open("/mnt/storage/pjohri/" + s_folder + "/rep" + str(repID) + "/output_gen" + str(last_gen) + ".fixed", 'r')
	for line in t_fixed:
		line1 = line.strip('\n')
		line2 = line1.split()
		if len(line2) == 9 and "m" in line:
			if line2[2] == mutn_type:
				gen1 = int(line2[7])
				gen2 = int(line2[8])
				if gen2 > gen_burnin:
					sel = float(line2[4])
					time = gen2 - gen1 + 1
					s_posn = line2[3]
					result.write("rep" + str(repID) + '\t' + str(2.0*Ne*sel) + '\t' + str(2.0*Ne*B*sel) + '\t' + str(time) + '\t' + str(float(time)/float(Ne)) + '\t' + str(gen2) + '\t' + s_posn + '\n')
	t_fixed.close()
	repID += 1
result.close()
print ("done")








