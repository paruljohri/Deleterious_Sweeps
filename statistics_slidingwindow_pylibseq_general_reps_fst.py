#Basic stats, sliding window, SLIm ms output
#adding divergence to this
#python statistics_slidingwindow_pylibseq_general_reps_fst.py -winSize 500 -stepSize 500 -regionLen 10000 -input_folder /scratch/ekhowell/projects/deleterious_sweeps/Fst/neutral_arguello -output_folder /scratch/pjohri1/DelSweeps/FST -output_prefix neutral_arguello
#ben_dfe_fitrescaled/  del_dfe_fitnotrescaled/  del_dfe_fitrescaled/  neutral_arguello/  neutral_bneck_li_stephan/

from __future__ import print_function
import libsequence
import sys
import pandas
import math
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about number of sliding windows and step size')
parser.add_argument('-winSize', dest = 'winSize', action='store', nargs = 1, default = 100, type = int, help = 'size of each sliding window in bp')#500 bp for small, 5000bp for big
parser.add_argument('-stepSize', dest = 'stepSize', action='store', nargs = 1, default = 100, type = int, help = 'size of step size in bp')#250 bp for small, 5000 bp for big
parser.add_argument('-regionLen', dest = 'regionLen', action='store', nargs = 1, type = int, help = 'length in bp of region simulated')#Length of coding region simulated
parser.add_argument('-numIndvPop1', dest = 'numIndvPop1', action='store', nargs = 1, type = int, help = 'number of haploid genomes from pop1')
parser.add_argument('-numIndvPop2', dest = 'numIndvPop2', action='store', nargs = 1, type = int, help = 'number of haploid genomes from pop2')
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
args = parser.parse_args()
chr_len =  args.regionLen[0]
win_size_bp = args.winSize[0]
win_size = args.winSize[0]/float(chr_len)
step_size = args.stepSize[0]/float(chr_len)
num_indv1 = args.numIndvPop1[0]
num_indv2 = args.numIndvPop2[0]
infolder = args.input_folder[0]
outfolder = args.output_folder[0]
prefix = args.output_prefix[0]

def read_fixed_mutations(f_fixed):
    d_subs = {}
    for line in f_fixed:
        line1 = line.strip('\n')
        line2 = line1.split()
        if line1[0]!="#" and line2[0]!="Mutations:":
            posn = float(line2[3])/float(chr_len)
            #if "grow10" in input_folder  or "red10" in input_folder:
            num_gen = line2.pop()
            #if int(num_gen) >= 100000:
            d_subs[posn] = d_subs.get(posn, 0) + 1
            #else:
            #    d_subs[posn] = d_subs.get(posn, 0) + 1
    return d_subs #return a dictionary with base positions as keys and number of fixed substitutions as values

def avg_divergence_win(d_subs, start, end):
    s_sum = 0
    for posn in d_subs.keys():
        if float(posn) <= end and float(posn) > start:
            s_sum = s_sum + 1
            
    return s_sum

def get_S(f_ms):
    for line in f_ms:
        line1 = line.strip('\n')
        if "segsites" in line1:
            S = line1.split()[1]
    f_ms.seek(0)
    return S
#result files:

result =  open(outfolder + "/" + prefix + "_" + str(args.winSize[0]) +  ".stats", 'w+')
result.write("repID" + '\t' + "posn" + '\t' + "S" + '\t' + "FST_HSM" + '\t' + "FST_Slatkin" + '\t' + "FST_HBK" + '\t' + "SNPs_shared" + '\t' + "SNPs_private1" + '\t' + "SNPs_private2" + '\t' + "Fixed_diff" + '\t' + "pop1_thetapi" + '\t' + "pop2_thetapi" + '\t' + "pop1_thetaw" + '\t' + "pop2_thetaw" + '\t' + "pop1_thetah" + '\t' + "pop2_thetah" + '\t' + "pop1_tajimasd" + '\t' + "pop2_tajimasd" + '\t' + "pop1_numsing" + '\t' + "pop2_numsing" + '\t' + "pop1_hapdiv" + '\t' + "pop2_hapdiv" + '\n')

#go through all simulation replicates and read data into pylibseq format
#addin the option of ignoring some files if they don't exist
os.system("ls " + infolder + "/*.ms > " + outfolder + "/" + prefix + ".list")
f_list = open(outfolder + "/" + prefix + ".list", 'r')
numsim = 1
s_absent = 0
for line in f_list:
    line1 = line.strip('\n')
    f_name = line1.split(".")[0]
    print ("Reading file:" + line1)
    #try:
    if numsim > 0:
        f_ms = open(f_name + ".ms", 'r')
        #f_subs = open(f_name + ".fixed", 'r')
        #d_subs = read_fixed_mutations(f_subs)
        #S = get_S(f_ms)
        l_Pos = [] #list of positions of SNPs
        l_Genos = [] #list of alleles
        d_tmp, d_tmp_pop1, d_tmp_pop2 = {}, {}, {}
        genome_count = 0 #To count the number of individuals in the ms file
        for line in f_ms:
            line1 = line.strip('\n')
            if "positions" in line1:
                line2 = line1.split()
                i = 0
                for x in line2:
                    if "position" not in x:
                        l_Pos.append(float(x))
                        d_tmp[str(i)] = ""
                        d_tmp_pop1[str(i)] = ""
                        d_tmp_pop2[str(i)] = ""
                        i = i + 1
            elif "//" not in line and "segsites" not in line and len(line1) > 0 and "\\" not in line:
                #print (line)
                #print (d_tmp)
                genome_count += 1
                i = 0
                while i < len(line1):
                    d_tmp[str(i)] = d_tmp[str(i)] + line1[i]
                    if genome_count <= num_indv1:
                        d_tmp_pop1[str(i)] = d_tmp_pop1[str(i)] + line1[i]
                    elif genome_count > num_indv1:
                        d_tmp_pop2[str(i)] = d_tmp_pop2[str(i)] + line1[i]
                    i = i + 1
        #print (d_tmp)
        l_data = []
        l_data_pop1, l_data_pop2 = [], []
        i = 0
        while i < len(l_Pos):
            l_Genos.append(d_tmp[str(i)])
            t_tmp = (l_Pos[i], d_tmp[str(i)])
            t_tmp_pop1 = (l_Pos[i], d_tmp_pop1[str(i)])
            t_tmp_pop2 = (l_Pos[i], d_tmp_pop2[str(i)])
            l_data.append(t_tmp)
            l_data_pop1.append(t_tmp_pop1)
            l_data_pop2.append(t_tmp_pop2)
            i = i + 1
        #print (l_Pos)
        #print (l_Genos)


        #assign object
        sd = libsequence.SimData(l_data)
        sd_pop1 = libsequence.SimData(l_data_pop1)
        sd_pop2 = libsequence.SimData(l_data_pop2)
        print (sd_pop1.size())
        print (sd_pop2.size())
        #sd.assign(l_Pos[10:100],l_Genos[10:100])

        #define sliding windows:
        w = libsequence.Windows(sd,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
        w_pop1 = libsequence.Windows(sd_pop1,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
        w_pop2 = libsequence.Windows(sd_pop2,window_size=win_size,step_len=step_size,starting_pos=0.0,ending_pos=1.0)
        #chromosome length = 30kb, window size = 5 kb
        num_win = len(w)

        #calculate summary statistic in sliding window:
        print ("calculating stats in windows")
        win_name = 1
        for i in range(len(w)):
            wi = w[i]
            w_pop1i = w_pop1[i]
            w_pop2i = w_pop2[i]
            #print (wi)
            psw_pop1i = libsequence.PolySIM(w_pop1i)
            psw_pop2i = libsequence.PolySIM(w_pop2i)
            try:
                fwi = libsequence.Fst(wi,[num_indv1,num_indv2])
                result.write("rep" + str(numsim) + '\t' + str(win_name) + '\t' + str(len(wi.pos())) + '\t' +  str(fwi.hsm()) + '\t' + str(fwi.slatkin()) + '\t' + str(fwi.hbk()) + '\t' + str(len(fwi.shared(0,1))) + '\t' + str(len(fwi.priv(0,1)[0])) + '\t' + str(len(fwi.priv(0,1)[1])) + '\t' + str(len(fwi.fixed(0,1))) + '\t')
            except:
                result.write("rep" + str(numsim) + '\t' + str(win_name) + '\t' + str(len(wi.pos())) + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t' + "NA" + '\t')
            result.write(str(psw_pop1i.thetapi()/float(win_size_bp)) + '\t' + str(psw_pop2i.thetapi()/float(win_size_bp)) + '\t' + str(psw_pop1i.thetaw()/float(win_size_bp)) + '\t' + str(psw_pop2i.thetaw()/float(win_size_bp)) + '\t' + str(psw_pop1i.thetah()/float(win_size_bp)) + '\t' + str(psw_pop2i.thetah()/float(win_size_bp)) + '\t' + str(psw_pop1i.tajimasd()/float(win_size)) + '\t' + str(psw_pop2i.tajimasd()) + '\t' + str(psw_pop1i.numexternalmutations()/float(win_size_bp)) + '\t' + str(psw_pop2i.numexternalmutations()/float(win_size_bp)) + '\t' + str(psw_pop1i.hapdiv()) + '\t' + str(psw_pop2i.hapdiv()) + '\n')
            #divergence:
            #s_start = (i)*(1.0/float(num_win))
            #s_end = s_start + win_size
            #result.write(str(avg_divergence_win(d_subs, s_start, s_end)) + '\n')
            win_name = win_name + 1
        
    #except:
    else:
        s_absent = s_absent + 1
        print ("This file does not exist or cannot be read or is empty")
    
    numsim = numsim + 1

print ("Number of files not read:" + '\t' + str(s_absent))
print ("Finished")






