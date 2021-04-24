#This is to separate .ms files by mutation types:
#m1 - f0
#m2 - f1
#m3 - f2
#m4 - f3
#m5 - beneficial
#How to run:
#python separate_mutn_types_fst.py -chromLen 10000 -mutnType m1 -input_folder /scratch/ekhowell/projects/deleterious_sweeps/Fst/neutral_arguello -output_folder /scratch/pjohri1/DelSweeps/FST/neutral_arguello
import sys
import math
import argparse
import os

#parsing user given constants
parser = argparse.ArgumentParser(description='Information about input/output folder and mutation type')
parser.add_argument('-chromLen', dest = 'chromLen', action='store', nargs = 1, type = int, help = 'chrom len')
parser.add_argument('-input_folder', dest = 'input_folder', action='store', nargs = 1, type = str, help = 'full path to folder with .ms files')
parser.add_argument('-output_folder', dest = 'output_folder', action='store', nargs = 1, type = str, help = 'full path to folder where you want to write the output')
#parser.add_argument('-output_prefix', dest = 'output_prefix', action='store', nargs = 1, type = str, help = 'full path to output file')
parser.add_argument('-mutnType', dest = 'mutnType', action='store', nargs = 1, type = str, help = 'type of mutation type- m1/m2/m5')
args = parser.parse_args()
chr_len =  args.chromLen[0]
infolder = args.input_folder[0]
outfolder = args.output_folder[0]
#prefix = args.output_prefix[0]
mutn_type = args.mutnType[0]

#read ms file:
def read_ms(f_ms):
    l_Pos = [] #list of positions of SNPs
    d_Alleles = {}#alleles at one site
    l_Int_Posn = []
    for line in f_ms:
        line1 = line.strip('\n')
        if "positions" in line1:
            line2 = line1.split()
            i = 0
            for x in line2:
                if "position" not in x:
                    l_Pos.append(float(x))
                    d_Alleles[str(i)] = ""
                    l_Int_Posn.append(i)
                    i = i + 1
        elif "//" not in line and "segsites" not in line:
            for j in l_Int_Posn:
                d_Alleles[str(j)] = d_Alleles[str(j)] + line1[j]
    f_ms.close()
    return (l_Pos, d_Alleles, l_Int_Posn)

#read full file:
def read_full_slim_output(f_full):
    d_Mutn_Type = {}
    mark = 0
    for line in f_full:
        line1 = line.strip('\n')
        if "Genomes:" in line1:
            mark = 0
        if mark == 1:
            line2 = line1.split()
            d_Mutn_Type[int(line2[3])] = line2[2]
        if "Mutations:" in line1:
            mark = 1
    f_full.close()
    return(d_Mutn_Type)
os.system("ls " + infolder + "/*.ms > " + outfolder + "/tmp.list")
f_list = open(outfolder + "/tmp.list", 'r')
numsim = 1
s_absent = 0
print(mutn_type)
for line in f_list:
    line1 = line.strip('\n')
    f_name = line1.split(".")[0]
    f_name1 = f_name.split("/").pop()
    print ("Reading file:" + line1)
    f_ms = open(f_name + ".ms", 'r')
    f_full = open(f_name + ".full", 'r')
    t_ms = read_ms(f_ms)
    l_posn = t_ms[0]
    d_alleles = t_ms[1]
    l_int_posn = t_ms[2]
    d_mutn_type = read_full_slim_output(f_full)
    #print (d_mutn_type)
    #finding relevant positions:
    l_posn_select = []
    l_int_posn_select = []
    for x in l_int_posn:
        #print (x)
        s_int_posn = round(float(l_posn[x])*int(chr_len))
        #print (s_int_posn)
        if d_mutn_type[s_int_posn] == mutn_type:
            l_posn_select.append(l_posn[x])
            l_int_posn_select.append(x)
    result = open(outfolder + "/" + f_name1 + "_" + mutn_type + ".ms", 'w+')
    result.write("\\" + '\n')
    result.write("segsites: " + str(len(l_int_posn_select)) + '\n')
    result.write("positions:")
    for y in l_posn_select:
        result.write(" " + str(y))
    result.write('\n')
    #write the genotypes:
    i = 0
    while i < len(d_alleles["0"]):
        for x in l_int_posn_select:
            result.write(d_alleles[str(x)][i])
        result.write('\n')
        i += 1
    result.close()

print ("done")


