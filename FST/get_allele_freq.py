# Purpose: for each segregating site contained within the .full SLiM output, print the allele frequency in the African and European population
# This script is used within the create_snp_list.sh script to match mutation type (m1, m2, m3, m4, m5) with allele frequency

# Import modules
import argparse

# Parse command line input
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='path of MS file containing polymorphic sites')
parser.add_argument('--output', required=True, help='path of output file')
args = parser.parse_args()

# Open input and output files
infile = open(str(args.input), "r")
outfile = open(str(args.output), "w")

# Initialize lists to hold names of segregating sites
pop1 = []
pop2 = []

# Initialize a dictionary to hold the name of the segregating site along with its frequency
pop1_dict = {}
pop2_dict = {}

# Create the list of sergegating site names and the dictionary with names as keys and frequencies as values
# "Duplicate" site names can arise when different mutation types hit the same spot
# Assign unique tag to duplicated segregating sites. Note: duplicates cannot exceed number of mutation types

# Start at line with positions
line = infile.readlines()[3]
record = line.split()
# Set index to 1
i = 1
# For each position in the line
while (i < len(record)):
	# Check to see if there are 4 duplicates
	if (((i + 3) < len(record)) and (record[i] == record[i + 3])):
		site1 = str(round((float(record[i]) * 10000))) + "_a"
		site2 = str(round((float(record[i + 1]) * 10000))) + "_b"
		site3 = str(round((float(record[i + 2]) * 10000))) + "_c"
		site4 = str(round((float(record[i + 3]) * 10000))) + "_d"
		site_list = site1, site2, site3, site4
		for site in site_list:
			pop1.append(site)
			pop2.append(site)
			pop1_dict[site] = 0
			pop2_dict[site] = 0
		i += 4
		continue
	# Check to see if there are 3 duplicates
	elif (((i + 2) < len(record)) and (record[i] == record[i + 2])):
		site1 = str(round((float(record[i]) * 10000))) + "_a"
		site2 = str(round((float(record[i + 1]) * 10000))) + "_b"
		site3 = str(round((float(record[i + 2]) * 10000))) + "_c"
		site_list = site1, site2, site3
		for site in site_list:
			pop1.append(site)
			pop2.append(site)
			pop1_dict[site] = 0
			pop2_dict[site] = 0
		i += 3
		continue
	# Check to see if there are 2 duplicates
	elif (((i + 1) < len(record)) and (record[i] == record[i + 1])):
		site1 = str(round((float(record[i]) * 10000))) + "_a"
		site2 = str(round((float(record[i + 1]) * 10000))) + "_b"
		site_list = site1, site2
		for site in site_list:
			pop1.append(site)
			pop2.append(site)
			pop1_dict[site] = 0
			pop2_dict[site] = 0
		i += 2
		continue
	# In case of no duplicates
	else:
		site = round((float(record[i]) * 10000))
		pop1.append(site)
		pop2.append(site)
		pop1_dict[str(site)] = 0
		pop2_dict[str(site)] = 0
		i += 1

# Add the occurances of each segregating site to the value of its key in the dictionary

# Set individual count to 0 to keep track of pop1 and pop2 individuals
indv_num = 0
# Set the file to line 0
infile.seek(0)
# For each line that contains genotypes
for line in infile:
	if ("//" not in line and "segsites" not in line and "positions" not in line and ("0" in line or "1" in line)):
		# Break it into elements
		elements = list(line.strip())
		# Set element counter to 0
		j = 0
		# If we are still on individuals in the first population
		if (indv_num < 100):	
			while (j < len(elements)):
				# Populate dictionary keys with their corresponding genotypes
				for x in pop1_dict:
					if (str(pop1[j]) == str(x)):
						pop1_dict[x] = pop1_dict[x] + int(elements[j])
						j += 1
		# If we are on individuals in the second population
		if (indv_num >=  100):	
			while (j < len(elements)):
				# Populate dictionary keys with their corresponding genotypes				
				for x in pop2_dict:
					if (str(pop2[j]) == str(x)):
						pop2_dict[x] = pop2_dict[x] + int(elements[j])
						j += 1
		# Increase individual counter by 1
		indv_num += 1

# Output the site name (as identified in .full output from SLiM), and the frequencies in each population
outfile.write("Site" + "\t" + "African" + "\t" + "European" + "\n")
for y in pop1_dict:
	outfile.write(str(y) + "\t" + str(pop1_dict[y]/100) + "\t" + str(pop2_dict[y]/100) + "\n")

# Output file has extension .allele_freq
outfile.close()

