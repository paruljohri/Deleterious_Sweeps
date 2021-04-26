# Purpose: combine the file containing SNP frequencies in the African and European population with the file containing mutation types
# This allows us to track the allele frequencies (and therefore outlier contribution) for segregating sites of different mutation types
# This script uses the get_allele_freq.py to match mutation type (m1, m2, m3, m4, m5) with allele frequency

# Run with: create_snp_list.sh
# Note: must be run in folder containing simulation output

# Filenames are specified as "ms_file_prefix_replicate_ms_file_suffix.ms" and "ms_file_prefix_replicate_ms_file_suffix.full"
# Prefix of SLiM output
FILE_PREFIX="ms_file_prefix" 
# Suffix of SLiM output
FILE_SUFFIX="ms_file_suffix"

# Index represents the replicate number (assumming here that 100 replicates have been run)
for i in {1..100}
do
	# Obtain allele frequencies for each SNP in the SLiM ".full" output
	python get_allele_freq.py --input "${FILE_PREFIX}${i}${FILE_SUFFIX}.ms" --output "${FILE_PREFIX}${i}${FILE_SUFFIX}.allele_freq"
	# Merge with original ".full" output to connect mutation type with allele frequencies
	# Outfile has extension _snp_list
	cat "${FILE_PREFIX}${i}${FILE_SUFFIX}.full" | awk 'FNR > 2 {if($0~/m/) print $3"\t"$4"\t"$9}' | sort -k2 -n | paste -d"\t" /dev/stdin "${FILE_PREFIX}${i}${FILE_SUFFIX}.allele_freq" | awk '{if ($3!="200") print $0}' > "${FILE_PREFIX}${i}${FILE_SUFFIX}_snp_list"
done
