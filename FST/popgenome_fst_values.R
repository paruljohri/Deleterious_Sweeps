##########################################################################
# Purpose: compare outlier and mean Fst values over certain
# window size and certain outlier threshold
# Note) Mutation types: m0:s=0, m1:f0, m2:f1, m3:f2, m4:f3, m5:beneficial
##########################################################################

# PACKAGES
library("PopGenome")
library("data.table")
library("ggplot2")
library("tidyverse")

# Set working directory to location of files
setwd("~/working_directory/")

# ASSIGN VARIABLES

# Initialize global list variables
avg_outlier <- list()
mean_slide_fst <- list()

# FUNCTION

# Iterate i number of times
my_func <- function(i,j,k,l){
  
  ##############################################################################################
  # INPUTS (Edit these!)
  
  # Simulation file names (index represents replicate number)
  mut_file <- paste0(c("simulation_prefix",as.character(i),"simulation_suffix"), collapse = "")
  
  # Designate window size for Fst mean and outlier calculation
  window <- 10
  
  # Designate outlier significance
  sig <- .01
  ##############################################################################################

  # Calculating Fst #
  
  # Read in the MS file for the polymorphic sites only
  ms <- readMS(paste0(c(as.character(mut_file),".ms"), collapse = ""), big.data = TRUE)
  
  # Create lists of the individuals belonging to the African and European populations
  afr <- as.character(seq(from = 1, to = 100))
  eur <- as.character(seq(from = 101, to = 200))
  
  # Assign populations in the MS file
  ms_pop <- set.populations(ms, list(afr,eur))
  
  # Cut the sites into sliding windows based on the number of SNPs in each window
  ms_pop_slide <- sliding.window.transform(ms_pop, width=window, jump=window, type=1, whole.data=TRUE)
  
  # Calculate the Fst value for each sliding window
  ms_pop_slide_fst <- F_ST.stats(ms_pop_slide)
  
  # Calculate mean sliding window Fst
  mean_slide_fst[i] <<- (mean(get.F_ST(ms_pop_slide_fst)[,2]))
  
  # Getting Mutation Types #
  
  # Create a table of Fst values for each SNP window
  fst <- ms_pop_slide_fst@nucleotide.F_ST
  snp_list <- data.table(cbind(fst, ms_pop_slide_fst@region.names), key="`FST (Nucleotide)`")
  
  # Determine the number of outliers in the top X%
  if (floor(length(snp_list$`FST (Nucleotide)`) * sig) == 0){
    outlier <- ceiling(length(snp_list$`FST (Nucleotide)`) * sig)
  } else {
    outlier <- floor(length(snp_list$`FST (Nucleotide)`) * sig)
  }
  fst_value <- tail(snp_list, outlier)$`FST (Nucleotide)`
  
  # Calculating Average Outlier Fst #
  
  # Determine the average Fst value of the outlier region
  avg_outlier[i] <<- (mean(as.numeric(unlist(fst_value))))
}

# APPLY FUNCTION

# Specify how many replicates to run this operation for 
num <- seq(from=1,to=100)
lapply(num, function(x){
  my_func(i)
})

# Print output to screen
paste0(c("Outlier: ",mean(unlist(avg_outlier)),", Mean: ",(mean(unlist(mean_slide_fst)))), collapse = "")
