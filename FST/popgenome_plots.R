#########################################################################
# Purpose:
# 1) Compare outlier and mean Fst values over certain
# window size and certain outlier threshold
# 2) Plot the allele frequencies for all specified replicates
# 3) Plot the outlier mutation type distribution for all
# specified replicates
# Mutation types: m0:s=0, m1:f0, m2:f1, m3:f2, m4:f3 m5:beneficial
# Requirements: uses _snp_list file generated for each replicate with
# create_snp_list.sh script. This matches mutation type with allele
# frequency
#########################################################################

# PACKAGES
library("PopGenome")
library("data.table")
library("ggplot2")
library("tidyverse")

# Set working directory to location of files
setwd("~/working_directory/")

#############################################################
# INPUTS (EDIT THESE!)

# Designate window size for Fst mean and outlier calculation
window <- 1

# Designate outlier significance
sig <- .01

# Number of replicates to run
rep <- 100
#############################################################

# ASSIGN VARIABLES

# Initialize global list variables
mut_tab <- data.frame()

# Initialize plot
par("mar" = c(5,5,5,5))
plot(x=0, ylim = c(0,1), xlim = c(0,1),
     ylab = "Population 1",
     xlab = "Population 2",
     cex.lab=2, cex.axis=1.75) #ylab:eur, xlab:afr

# FUNCTION

# Iterate i number of times
my_func <- function(i){
  
  # Simulation file names (index represents replicate number)
  mut_file <- paste0(c("simulation_prefix",as.character(i),"simulation_suffix"), collapse = "")
  
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
  
  # Getting Mutation Types #
  
  # Create a table of Fst values for each SNP window
  fst <- ms_pop_slide_fst@nucleotide.F_ST
  snp_list <- data.table(cbind(fst, ms_pop_slide_fst@region.names), key="`FST (Nucleotide)`")
  
  # Determine the number of outliers in the top XX%
  if (floor(length(snp_list$`FST (Nucleotide)`) * sig) == 0){
    outlier <- ceiling(length(snp_list$`FST (Nucleotide)`) * sig)
  } else {
    outlier <- floor(length(snp_list$`FST (Nucleotide)`) * sig)
  }

  # Grab out the SNP numbers for the outlier regions
  regions <- tail(snp_list, outlier)$V2
  regions <- as.list(strsplit(gsub('-',',', gsub(':','', gsub(' ','', regions))), ','))
  
  # Create a table with the SNP intervals for outlier regions
  entries <- seq(from=1, to=length(regions))
  from <- lapply(entries, function(x) {from <- (regions[[x]][1])})
  to <- lapply(entries, function(x) {from <- (regions[[x]][2])})
  outlier_df <- cbind(from,to)
  
  # Read in the list of mutations and their positions for this sim run
  muts <- read.delim(paste0(c(as.character(mut_file),"_snp_list"), collapse=""))
  
  # Get the mutation types for all SNPs contributing to the outlier regions
  my_func <- function(x)
  {
    rbind(lapply(seq(from=as.numeric(outlier_df[x,1]), to=as.numeric(outlier_df[x,2])), function(y){as.character.numeric_version(muts[y,]["X"])}))
  }
  
  # Create a table with the Fst, SNP interval, and list of contributing mutations
  mutations <- as.character(as.list(lapply(entries, function(z) {unlist(my_func(z))})))
  new <- cbind(outlier_df, as.data.frame(mutations))
  fst_value <- tail(snp_list, outlier)$`FST (Nucleotide)`
  full_table <- cbind(fst_value,new)
  
  # Create a list of all the mutations contributing to the outliers in this run
  mut_list <- unlist(as.list(lapply(entries, function(z) {unlist(my_func(z))})))
  dat <- data.frame(table(sapply(mut_list, function(x) x)))
  dat$Prop <- (dat$Freq)/length(mut_list)
  mut_tab <<- rbind(dat, mut_tab)
  
  # Plotting Allele Frequencies #
  
  # Get a list of the site positions for mutations in outlier regions
  my_func2 <- function (x) {
    lapply(seq(from=as.numeric(outlier_df[x,1]), to=as.numeric(outlier_df[x,2])), function(y){muts[y,]$Site})
  }
  
  # Create a list of the site positions for mutations in outlier regions
  site_list <- as.character(unlist(lapply(entries, function(x) {my_func2(x)})))

  # Create a new dataframe with the entries for outlier mutations
  subset <- subset(muts, muts$Site %in% site_list)
  
  # Create a new dataframe with the entries for non-outlier mutations
  diff <- setdiff(muts,subset)
  
  # Plot the data from this run onto the full graph
  points(diff$African, diff$European, pch = 20, col = "lightgray")
  
  # Add in colors to the outlier mutations based on mutation type
  points(subset$African, subset$European,
         col=ifelse(subset$X == "m1", "green", ifelse(subset$X == "m2", "red", ifelse(subset$X == "m3", "orange", ifelse(subset$X == "m4", "yellow", ifelse(subset$X == "m5", "blue", "lightgray"))))),
         pch=ifelse(subset$X == "m1", 1, ifelse(subset$X == "m2", 1, ifelse(subset$X == "m3", 1, ifelse(subset$X == "m4", 1, ifelse(subset$X == "m5", 1, 20))))))
  abline(a=0, b=1)
  legend("bottomright", title="Outlier Mutation Types",
         c("0","1","2", "3", "beneficial"), fill=c("green", "red", "orange", "yellow", "blue"),
         horiz = FALSE, cex = 1.5, xjust = 0.5, x.intersp = 0.15,
         bg = "transparent")
}

# APPLY FUNCTION

# Run this a specified number of times
list <- seq(from=1,to=rep)
lapply(list, function(x){
  my_func(x)
})

# MUTATION TYPE DISTRIBUTION

# Calculate mean for each mutation type present
mean_fun <- function(x){
  sum(x)/rep
}

# Calculate standard deviation for each mutation type present
sd_fun <- function(x){
  mean <- sum(x)/rep
  sqrt(sum(((x-mean)^2)/(rep-1)))
}

# Merge tables for standard deviation and mean
sd_tab <- aggregate(Prop ~ Var1, data=mut_tab, FUN=sd_fun)
colnames(sd_tab) <- c("mut", "sd")
mean_tab <- aggregate(Prop ~ Var1, data=mut_tab, FUN=mean_fun)
colnames(mean_tab) <- c("mut", "mean")
table <- merge(sd_tab, mean_tab, by="mut")

# Plot the frequency of each mutation type
ggplot(table, aes(x = mut, y = mean)) +
  geom_bar(stat="identity", fill = c("green", "red", "orange", "yellow", "blue")) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  ylab("Proportion") + xlab("Mutation Type") +
  theme_linedraw(base_size = 40) +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 38),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0))
