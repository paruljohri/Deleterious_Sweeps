#Amova for the different evolutionary scenarios:
library("ape", lib.loc="~/R/win-library/3.5")
library("adegenet", lib.loc="~/R/win-library/3.5")
library("pegas", lib.loc="~/R/win-library/3.5")
setwd("../Work/DelSweeps/Amova/")
l_folders <- c("neutral_bneck_li_stephan", "strong_ben_dfe_fitrescaled", "weak_ben_dfe_fitrescaled", "del_dfe_fitrescaled", "del_dfe_fitnotrescaled", "neutral_arguello")
l_num_reps <- c(10, 20, 20, 100, 100, 100)
l_win_sizes <- c(500, 1000, 2000)
region_size <- 10000
#Run through different folders/ scenarios:
t_means <- c()
for (i in c(1:6)){
	s_folder <- l_folders[i]
	v_replicates <- c(1:l_num_reps[i])
	print(paste("folder name:", s_folder))
	#Run through different window sizes
	for (win_size in l_win_sizes){
		print(paste("window size:", win_size))
		v_windows <- c(1:(region_size/win_size))
		print (v_windows)
		
		#unzip folder:
		print("unzipping the folder...")
		unzip(paste("C:/Users/Parul Johri/Work/DelSweeps/Amova/", s_folder, "_FASTA_", win_size, "bp.zip", sep=""), exdir="C:/Users/Parul Johri/Work/DelSweeps/Amova")

		#Get phiST values for every window
		t_phi <- c()
		#Through all replicates
		for (repID in v_replicates){
			print(paste("replicate:", repID))
			#Through all sliding windows
			for (win in v_windows){
				m_fasta <- read.dna(paste("../Amova/", s_folder, "_FASTA_", win_size, "bp/rep", repID, "_win", win, ".fasta", sep=""), format = "fasta", skip = 0, nlines = 0, comment.char = "#", as.character = FALSE, as.matrix = TRUE)
				d_fasta <- dist.dna(m_fasta)
				p_fasta <- factor(c(rep("pop1", 100), rep("pop2", 100)))
				d_pop <- amova(d_fasta ~ p_fasta, nperm = 100)
				phiST <- getPhi(d_pop$varcomp$sigma2)[1,1]
				pval <- d_pop$varcomp$P.value[1]
				t_phi <- rbind(t_phi, c(repID, win, phiST, pval))
			}
		}
		#Write the file with phiST values
		colnames(t_phi) <- c("replicate", "window", "phiST", "pval")
		write.table(t_phi, file=paste(s_folder, "_", win_size, ".phiST", sep=""), append=F, sep='\t', row.names=F)
		
		#Look for outliers and write those out in another file
		t_phi_sorted <- t_phi[order(t_phi[,3], decreasing=T),]
		num_outliers = round((1/100)*length(t_phi_sorted))
		t_outliers <- t_phi_sorted[c(1:num_outliers),]
		write.table(t_outliers, file=paste(s_folder, "_outliers_1percent_", win_size,".phiST", sep=""), append=F, sep='\t', row.names=F)
		#Save the means:
		t_means <- rbind(t_means, c(s_folder, win_size, mean(t_phi[,3], na.rm=T), mean(t_outliers[,3]), mean(t_outliers[,3])/mean(t_phi[,3], na.rm=T)))
	}
}
colnames(t_means) <- c("folder", "window", "mean_phiST", "outlier_phiST", "outlier_by_mean")
write.table(t_means, file="mean_outlier_phiST.summary", append=F, sep='\t', row.names=F)



