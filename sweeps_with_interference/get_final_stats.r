
#This is to get final stats:
args = commandArgs(trailingOnly=TRUE)
s_folder <- args[1]

N <- 10000
num_gen <- (100*N) - (10*N)

if (grepl("Droso", s_folder, fixed = TRUE)==T){
	#For Drosophila:
	num_reps <- 100
	u <- 5.58e-7
	num_exons <- 5
	num_introns <- 4
	num_genes <- 2
	exon_len <- 300
	intron_len <- 100
	inter_len <- 4000
	f0 <- 0.25
	f1 <- 0.49
	}
if (grepl("Human", s_folder, fixed = TRUE)==T){
	#For humans:
	num_reps <- 300
	u <- 1.2e-8
	num_exons <- 5
        num_introns <- 4
        num_genes <- 2
        exon_len <- 300
        intron_len <- 2000
        inter_len <- 15000
        f0 <- 0.51
        f1 <- 0.14
	}

#define some more constants:
tot_exon_len <- exon_len*num_exons*num_genes
tot_nonexon_len <- (inter_len*(num_genes+1)) + (intron_len*num_introns*num_genes)

t_time <- c()
t_pfix <- c()

#For tehoretical calcualtions:
Pfix <- function(x) {(x)/(exp(2*N*x) - 1)}

#For neutral coding:
t <- read.table(paste("/home/pjohri/DelSweeps/", s_folder, "/time_to_fixation_m1.txt", sep=""),h=T)
t_mean <- aggregate((t$time_Ngen/2)~t$repID, FUN=mean)
t_time <- rbind(t_time, c("0-1syn", mean(t_mean[,2]), quantile(t_mean[,2], 0.25)))

gamma1 <- 0
gamma2 <- 1
s1 <- gamma1/(2*N)
s2 <- gamma2/(2*N)
Integral <- integrate(Pfix, lower = s1, upper = s2)
MeanProb <- Integral$value*(1/(s2 - s1))

t_num_sub <- aggregate(t$time_Ngen~t$repID, FUN=length)
t_prob_fix <- t_num_sub[,2]/(2*N*u*tot_exon_len*num_gen*f0)
s_num_subs <- length(t$time_Ngen)
s_prob_fix_pooled <- s_num_subs/(2*N*u*tot_exon_len*num_gen*f0*num_reps)
t_pfix <- rbind(t_pfix, c("0-1syn", MeanProb, s_num_subs, s_prob_fix_pooled, mean(t_prob_fix), sd(t_prob_fix)))


#For neutral intergenic:
t <- read.table(paste("/home/pjohri/DelSweeps/", s_folder,"/time_to_fixation_m5.txt", sep=""),h=T)
t_mean <- aggregate((t$time_Ngen/2)~t$repID, FUN=mean)
t_time <- rbind(t_time, c("0-1int", mean(t_mean[,2]), quantile(t_mean[,2], 0.25)))

gamma1 <- 0
gamma2 <- 1
s1 <- gamma1/(2*N)
s2 <- gamma2/(2*N)
Integral <- integrate(Pfix, lower = s1, upper = s2)
MeanProb <- Integral$value*(1/(s2 - s1))

t_num_sub <- aggregate(t$time_Ngen~t$repID, FUN=length)
t_prob_fix <- t_num_sub[,2]/(2*N*u*tot_nonexon_len*num_gen)
s_num_subs <- length(t$time_Ngen)
s_prob_fix_pooled <- s_num_subs/(2*N*u*tot_nonexon_len*num_gen*num_reps)
t_pfix <- rbind(t_pfix, c("0-1", MeanProb, s_num_subs, s_prob_fix_pooled, mean(t_prob_fix), sd(t_prob_fix)))


#For mildly deleterious coding:
t <- read.table(paste("/home/pjohri/DelSweeps/", s_folder, "/time_to_fixation_m2.txt", sep=""),h=T)
gamma_start <- 1
gamma_end <- 10
s_start <- gamma_start/(2*N)
s_end <- gamma_end/(2*N)

gamma1 <- 1
gamma2 <- 2
while (gamma1 < 10){
	s1 <- gamma1/(2*N)
	s2 <- gamma2/(2*N)

	Integral <- integrate(Pfix, lower = s1, upper = s2)
        MeanProb <- Integral$value*(1/(s2 - s1))
	
	s_num_subs <- length(t$time_Ngen[which(t$gamma>-1.0*gamma2 & t$gamma <= -1.0*gamma1)])
	if (s_num_subs > 0){
		t_mean <- aggregate((t$time_Ngen[which(t$gamma>-1.0*gamma2 & t$gamma <= -1.0*gamma1)]/2)~t$repID[which(t$gamma>-1.0*gamma2 & t$gamma <= -1.0*gamma1)], FUN=mean)
		t_time <- rbind(t_time, c(paste(gamma1, "-", gamma2, sep=""), mean(t_mean[,2]), quantile(t_mean[,2], 0.25)))
		
		t_num_sub <- aggregate(t$time_Ngen[which(t$gamma>-1.0*gamma2 & t$gamma <= -1.0*gamma1)]~t$repID[which(t$gamma>-1.0*gamma2 & t$gamma <= -1.0*gamma1)], FUN=length)
		t_prob_fix <- t_num_sub[,2]/(2*N*u*tot_exon_len*num_gen*f1*((s2-s1)/(s_end-s_start)))
		s_prob_fix_pooled <- s_num_subs/(num_reps*2*N*u*tot_exon_len*num_gen*f1*((s2-s1)/(s_end-s_start)))
		t_pfix <- rbind(t_pfix, c(paste(gamma1, "-", gamma2, sep=""), MeanProb, s_num_subs, s_prob_fix_pooled, mean(t_prob_fix), sd(t_prob_fix)))
		}
	else{
		t_time <- rbind(t_time, c(paste(gamma1, "-", gamma2, sep=""), "NA", "NA"))
		t_pfix <- rbind(t_pfix, c(paste(gamma1, "-", gamma2, sep=""), MeanProb, s_num_subs, "NA", "NA", "NA"))
		}
	gamma1 <- gamma1 + 1
	gamma2 <- gamma2 + 1
	}

colnames(t_time) <- c("gamma", "mean", "first_quartile")
colnames(t_pfix) <- c("gamma", "pfix_kimura", "num_subs", "pfix_total", "pfix_mean", "pfix_sd")
write.table(t_time, file=paste("/home/pjohri/DelSweeps/", s_folder, "/time_to_fixation.stats", sep=""), append=F, sep='\t', row.names=F)
write.table(t_pfix, file=paste("/home/pjohri/DelSweeps/", s_folder, "/prob_of_fixation.stats", sep=""), append=F, sep='\t', row.names=F)


