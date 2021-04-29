#Plotting:
setwd("../Work/DelSweeps/single_sweeps_dfe/")
library("ggplot2", lib.loc="~/R/win-library/3.5")

#For Drosophila:

#Read Tajima-simulation predictions:
t_tajima <- read.table("../Theory/pi_reduction_sweeps_Tajimasims_2.txt", h=T)

#For average rec:
t_pi <- read.table("Droso_N10k_mean_rec/pi_postfixation_10.txt", h=T)
rec <- (1e-8)*195
pi_neu <- 0.01899446

#For 0.5 x average rec
t_pi <- read.table("Droso_N10k_half_rec/pi_postfixation_10.txt", h=T)
rec <- (1e-8)*195*0.5
pi_neu <- 0.01851208

#For 0.1 x average rec
t_pi <- read.table("Droso_N10k_tenth_rec/pi_postfixation_10.txt", h=T)
rec <- (1e-8)*195*0.1
pi_neu <- 0.01563094

#For all recombination rates:
u <- (3e-9)*195
t_posn <- seq(-495, 500,10)
N <- 1e4
t_rho <- 2*N*rec*t_posn
B <- pi_neu/(4*N*u)
yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))


#gamma=1-2 with Tajima sims:
t_pi_1 <- subset(t_pi, gamma > -2.0 & gamma < -1.0)
t_mean <- colMeans(t_pi_1[, c(4:dim(t_pi)[2])]/pi_neu, na.rm=TRUE)
t_se <- apply(t_pi_1[, c(4:dim(t_pi)[2])]/pi_neu,2,sd)/sqrt(dim(t_pi_1)[1])

t_data <- data.frame(t_posn, t_rho, t_mean, t_se)
colnames(t_data) <- c("bases", "rho", "thetapi_m", "thetapi_se")

t_data$min_se <- t_data$thetapi_m - t_data$thetapi_se
t_data$max_se <- t_data$thetapi_m + t_data$thetapi_se

gg <- ggplot(t_data) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "khaki3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="red4") + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_0_semi, ymin=1.0-thetapi_m_gamma_0_semi-thetapi_se_gamma_5_semi, ymax=1.0-thetapi_m_gamma_0_semi+thetapi_se_gamma_0_semi), color="red3", shape=4, size=1) + coord_cartesian(xlim=c(min(t_data$rho, na.rm=T),max(t_data$rho, na.rm=T)), ylim=c(0.5, 1.1))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)

#gamma=-5 with Tajima's sims:
t_pi_5 <- subset(t_pi, gamma > -5.5 & gamma < -4.5)
t_mean <- colMeans(t_pi_5[, c(4:dim(t_pi)[2])]/pi_neu, na.rm=TRUE)
t_se <- apply(t_pi_5[, c(4:dim(t_pi)[2])]/pi_neu,2,sd)/sqrt(dim(t_pi_5)[1])

t_data <- data.frame(t_posn, t_rho, t_mean, t_se)
colnames(t_data) <- c("bases", "rho", "thetapi_m", "thetapi_se")

t_data$min_se <- t_data$thetapi_m - t_data$thetapi_se
t_data$max_se <- t_data$thetapi_m + t_data$thetapi_se

gg <- ggplot(t_data) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "khaki3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="red4") + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_5_semi, ymin=1.0-thetapi_m_gamma_5_semi-thetapi_se_gamma_5_semi, ymax=1.0-thetapi_m_gamma_5_semi+thetapi_se_gamma_5_semi), color="red3", shape=4, size=1) + coord_cartesian(xlim=c(min(t_data$rho, na.rm=T),max(t_data$rho, na.rm=T)), ylim=c(0.5, 1.2))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)

#For Humans:
rm(list = ls())

#For average rec:
t_pi <- read.table("Human_N10k_mean_rec/pi_postfixation_del_200.txt", h=T)
rec <- 1.0e-8
pi_neu <- 0.0004174326

#For 0.5 x average rec
t_pi <- read.table("Human_N10k_half_rec/pi_postfixation_del_200.txt", h=T)
rec <- (1.0e-8)*0.5
pi_neu <- 0.0004082227

#For 0.1 x average rec
t_pi <- read.table("Human_N10k_tenth_rec/pi_postfixation_del_200.txt", h=T)
rec <- (1.0e-8)*0.1
pi_neu <- 0.0003982966

#For all recombination rates:
u <- 1.2e-8
t_posn <- seq(-39900, 40000,200)
N <- 1e4
B <- pi_neu/(4*N*u)
t_rho <- 2*N*rec*t_posn
yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))

#Plotting gamma=1-2
t_pi_sub <- subset(t_pi, gamma >= -2 & gamma < -1.0)
t_pi_sub1 <- t_pi_sub[, c(4:dim(t_pi)[2])]
t_mean <- colMeans(t_pi_sub1/pi_neu, na.rm=TRUE)
t_sd <- apply(t_pi_sub1/pi_neu,2,sd, na.rm=T)
t_len <- c()
i <- 1
while (i<=dim(t_pi_sub1)[2]){
	t_len <- c(t_len, length(na.omit(t_pi_sub1[,i])))
	i <- i + 1}
t_se <- t_sd/sqrt(t_len)

t_data <- data.frame(t_posn, t_rho, t_mean, t_se)
colnames(t_data) <- c("bases", "rho", "thetapi_m", "thetapi_se")

t_data$min_se <- t_data$thetapi_m - t_data$thetapi_se
t_data$max_se <- t_data$thetapi_m + t_data$thetapi_se

gg <- ggplot(t_data) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "khaki3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="red4") + coord_cartesian(ylim=c(0.2, 1.8))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)

#Plotting gamma<4:
t_pi_sub <- subset(t_pi, gamma < -4.0)
t_pi_sub1 <- t_pi_sub[, c(4:dim(t_pi)[2])]
t_mean <- colMeans(t_pi_sub1/pi_neu, na.rm=TRUE)
t_sd <- apply(t_pi_sub1/pi_neu,2,sd, na.rm=T)
t_len <- c()
i <- 1
while (i<=dim(t_pi_sub1)[2]){
	t_len <- c(t_len, length(na.omit(t_pi_sub1[,i])))
	i <- i + 1}
t_se <- t_sd/sqrt(t_len)

t_data <- data.frame(t_posn, t_rho, t_mean, t_se)
colnames(t_data) <- c("bases", "rho", "thetapi_m", "thetapi_se")

t_data$min_se <- t_data$thetapi_m - t_data$thetapi_se
t_data$max_se <- t_data$thetapi_m + t_data$thetapi_se

gg <- ggplot(t_data) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "khaki3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="red4") + coord_cartesian(ylim=c(0.0, 3.0))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)

