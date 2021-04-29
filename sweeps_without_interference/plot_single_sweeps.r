#To plot single sweep results
#Use geom_line() if you can
setwd("../Work/DelSweeps/single_sweeps/")
library("ggplot2", lib.loc="~/R/win-library/3.5")

#gamma=+5
t_pi <- read.table("gamma_pos_five_10.final", h=T)
t_th <- read.table("../Theory/pi_reduction_sweeps_brian.txt", h=T)
t_tajima <- read.table("../Theory/pi_reduction_sweeps_Tajimasims_2.txt", h=T)

t_pi$min_se <- t_pi$thetapi_m - t_pi$thetapi_se
t_pi$max_se <- t_pi$thetapi_m + t_pi$thetapi_se
yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))
gg <- ggplot(t_pi) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "darkseagreen3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="blue") + geom_point(data=t_th, aes(x=rho, y=1.0-thetapi_m), color="dark blue", shape=16, size=3) + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_5_semi, ymin=1.0-thetapi_m_gamma_5_semi-thetapi_se_gamma_5_semi, ymax=1.0-thetapi_m_gamma_5_semi+thetapi_se_gamma_5_semi), color="dark blue", shape=4, size=1) + coord_cartesian(ylim=c(0.3, max(t_pi$thetapi_m)), xlim=c(-13,13))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)
>>save as 8.27 x 5
#+ geom_point(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_5_semi), color="dark blue")
#gamma=-5
t_pi <- read.table("gamma_neg_five_10.final", h=T)
t_th <- read.table("../Theory/pi_reduction_sweeps_brian.txt", h=T)
t_tajima <- read.table("../Theory/pi_reduction_sweeps_Tajimasims_2.txt", h=T)

t_pi$min_se <- t_pi$thetapi_m - t_pi$thetapi_se
t_pi$max_se <- t_pi$thetapi_m + t_pi$thetapi_se
yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))
gg <- ggplot(t_pi) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "khaki3", color = "transparent") +geom_line(aes(x=rho, y=thetapi_m), color="red4") + geom_point(data=t_th, aes(x=rho, y=1.0-thetapi_m), color="red3", shape=16, size=3) + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_5_semi, ymin=1.0-thetapi_m_gamma_5_semi-thetapi_se_gamma_5_semi, ymax=1.0-thetapi_m_gamma_5_semi+thetapi_se_gamma_5_semi), color="red3", shape=4, size=1) + coord_cartesian(ylim=c(0.3, max(t_pi$thetapi_m)), xlim=c(-13,13))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)
>>save as 8.27 x 5
#geom_point(data=t_tajima, aes(x=rho, y=1.0-thetapi_m), color="red3")
#gamma=0
t_pi <- read.table("gamma_zero_10.final", h=T)
t_tajima <- read.table("../Theory/pi_reduction_sweeps_Tajimasims_2.txt", h=T)
t_pi$min_se <- t_pi$thetapi_m - t_pi$thetapi_se
t_pi$max_se <- t_pi$thetapi_m + t_pi$thetapi_se

yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))
gg <- ggplot(t_pi) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "grey60", color = "transparent") +geom_line(aes(x=rho, y=thetapi_m), color="black") + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_0_semi, ymin=1.0-thetapi_m_gamma_0_semi-thetapi_se_gamma_0_semi, ymax=1.0-thetapi_m_gamma_0_semi+thetapi_se_gamma_0_semi), color="black", shape=4, size=1) + coord_cartesian(ylim=c(0.3, max(t_pi$thetapi_m)), xlim=c(-13,13))
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)


#gamma=+5, h=0
t_pi <- read.table("gamma_pos_five_dom_zero_10.final", h=T)
t_tajima <- read.table("../Theory/pi_reduction_sweeps_Tajimasims_2.txt", h=T)
t_pi$min_se <- t_pi$thetapi_m - t_pi$thetapi_se
t_pi$max_se <- t_pi$thetapi_m + t_pi$thetapi_se

yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))
gg <- ggplot(t_pi) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "darkseagreen3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="blue") + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_5_h0, ymin=1.0-thetapi_m_gamma_5_h0-thetapi_se_gamma_5_h0, ymax=1.0-thetapi_m_gamma_5_h0+thetapi_se_gamma_5_h0), color="dark blue", shape=4, size=1) + coord_cartesian(ylim=c(0.3, max(t_pi$thetapi_m)), xlim=c(-13,13))#+ geom_point(data=t_th, aes(x=rho, y=1.0-thetapi_m), color="dark blue", shape=16, size=3)
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)


#gamma=+5, h=1
t_pi <- read.table("gamma_pos_five_dom_one_10.final", h=T)
t_tajima <- read.table("../Theory/pi_reduction_sweeps_Tajimasims_2.txt", h=T)
t_pi$min_se <- t_pi$thetapi_m - t_pi$thetapi_se
t_pi$max_se <- t_pi$thetapi_m + t_pi$thetapi_se

yname=expression(paste(italic(pi), "/", italic(pi)[0]))
xname=expression(italic(rho))
gg <- ggplot(t_pi) + geom_ribbon(aes(x=rho, ymin = min_se, ymax = max_se), alpha = 0.5, fill = "darkseagreen3", color = "transparent") + geom_line(aes(x=rho, y=thetapi_m), color="blue") + geom_pointrange(data=t_tajima, aes(x=rho, y=1.0-thetapi_m_gamma_5_h1, ymin=1.0-thetapi_m_gamma_5_h1-thetapi_se_gamma_5_h1, ymax=1.0-thetapi_m_gamma_5_h1+thetapi_se_gamma_5_h1), color="dark blue", shape=4, size=1) + coord_cartesian(ylim=c(0.3, max(t_pi$thetapi_m)), xlim=c(-13,13)) #+ geom_point(data=t_th, aes(x=rho, y=1.0-thetapi_m), color="dark blue", shape=16, size=3)
gg1 <- gg + theme_classic() + labs(x=xname, y=yname) + theme(axis.text.x=element_text(size=15),axis.text.y=element_text(size=15),axis.title.x=element_text(size=25),axis.title.y=element_text(size=25), legend.position="right")
gg2 <- gg1 + theme(panel.background = element_rect(fill = "white")) 
print(gg2)
