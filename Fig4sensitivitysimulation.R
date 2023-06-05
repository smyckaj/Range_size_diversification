#!/usr/bin/Rscript
#setwd("~/Desktop/honza/Diversification_and_range_size")
setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)
library(phangorn)
library(diversitree)
library(secsse)
library(parallel)
library(picante)
library(phylolm)
library(vioplot)
library(scales)
library(dplyr)




load("medianSecSSEmodeliii.RData")
row.names(bakedranges)=bakedranges$species


#real data
bakedranges$DR=1/evol.distinct(mamphylo)$w

vioplot(log10(bakedranges$DR)~bakedranges$abovemedian,xaxt="n",ylim=c(-2,1),  xlab="", ylab="log(DR metric)")

linmod=lm(log10(bakedranges$DR)~bakedranges$abovemedian)
summary(linmod)
abline(linmod,lwd=3, col="orange")

linmodp=phylolm(log10(DR)~abovemedian, data = bakedranges, phy=mamphylo, model="lambda")
summary(linmodp)
abline(linmodp,lwd=3, col="green")

legend("topleft",c("standard linear", "phylogenetically corrected"), pch=19, col=c("orange","green"),inset=0.0)
gplots::angleAxis(1,at=1:2,labels=c("small-ranged", "large-ranged"))




#simulated data
parlist=starting.point.classe(mamphylo,4)
parlist[1:length(parlist)]=rep(0, length(parlist))
fitsecsse
parlist[names(parlist)=="lambda111"]=fitsecsse$MLpars[[1]][[1]][1,1]
parlist[names(parlist)=="lambda212"]=fitsecsse$MLpars[[1]][[2]][2,1]
parlist[names(parlist)=="lambda222"]=fitsecsse$MLpars[[1]][[2]][2,2]
parlist[names(parlist)=="lambda333"]=fitsecsse$MLpars[[1]][[3]][3,3]
parlist[names(parlist)=="lambda434"]=fitsecsse$MLpars[[1]][[4]][4,3]
parlist[names(parlist)=="lambda444"]=fitsecsse$MLpars[[1]][[4]][4,4]
parlist[names(parlist)=="mu1"]=fitsecsse$MLpars[[2]][1]
parlist[names(parlist)=="mu2"]=fitsecsse$MLpars[[2]][2]
parlist[names(parlist)=="mu3"]=fitsecsse$MLpars[[2]][3]
parlist[names(parlist)=="mu4"]=fitsecsse$MLpars[[2]][4]
parlist[names(parlist)=="q12"]=fitsecsse$MLpars[[3]][1,2]
parlist[names(parlist)=="q13"]=fitsecsse$MLpars[[3]][1,3]
parlist[names(parlist)=="q21"]=fitsecsse$MLpars[[3]][2,1]
parlist[names(parlist)=="q24"]=fitsecsse$MLpars[[3]][2,4]
parlist[names(parlist)=="q31"]=fitsecsse$MLpars[[3]][3,1]
parlist[names(parlist)=="q34"]=fitsecsse$MLpars[[3]][3,4]
parlist[names(parlist)=="q42"]=fitsecsse$MLpars[[3]][4,2]
parlist[names(parlist)=="q43"]=fitsecsse$MLpars[[3]][4,3]

#one test
repeat{
    simtr=tree.classe(parlist , max.taxa = 5000, max.t = 200)
    if (inherits(simtr, "phylo")) 
        break
}
plot(simtr)
axisPhylo()

simbakedranges=data.frame(abovemedian=as.numeric(!(simtr$tip.state%%2)), DR=1/evol.distinct(simtr)$w)
row.names(simbakedranges)=names(simtr$tip.state)


vioplot(log10(simbakedranges$DR)~simbakedranges$abovemedian)

simlinmod=lm(log10(simbakedranges$DR)~simbakedranges$abovemedian)
summary(simlinmod)
abline(simlinmod,lwd=3)
simlinmod$coefficients[1]
simlinmod$coefficients[2]
summary(simlinmod)$coefficients[2,4]

simlinmodp=phylolm(log10(DR)~abovemedian, data = simbakedranges, phy=simtr, model="lambda")
summary(simlinmodp)
abline(simlinmodp,lwd=3)
simlinmodp$coefficients[1]
simlinmodp$coefficients[2]
summary(simlinmodp)$coefficients[2,4]



#100 datasets
simtrees=list()
simdata=list()
simintercept=rep(0,100)
siminterceptp=rep(0,100)
simslope=rep(0,100)
simslopep=rep(0,100)
simpval=rep(0,100)
simvalp=rep(0,100)

plot(NULL, xlim=c(0,1), ylim=c(-1.5,0))

for(i in 1:100){
    repeat{
        simtr=tree.classe(parlist , max.taxa = 5000, max.t = 200)
        if (inherits(simtr, "phylo")) 
            break
    }
    simtrees[[i]]=simtr
    
    simbakedranges=data.frame(abovemedian=as.numeric(!(simtr$tip.state%%2)), DR=1/evol.distinct(simtr)$w)
    row.names(simbakedranges)=names(simtr$tip.state)
    simdata[[i]]=simbakedranges
    
    simlinmod=lm(log10(simbakedranges$DR)~simbakedranges$abovemedian)
    simintercept[i]=simlinmod$coefficients[1]
    simslope[i]=simlinmod$coefficients[2]
    simpval[i]=summary(simlinmod)$coefficients[2,4]
    abline(simlinmod, col="orange")
    
    simlinmodp=phylolm(log10(DR)~abovemedian, data = simbakedranges, phy=simtr, model="lambda")
    siminterceptp[i]=simlinmodp$coefficients[1]
    simslopep[i]=simlinmodp$coefficients[2]
    simvalp[i]=summary(simlinmodp)$coefficients[2,4]
    abline(simlinmodp, col="green")
    
    
}

#regression slope plot
vioplot(simslope, simslopep, col=c(alpha("orange",0.3), alpha("green",0.3)), xaxt="n", ylab="regression slope")
points(1:2, c(linmod$coefficients[2],linmodp$coefficients[2]),col=c("orange","green"),pch=18,cex=3 )

# abline(linmod$coefficients[2],0, col=alpha("orange",0.5), lwd=3)
# abline(linmodp$coefficients[2],0, col=alpha("green",0.5), lwd=3)
legend("topleft",c("standard linear", "phylogenetically corrected"), pch=19, col=c("orange","green"),inset=0.0)


#p value plot
vioplot(simpval, simvalp, col=c(alpha("orange",0.3), alpha("green",0.3)), xaxt="n", ylab="p value")
points(1:2, c(summary(linmod)$coefficients[2,4],summary(linmodp)$coefficients[2,4]),col=c("orange","green"),pch=18,cex=3 )

# abline(summary(linmod)$coefficients[2,4],0, col=alpha("orange",0.5), lwd=3)
# abline(summary(linmodp)$coefficients[2,4],0, col=alpha("green",0.5), lwd=3)
legend("topleft",c("standard linear", "phylogenetically corrected"), pch=19, col=c("orange","green"),inset=0.0)

#lines plot
vioplot(log10(bakedranges$DR)~bakedranges$abovemedian,xaxt="n",ylim=c(-2,1),  xlab="", ylab="log(DR metric)")
abline(linmod,lwd=3, col="orange")
abline(linmodp,lwd=3, col="green")

for (i in 1:100){
    abline(simintercept[i],simslope[i], col=alpha("orange",0.05))
    abline(siminterceptp[i],simslopep[i], col=alpha("green",0.05))
}

legend("topleft",c("standard linear", "phylogenetically corrected"), pch=19, col=c("orange","green"),inset=0.0)
gplots::angleAxis(1,at=1:2,labels=c("small-ranged", "large-ranged"))


#examples
sam=sample(1:100,8)

par(mfrow=c(3,3))
vioplot(log10(bakedranges$DR)~bakedranges$abovemedian,xaxt="n",ylim=c(-2,1),  xlab="", ylab="log(DR metric)", main="Real")
abline(linmod,lwd=3, col="orange")
abline(linmodp,lwd=3, col="green")
gplots::angleAxis(1,at=1:2,labels=c("small-ranged", "large-ranged"))

for (j in sam){
    vioplot(log10(simdata[[j]]$DR)~simdata[[j]]$abovemedian,xaxt="n",ylim=c(-2,1),  xlab="", ylab="log(DR metric)", main=paste("Simulation",j))
    abline(simintercept[i],simslope[i],lwd=3, col="orange")
    abline(siminterceptp[i],simslopep[i],lwd=3, col="green")
    gplots::angleAxis(1,at=1:2,labels=c("small-ranged", "large-ranged"))
}

#save.image("~/Dropbox/CTS postdoc/Diversification_and_range_size/Fig4sensitivitysimulation.RData")

svg("fig4.svg", 10.5/2.54, 9/2.54, pointsize=8)
par(mfrow=c(1,2))
vioplot(simslope, col=alpha("black",0.3),ylim=c(-0.1,0),xlim=c(0.5,2.5), yaxt="n", xaxt="n", ylab="regression slope", main="")
vioplot(at=2,simslopep, col=alpha("black",0.3), lty=1, add=T)

points(1:2, c(linmod$coefficients[2],linmodp$coefficients[2]),pch=4,cex=5 )
abline(0,0,lty=2)
axis(side=2, at=c(-0.1, -0.05, 0))
gplots::angleAxis(1,at=1:2,labels=c("standard linear", "phylo. corr."))

#legend("topleft",c("standard linear", "phylogenetically corrected"), lty=c(1,2), lwd=2, inset=0.0)

#p value plot
vioplot(simpval, col=alpha("black",0.3),ylim=c(0,0.1),xlim=c(0.5,2.5), yaxt="n", xaxt="n", ylab="p value", main="")
vioplot(at=2,simvalp, col=alpha("black",0.3), lty=1, add=T)

points(1:2, c(summary(linmod)$coefficients[2,4],summary(linmodp)$coefficients[2,4]),pch=4,cex=5 )
abline(0.05,0,lty=2)
axis(side=2, at=c(0, 0.05,  0.1))
gplots::angleAxis(1,at=1:2,labels=c("standard linear", "phylo. corr."))
dev.off()



