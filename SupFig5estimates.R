#!/usr/bin/Rscript
#setwd("~/Desktop/honza/Diversification_and_range_size")
setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)
library(phangorn)
library(diversitree)
library(secsse)
library(parallel)

#load secsse estimates
load("medianSecSSEmodeliii.RData")

#png("figSF2.png",3000,1500,res=300)
svg("sfestimates.svg", 15/2.54, 18/2.54, pointsize=10)
par(mfrow=c(4,2))
par(mar = c(0.1, 4.1, 4.1, 2.1))

a=2
b=6
c=2

plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])
plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="all together", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

par(mar = c(5, 4.1, 4.1, 2.1))
plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="all together", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
gplots::angleAxis(1,at=1:5,labels=c("lambdaL", "lambdaLS", "lambdaS", "muL", "muS"))
legend("topright",c("regime 1", "regime 2"), pch=c(a,b),inset=0.05,pt.cex=c)
par(mar = c(0.1, 4.1, 4.1, 2.1))


fitsecsse=orderresult[[1]]
plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])



plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="Carnivora", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

fitsecsse=orderresult[[2]]
plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])

plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="Cetartiodactyla", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

fitsecsse=orderresult[[3]]
plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])

plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="Chiroptera", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

fitsecsse=orderresult[[4]]
plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])

plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="Eulipotyphla", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

fitsecsse=orderresult[[5]]
plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])

plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="Primates", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

fitsecsse=orderresult[[6]]
plotvect=c(fitsecsse$MLpars[[1]][[2]][2,2], fitsecsse$MLpars[[1]][[2]][2,1], fitsecsse$MLpars[[1]][[1]][1,1],
           fitsecsse$MLpars[[1]][[4]][4,4], fitsecsse$MLpars[[1]][[4]][4,3], fitsecsse$MLpars[[1]][[3]][3,3],
           fitsecsse$MLpars[[2]][2], fitsecsse$MLpars[[2]][1], fitsecsse$MLpars[[2]][4], fitsecsse$MLpars[[2]][3],
           fitsecsse$MLpars[[3]][2,1], fitsecsse$MLpars[[3]][1,2], fitsecsse$MLpars[[3]][4,3], fitsecsse$MLpars[[3]][3,4],
           fitsecsse$MLpars[[3]][3,1], fitsecsse$MLpars[[3]][1,3])

plot(c(1,2,3,1,2,3,4,5,4,5),plotvect[1:10], 
     xlim=c(0,6),ylim=c(0,1.1),xaxt="n",xlab="",ylab="rate", main="Rodentia", pch=c(a,a,a,b,b,b,a,a,b,b),cex=c)
abline(h=c(0,0.2,0.4,0.6,0.8,1),col = "lightgray", lty = "dotted")
abline(v=1:5,col = "lightgray", lty = "dotted")

dev.off()

