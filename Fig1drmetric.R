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
library(scales)
library(phylolm)
library(rphylopic)
library(gplots)




#prepare data
########################################################
#read mammalian phylogeny
mamphylo=read.tree("vertlife_mammaltree.txt")
mamphylo=force.ultrametric(mamphylo)

#load range data
load("cleanmainlanddf.RData")

#reshuffle them by species
bakedranges=data.frame(species=mamphylo$tip.label)
bakedranges$order=rangesmergeddf$order_[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$family=rangesmergeddf$family[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$genus=rangesmergeddf$genus[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$area=rangesmergeddf$area[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$mainlandarea=rangesmergeddf$mainlandarea[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$islandendemic=as.numeric(bakedranges$mainlandarea==0)
bakedranges$mainlandendemic=as.numeric(bakedranges$mainlandarea==bakedranges$area)
bakedranges=bakedranges[-which(is.na(bakedranges$area)),] #drop species with no range info

tipstodelete=setdiff(mamphylo$tip.label,bakedranges$species)
mamphylo=drop.tip(mamphylo,tipstodelete)

row.names(bakedranges)=bakedranges$species

#terminal branch length and DR calculation
tbl=function(tree){
  n<-length(tree$tip.label)
  ee<-setNames(tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])],tree$tip.label)
  return(ee)
}


bakedranges$DR=1/evol.distinct(mamphylo)$w
bakedranges$tbl=tbl(mamphylo)

# log values
bakedranges$log10DR=log10(bakedranges$DR)
bakedranges$log10tbl=log10(bakedranges$tbl)
bakedranges$log10area=log10(bakedranges$area)


bakedrangescarn=bakedranges[which(bakedranges$order=="CARNIVORA"),]
mamphylocarn=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangescarn$species))
bakedrangescarn$DR=1/evol.distinct(mamphylocarn)$w
bakedrangescarn$tbl=tbl(mamphylocarn)
bakedrangescarn$log10DR=log10(bakedrangescarn$DR)
bakedrangescarn$log10tbl=log10(bakedrangescarn$tbl)

bakedrangesceta=bakedranges[which(bakedranges$order=="CETARTIODACTYLA"),]
mamphyloceta=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesceta$species))
bakedrangesceta$DR=1/evol.distinct(mamphyloceta)$w
bakedrangesceta$tbl=tbl(mamphyloceta)
bakedrangesceta$log10DR=log10(bakedrangesceta$DR)
bakedrangesceta$log10tbl=log10(bakedrangesceta$tbl)

bakedrangeschiro=bakedranges[which(bakedranges$order=="CHIROPTERA"),]
mamphylochiro=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeschiro$species))
bakedrangeschiro$DR=1/evol.distinct(mamphylochiro)$w
bakedrangeschiro$tbl=tbl(mamphylochiro)
bakedrangeschiro$log10DR=log10(bakedrangeschiro$DR)
bakedrangeschiro$log10tbl=log10(bakedrangeschiro$tbl)

bakedrangeseuli=bakedranges[which(bakedranges$order=="EULIPOTYPHLA"),]
mamphyloeuli=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeseuli$species))
bakedrangeseuli$DR=1/evol.distinct(mamphyloeuli)$w
bakedrangeseuli$tbl=tbl(mamphyloeuli)
bakedrangeseuli$log10DR=log10(bakedrangeseuli$DR)
bakedrangeseuli$log10tbl=log10(bakedrangeseuli$tbl)

bakedrangesprim=bakedranges[which(bakedranges$order=="PRIMATES"),]
mamphyloprim=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesprim$species))
bakedrangesprim$DR=1/evol.distinct(mamphyloprim)$w
bakedrangesprim$tbl=tbl(mamphyloprim)
bakedrangesprim$log10DR=log10(bakedrangesprim$DR)
bakedrangesprim$log10tbl=log10(bakedrangesprim$tbl)

bakedrangesrode=bakedranges[which(bakedranges$order=="RODENTIA"),]
mamphylorode=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesrode$species))
bakedrangesrode$DR=1/evol.distinct(mamphylorode)$w
bakedrangesrode$tbl=tbl(mamphylorode)
bakedrangesrode$log10DR=log10(bakedrangesrode$DR)
bakedrangesrode$log10tbl=log10(bakedrangesrode$tbl)

#plots
svg("fig1.svg", 10/2.54, 18/2.54, pointsize=8)

par(mfrow=c(4,2))

#all together
modallp=phylolm(log10DR~log10area, data = bakedranges, phy=mamphylo, model="lambda")
summary(modallp)
modall=lm(log10DR~log10area, data = bakedranges)
summary(modall)

par(mar=c(4,4,0.1,0))
plot(bakedranges$log10DR~bakedranges$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modall, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" < 0.001"),adj=0)
text(-4,0.75,expression("p"[phy]*" < 0.001"),adj=0)

plot(NULL)
legend("center",c("mean estimate of standard linear model","95% confidence intervals"),lty=c(1,2),inset=0.0)

#carnivora
modcarnp=phylolm(log10DR~log10area, data = bakedrangescarn, phy=mamphylocarn, model="lambda")
summary(modcarnp)
modcarn=lm(log10DR~log10area, data = bakedrangescarn)
summary(modcarn)

plot(bakedrangescarn$log10DR~bakedrangescarn$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modcarn, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" = 0.613"),adj=0)
text(-4,0.75,expression("p"[phy]*" = 0.109"),adj=0)

#cetartiodactyla
modcetap=phylolm(log10DR~log10area, data = bakedrangesceta, phy=mamphyloceta, model="lambda")
summary(modcetap)
modceta=lm(log10DR~log10area, data = bakedrangesceta)
summary(modceta)

plot(bakedrangesceta$log10DR~bakedrangesceta$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modceta, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" = 0.003"),adj=0)
text(-4,0.75,expression("p"[phy]*" = 0.045"),adj=0)

#chiroptera
modchirop=phylolm(log10DR~log10area, data = bakedrangeschiro, phy=mamphylochiro, model="lambda")
summary(modchirop)
modchiro=lm(log10DR~log10area, data = bakedrangeschiro)
summary(modchiro)

plot(bakedrangeschiro$log10DR~bakedrangeschiro$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modchiro, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" < 0.001"),adj=0)
text(-4,0.75,expression("p"[phy]*" = 0.037"),adj=0)

#eulipotyphla
modeulip=phylolm(log10DR~log10area, data = bakedrangeseuli, phy=mamphyloeuli, model="lambda")
summary(modeulip)
modeuli=lm(log10DR~log10area, data = bakedrangeseuli)
summary(modeuli)

plot(bakedrangeseuli$log10DR~bakedrangeseuli$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modeuli, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" = 0.023"),adj=0)
text(-4,0.75,expression("p"[phy]*" = 0.182"),adj=0)

#primates
modprimp=phylolm(log10DR~log10area, data = bakedrangesprim, phy=mamphyloprim, model="lambda")
summary(modprimp)
modprim=lm(log10DR~log10area, data = bakedrangesprim)
summary(modprim)

plot(bakedrangesprim$log10DR~bakedrangesprim$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modprim, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" = 0.010"),adj=0)
text(-4,0.75,expression("p"[phy]*" = 0.126"),adj=0)

#rodentia
modrodep=phylolm(log10DR~log10area, data = bakedrangesrode, phy=mamphylorode, model="lambda")
summary(modrodep)
modrode=lm(log10DR~log10area, data = bakedrangesrode)
summary(modrode)

plot(bakedrangesrode$log10DR~bakedrangesrode$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modrode, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1])
lines(newdata$log10area,ci[,2],lty=2)
lines(newdata$log10area,ci[,3],lty=2)

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

text(-4,1,expression("p"[std]*" < 0.001"), adj=0)
text(-4,0.75,expression("p"[phy]*" = 0.174"), adj=0)


dev.off()




