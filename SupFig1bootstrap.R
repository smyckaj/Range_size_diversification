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



svg("sfbootstrap.svg", 10/2.54, 18/2.54, pointsize=8)

bootpredict=function(bootmatrix, newdata, lty, tr, col){
  for (i in 1:dim(bootmatrix)[1]){
    cib=bootmatrix[i,1]+bootmatrix[i,2]*newdata$log10area
    lines(newdata$log10area,cib, col=alpha(col,tr))
  }
  
}

par(mfrow=c(4,2))

#all together
modallp=phylolm(log10DR~log10area, data = bakedranges, phy=mamphylo, model="lambda", boot = 100)
modall=phylolm(log10DR~log10area, data = bakedranges, phy=mamphylo, model="lambda", lower.bound = 0, upper.bound = 0,
               starting.value = 0, boot = 100)

par(mar=c(4,4,0.1,0))
plot(bakedranges$log10DR~bakedranges$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modall, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modallp, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modall$bootstrap, newdata,1,0.05, col="green")
bootpredict(modallp$bootstrap, newdata,1,0.05, col="orange")




axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

plot(NULL)
legend("center",c("mean estimate of phylogenetic model","100 bootstrap estimates of phylogenetic model",
                  "mean estimate of standard model", "100 bootstrap estimates of standard model"),
       col=c("orange",alpha("orange", 0.3), "green", alpha("green",0.3)),lty=1,inset=0.0)

#carnivora
modcarnp=phylolm(log10DR~log10area, data = bakedrangescarn, phy=mamphylocarn, model="lambda", boot = 100)
modcarn=phylolm(log10DR~log10area, data = bakedrangescarn, phy=mamphylocarn, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0, boot = 100)

plot(bakedrangescarn$log10DR~bakedrangescarn$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modcarn, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modcarnp, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modcarn$bootstrap, newdata,1,0.05, col="green")
bootpredict(modcarnp$bootstrap, newdata,1,0.05, col="orange")

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#cetartiodactyla
modcetap=phylolm(log10DR~log10area, data = bakedrangesceta, phy=mamphyloceta, model="lambda", boot = 100)
modceta=phylolm(log10DR~log10area, data = bakedrangesceta, phy=mamphyloceta, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0, boot = 100)

plot(bakedrangesceta$log10DR~bakedrangesceta$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modceta, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modcetap, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modceta$bootstrap, newdata,1,0.05, col="green")
bootpredict(modcetap$bootstrap, newdata,1,0.05, col="orange")

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#chiroptera
modchirop=phylolm(log10DR~log10area, data = bakedrangeschiro, phy=mamphylochiro, model="lambda", boot = 100)
modchiro=phylolm(log10DR~log10area, data = bakedrangeschiro, phy=mamphylochiro, model="lambda", lower.bound = 0, upper.bound = 0,
                 starting.value = 0, boot = 100)

plot(bakedrangeschiro$log10DR~bakedrangeschiro$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modchiro, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modchirop, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modchiro$bootstrap, newdata,1,0.05, col="green")
bootpredict(modchirop$bootstrap, newdata,1,0.05, col="orange")

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#eulipotyphla
modeulip=phylolm(log10DR~log10area, data = bakedrangeseuli, phy=mamphyloeuli, model="lambda", boot = 100)
modeuli=phylolm(log10DR~log10area, data = bakedrangeseuli, phy=mamphyloeuli, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0, boot = 100)

plot(bakedrangeseuli$log10DR~bakedrangeseuli$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modeuli, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modeulip, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modeuli$bootstrap, newdata,1,0.05, col="green")
bootpredict(modeulip$bootstrap, newdata,1,0.05, col="orange")

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#primates
modprimp=phylolm(log10DR~log10area, data = bakedrangesprim, phy=mamphyloprim, model="lambda", boot = 100)
modprim=phylolm(log10DR~log10area, data = bakedrangesprim, phy=mamphyloprim, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0, boot = 100)

plot(bakedrangesprim$log10DR~bakedrangesprim$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modprim, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modprimp, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modprim$bootstrap, newdata,1,0.05, col="green")
bootpredict(modprimp$bootstrap, newdata,1,0.05, col="orange")

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

#rodentia
modrodep=phylolm(log10DR~log10area, data = bakedrangesrode, phy=mamphylorode, model="lambda", boot = 100)
modrode=phylolm(log10DR~log10area, data = bakedrangesrode, phy=mamphylorode, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0, boot = 100)

plot(bakedrangesrode$log10DR~bakedrangesrode$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "DR metric", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(log10(0.01),log10(10)),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modrode, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")
cip=predict(modrodep, newdata, interval = "confidence")
lines(newdata$log10area,cip[,1],col="orange")

bootpredict(modrode$bootstrap, newdata,1,0.05, col="green")
bootpredict(modrodep$bootstrap, newdata,1,0.05, col="orange")

axis(2, at=c(-2,-1,0,1),labels=c(0.01,0.1,1,10))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))



dev.off()