setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")
library(coda)
library(BAMMtools)
library(ape)
library(picante)
library(phylolm)
library(gplots)
library(scales)


mcmc=read.csv("BAMMmcmcfinal_out.txt",header = T)
dim(mcmc)

#convergence check
# plot(mcmc$logLik~mcmc$generation)
# plot(mcmc$N_shifts~mcmc$generation)
# 
mcmcbi=mcmc[50000:nrow(mcmc),]
# plot(mcmcbi$logLik~mcmcbi$generation)

effectiveSize(mcmcbi$logLik)
effectiveSize(mcmcbi$N_shifts)


#summary
mamphylo=read.tree("vertlife_mammaltree_ultraforce.txt")
ed=getEventData(mamphylo, "BAMMeventfinal_data.txt", burnin=0.5)

# plot.bammdata(ed, legend=T, spex='netdiv',labels=T)
# plot.bammdata(ed, legend=T, spex='s')
# plot.bammdata(ed, legend=T, spex='e')

# summary(ed)

# tip rates
tr=getTipRates(ed, returnNetDiv = T,statistic = 'mean')
# hist(tr$netdiv.avg)


# rangen size dependence
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
mamphylopr=drop.tip(mamphylo,tipstodelete)

row.names(bakedranges)=bakedranges$species
bakedranges$DR=1/evol.distinct(mamphylopr)$w

bakedranges$log10DR=log10(bakedranges$DR)
bakedranges$log10area=log10(bakedranges$area)

m=match(bakedranges$species,names(tr$netdiv.avg))
bakedranges$BAMM=tr$netdiv.avg[m]

#clean memory
rm(ed)
rm(tr)
gc()

#plots
svg("sfBAMM.svg", 10/2.54, 18/2.54, pointsize=8)

par(mfrow=c(4,2))
par(mar=c(4,4,0.1,0))

modallp=phylolm(BAMM~log10area, data = bakedranges, phy=mamphylopr, model="lambda")
summary(modallp)
modall=phylolm(BAMM~log10area, data = bakedranges, phy=mamphylopr, model="lambda", lower.bound = 0, upper.bound = 0,
               starting.value = 0)
summary(modall)

plot(bakedranges$BAMM~bakedranges$log10area,xaxt="n", yaxt="n", xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modall, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")

axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

plot(NULL)
legend("center",c("mean estimate of phylogenetic model","mean estimate of standard model"),
       col=c("orange", "green"),lty=1,inset=0.0)

#orders
bakedrangescarn=bakedranges[which(bakedranges$order=="CARNIVORA"),]
mamphylocarn=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangescarn$species))

bakedrangesceta=bakedranges[which(bakedranges$order=="CETARTIODACTYLA"),]
mamphyloceta=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesceta$species))

bakedrangeschiro=bakedranges[which(bakedranges$order=="CHIROPTERA"),]
mamphylochiro=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeschiro$species))

bakedrangeseuli=bakedranges[which(bakedranges$order=="EULIPOTYPHLA"),]
mamphyloeuli=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeseuli$species))

bakedrangesprim=bakedranges[which(bakedranges$order=="PRIMATES"),]
mamphyloprim=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesprim$species))

bakedrangesrode=bakedranges[which(bakedranges$order=="RODENTIA"),]
mamphylorode=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesrode$species))

#carnivora
modcarnp=phylolm(BAMM~log10area, data = bakedrangescarn, phy=mamphylocarn, model="lambda")
summary(modcarnp)
modcarn=phylolm(BAMM~log10area, data = bakedrangescarn, phy=mamphylocarn, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0)
summary(modcarn)

plot(bakedrangescarn$BAMM~bakedrangescarn$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))

axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#cetartiodactyla
modcetap=phylolm(BAMM~log10area, data = bakedrangesceta, phy=mamphyloceta, model="lambda")
summary(modcetap)
modceta=phylolm(BAMM~log10area, data = bakedrangesceta, phy=mamphyloceta, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0)
summary(modceta)

plot(bakedrangesceta$BAMM~bakedrangesceta$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))


axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#chiroptera
modchirop=phylolm(BAMM~log10area, data = bakedrangeschiro, phy=mamphylochiro, model="lambda")
summary(modchirop)
modchiro=phylolm(BAMM~log10area, data = bakedrangeschiro, phy=mamphylochiro, model="lambda",lower.bound = 0, upper.bound = 0,
                 starting.value = 0)
summary(modchiro)

plot(bakedrangeschiro$BAMM~bakedrangeschiro$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modchiro, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")


axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#eulipotyphla
modeulip=phylolm(BAMM~log10area, data = bakedrangeseuli, phy=mamphyloeuli, model="lambda")
summary(modeulip)
modeuli=phylolm(BAMM~log10area, data = bakedrangeseuli, phy=mamphyloeuli, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0)
summary(modeuli)

plot(bakedrangeseuli$BAMM~bakedrangeseuli$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))

axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))

#angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#primates
modprimp=phylolm(BAMM~log10area, data = bakedrangesprim, phy=mamphyloprim, model="lambda")
summary(modprimp)
modprim=phylolm(BAMM~log10area, data = bakedrangesprim, phy=mamphyloprim, model="lambda", lower.bound = 0, upper.bound = 0,
                starting.value = 0)
summary(modprim)

plot(bakedrangesprim$BAMM~bakedrangesprim$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modprim, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")

axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))
angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))


#rodentia
modrodep=phylolm(BAMM~log10area, data = bakedrangesrode, phy=mamphylorode, model="lambda")
summary(modrodep)
modrode=phylolm(BAMM~log10area, data = bakedrangesrode, phy=mamphylorode, model="lambda",lower.bound = 0, upper.bound = 0,
                starting.value = 0)
summary(modrode)

plot(bakedrangesrode$BAMM~bakedrangesrode$log10area,xaxt="n", yaxt="n",xlab=parse(text="area (km^2)"),ylab = "BAMM net rate", main="",xlim=c(log10(0.0001),log10(1e+08)),ylim=c(0,1),pch=20, col=alpha(1,0.1))

newdata=data.frame(log10area=seq(-4, 8, 0.01))
ci=predict(modrode, newdata, interval = "confidence")
lines(newdata$log10area,ci[,1], col="green")

axis(2, at=c(0,0.2,0.4,0.6,0.8,1),labels=c(0,0.2,0.4,0.6,0.8,1))
axis(1,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))
angleAxis(1,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

dev.off()

