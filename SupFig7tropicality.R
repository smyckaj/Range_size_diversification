setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)
library(phangorn)
library(diversitree)
library(secsse)
library(parallel)
library(popbio)
library(vioplot)




#prepare data
########################################################
#read mammalian phylogeny
mamphylo=read.tree("vertlife_mammaltree.txt")
mamphylo=force.ultrametric(mamphylo)

#load range data
load("cleanmainlanddf.RData")
load("latdatamammals.RData")

#reshuffle them by species
bakedranges=data.frame(species=mamphylo$tip.label)
bakedranges$order=rangesmergeddf$order_[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$family=rangesmergeddf$family[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$genus=rangesmergeddf$genus[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$area=rangesmergeddf$area[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$mainlandarea=rangesmergeddf$mainlandarea[match(bakedranges$species,rangesmergeddf$binomialul)]
bakedranges$tropicality=latdata$tropicality[match(bakedranges$species,latdata$binomialul)]
bakedranges$islandendemic=as.numeric(bakedranges$mainlandarea==0)
bakedranges$mainlandendemic=as.numeric(bakedranges$mainlandarea==bakedranges$area)
bakedranges=bakedranges[-which(is.na(bakedranges$area)),] #drop species with no range info


#detele tips with no range info, island endemics and chiroptera
tipstodelete=setdiff(mamphylo$tip.label,bakedranges$species)
mamphylo=drop.tip(mamphylo,tipstodelete)

svg("sftropicality.svg", 9/2.54, 9/2.54, pointsize=8)

plot(NULL, xlim=c(0,4),ylim=c(log10(0.0001),log10(1e+08)),xaxt="n",xlab="",ylab="area (km^2)", yaxt="n")


vioplot(log10(bakedranges$area[bakedranges$tropicality=="trop"]), at=1,col="gray30",add=T)
vioplot(log10(bakedranges$area[bakedranges$tropicality=="both"]), at=2,col="gray50",add=T)
vioplot(log10(bakedranges$area[bakedranges$tropicality=="temp"]), at=3,col="gray70",add=T)

axis(2,at=c(-4,-2,0,2,4,6,8),labels=rep("",7))
gplots::angleAxis(2,at=c(-4,-2,0,2,4,6,8),labels=c("0.0001","0.01", "1", "100", "10 000","1 000 000", "100 000 000"))

gplots::angleAxis(1, at=c(1,2,3),labels=c("tropical", "both", "temperate"))

dev.off()
