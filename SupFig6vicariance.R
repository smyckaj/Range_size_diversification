#!/usr/bin/Rscript
#setwd("~/Desktop/honza/Diversification_and_range_size")
setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)

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

#simulated vicariances
reps=10000

eventsbs=rep(0,reps)
for (i in 1:reps){
  size=sample(bakedranges$area,1)
  a=runif(1,0,size)
  b=size-a
  if(size>median(bakedranges$area)){
    if(a>median(bakedranges$area) & b>median(bakedranges$area)){eventsbs[i]="LtoLL"}
    if(a<median(bakedranges$area) & b>median(bakedranges$area)){eventsbs[i]="LtoLS"}
    if(a>median(bakedranges$area) & b<median(bakedranges$area)){eventsbs[i]="LtoLS"}
    if(a<median(bakedranges$area) & b<median(bakedranges$area)){eventsbs[i]="LtoSS"}
  }
  if(size<=median(bakedranges$area)){eventsbs[i]="StoSS"}
}

tabbs=table(eventsbs)
tabbs

eventshf=rep(0,reps)
for (i in 1:reps){
  size=sample(bakedranges$area,1)
  a=size*0.5
  b=size-a
  if(size>median(bakedranges$area)){
    if(a>median(bakedranges$area) & b>median(bakedranges$area)){eventshf[i]="LtoLL"}
    if(a<median(bakedranges$area) & b>median(bakedranges$area)){eventshf[i]="LtoLS"}
    if(a>median(bakedranges$area) & b<median(bakedranges$area)){eventshf[i]="LtoLS"}
    if(a<median(bakedranges$area) & b<median(bakedranges$area)){eventshf[i]="LtoSS"}
  }
  if(size<=median(bakedranges$area)){eventshf[i]="StoSS"}
}

tabhf=table(eventshf)
tabhf=append(tabhf,0, after=1)

eventsgp=rep(0,reps)
probs=sqrt(bakedranges$area)/sum(sqrt(bakedranges$area))
for (i in 1:reps){
  size=sample(bakedranges$area,1,prob=probs)
  a=runif(1,0,size)
  b=size-a
  if(size>median(bakedranges$area)){
    if(a>median(bakedranges$area) & b>median(bakedranges$area)){eventsgp[i]="LtoLL"}
    if(a<median(bakedranges$area) & b>median(bakedranges$area)){eventsgp[i]="LtoLS"}
    if(a>median(bakedranges$area) & b<median(bakedranges$area)){eventsgp[i]="LtoLS"}
    if(a<median(bakedranges$area) & b<median(bakedranges$area)){eventsgp[i]="LtoSS"}
  }
  if(size<=median(bakedranges$area)){eventsgp[i]="StoSS"}
}

tabgp=table(eventsgp)
tabgp

#plots
svg("sfvicariance.svg", 9/2.54, 10/2.54, pointsize=8)

barplot(t(rbind(rev(tabbs)/reps, rev(tabhf)/reps, rev(tabgp)/reps)),ylab="proportion of simulations", col=c("grey30", "grey50", "grey70", "grey90") )
gplots::angleAxis(1,at=c(0.7,1.9,3.1),labels=c("broken stick", "halves", "size dependent"))
legend("topright",c("large -> large, large","large -> large, small", "large -> small, small", "small -> small, small"), col=c("grey90", "grey70", "grey50", "grey30"),pch=19,inset=0.0)

dev.off()













