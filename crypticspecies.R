#!/usr/bin/Rscript
#setwd("~/Desktop/honza/Diversification_and_range_size")
setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)
library(phangorn)
library(diversitree)
library(secsse)
library(parallel)
library(popbio)
library(ade4)

source("./stickTips.R")
source("./stickTips_1.2.r")


set.seed(10)


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
bakedranges=bakedranges[-which(is.na(bakedranges$area)),] #drop species with no range info


#delete tips with no range info
tipstodelete=setdiff(mamphylo$tip.label,bakedranges$species)
mamphylo=drop.tip(mamphylo,tipstodelete)

#ordered bakedranges
bakedrangessort=bakedranges[order(bakedranges$area,decreasing = T),]

### 1%
perc1=round(length(mamphylo$tip.label)/100)
to_be_splitted1=as.character(bakedrangessort$species[1:perc1])

##phylogeny
cortable1=data.frame(new.tips=as.character(mamphylo$tip.label), old.tips=as.character(mamphylo$tip.label))
cortable1$new.tips=as.character(cortable1$new.tips)
cortable1$old.tips=as.character(cortable1$old.tips)
m=match(to_be_splitted1,cortable1$new.tips)
cortable1$new.tips[m]=paste(cortable1$new.tips[m], "1",sep = "_")

newcomerct=data.frame(new.tips=paste(to_be_splitted1, "2", sep="_"), old.tips=to_be_splitted1)
newcomerct$new.tips=as.character(newcomerct$new.tips)
newcomerct$old.tips=as.character(newcomerct$old.tips)

cortable1=rbind(cortable1, newcomerct)

#polytomy hardwires the split in the middle
mamphylo1=stickTips(mamphylo,cortable1, subtree.type = "polytomy")

#check
plot(extract.clade(mamphylo,node = getMRCA(mamphylo,c("Canis_lupus","Vulpes_corsac"))))
plot(extract.clade(mamphylo1,node = getMRCA(mamphylo1,c("Canis_lupus_1","Vulpes_corsac"))))

##table
bakedranges1=bakedranges
bakedranges1$species=as.character(bakedranges1$species)
n=match(to_be_splitted1,bakedranges1$species)
bakedranges1$species[n]=paste(bakedranges1$species[n], "1", sep="_")
bakedranges1$area[n]=runif(n,0,bakedranges$area[n])

newcomerdf=bakedranges[n,]
newcomerdf$species=paste(as.character(newcomerdf$species), "2", sep="_")
newcomerdf$area=bakedranges$area[n]-bakedranges1$area[n]

bakedranges1=rbind(bakedranges1, newcomerdf)
bakedranges1=bakedranges1[match(mamphylo1$tip.label, bakedranges1$species),]

save(mamphylo1,bakedranges1,file="cryptic1percent.RData")



### 5%
perc5=round(length(mamphylo$tip.label)/20)
to_be_splitted5=as.character(bakedrangessort$species[1:perc5])

##phylogeny
cortable5=data.frame(new.tips=as.character(mamphylo$tip.label), old.tips=as.character(mamphylo$tip.label))
cortable5$new.tips=as.character(cortable5$new.tips)
cortable5$old.tips=as.character(cortable5$old.tips)
m=match(to_be_splitted5,cortable5$new.tips)
cortable5$new.tips[m]=paste(cortable5$new.tips[m], "1",sep = "_")

newcomerct=data.frame(new.tips=paste(to_be_splitted5, "2", sep="_"), old.tips=to_be_splitted5)
newcomerct$new.tips=as.character(newcomerct$new.tips)
newcomerct$old.tips=as.character(newcomerct$old.tips)

cortable5=rbind(cortable5, newcomerct)

#polytomy hardwires the split in the middle
mamphylo5=stickTips(mamphylo,cortable5, subtree.type = "polytomy")


##table
bakedranges5=bakedranges
bakedranges5$species=as.character(bakedranges5$species)
n=match(to_be_splitted5,bakedranges5$species)
bakedranges5$species[n]=paste(bakedranges5$species[n], "1", sep="_")
bakedranges5$area[n]=runif(n,0,bakedranges$area[n])

newcomerdf=bakedranges[n,]
newcomerdf$species=paste(as.character(newcomerdf$species), "2", sep="_")
newcomerdf$area=bakedranges$area[n]-bakedranges5$area[n]

bakedranges5=rbind(bakedranges5, newcomerdf)
bakedranges5=bakedranges5[match(mamphylo5$tip.label, bakedranges5$species),]

save(mamphylo5,bakedranges5,file="cryptic5percent.RData")


### 10% 
perc10=round(length(mamphylo$tip.label)/10)
to_be_splitted10=as.character(bakedrangessort$species[1:perc10])

##phylogeny
cortable10=data.frame(new.tips=as.character(mamphylo$tip.label), old.tips=as.character(mamphylo$tip.label))
cortable10$new.tips=as.character(cortable10$new.tips)
cortable10$old.tips=as.character(cortable10$old.tips)
m=match(to_be_splitted10,cortable10$new.tips)
cortable10$new.tips[m]=paste(cortable10$new.tips[m], "1",sep = "_")

newcomerct=data.frame(new.tips=paste(to_be_splitted10, "2", sep="_"), old.tips=to_be_splitted10)
newcomerct$new.tips=as.character(newcomerct$new.tips)
newcomerct$old.tips=as.character(newcomerct$old.tips)

cortable10=rbind(cortable10, newcomerct)

#polytomy hardwires the split in the middle
mamphylo10=stickTips(mamphylo,cortable10, subtree.type = "polytomy")


##table
bakedranges10=bakedranges
bakedranges10$species=as.character(bakedranges10$species)
n=match(to_be_splitted10,bakedranges10$species)
bakedranges10$species[n]=paste(bakedranges10$species[n], "1", sep="_")
bakedranges10$area[n]=runif(n,0,bakedranges$area[n])

newcomerdf=bakedranges[n,]
newcomerdf$species=paste(as.character(newcomerdf$species), "2", sep="_")
newcomerdf$area=bakedranges$area[n]-bakedranges10$area[n]

bakedranges10=rbind(bakedranges10, newcomerdf)
bakedranges10=bakedranges10[match(mamphylo10$tip.label, bakedranges10$species),]

save(mamphylo10,bakedranges10,file="cryptic10percent.RData")


### 30% 
perc30=round(3*length(mamphylo$tip.label)/10)
to_be_splitted30=as.character(bakedrangessort$species[1:perc30])

##phylogeny
cortable30=data.frame(new.tips=as.character(mamphylo$tip.label), old.tips=as.character(mamphylo$tip.label))
cortable30$new.tips=as.character(cortable30$new.tips)
cortable30$old.tips=as.character(cortable30$old.tips)
m=match(to_be_splitted30,cortable30$new.tips)
cortable30$new.tips[m]=paste(cortable30$new.tips[m], "1",sep = "_")

newcomerct=data.frame(new.tips=paste(to_be_splitted30, "2", sep="_"), old.tips=to_be_splitted30)
newcomerct$new.tips=as.character(newcomerct$new.tips)
newcomerct$old.tips=as.character(newcomerct$old.tips)

cortable30=rbind(cortable30, newcomerct)

#polytomy hardwires the split in the middle
mamphylo30=stickTips(mamphylo,cortable30, subtree.type = "polytomy")

#check
plot(extract.clade(mamphylo,node = getMRCA(mamphylo,c("Canis_lupus","Vulpes_corsac"))))
plot(extract.clade(mamphylo30,node = getMRCA(mamphylo30,c("Canis_lupus_2","Vulpes_corsac_1"))))

##table
bakedranges30=bakedranges
bakedranges30$species=as.character(bakedranges30$species)
n=match(to_be_splitted30,bakedranges30$species)
bakedranges30$species[n]=paste(bakedranges30$species[n], "1", sep="_")
bakedranges30$area[n]=runif(n,0,bakedranges$area[n])

newcomerdf=bakedranges[n,]
newcomerdf$species=paste(as.character(newcomerdf$species), "2", sep="_")
newcomerdf$area=bakedranges$area[n]-bakedranges30$area[n]

bakedranges30=rbind(bakedranges30, newcomerdf)
bakedranges30=bakedranges30[match(mamphylo30$tip.label, bakedranges30$species),]

save(mamphylo30,bakedranges30,file="cryptic30percent.RData")
