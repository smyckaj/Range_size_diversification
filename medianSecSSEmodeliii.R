#!/usr/bin/Rscript
setwd("~/Desktop/honza/Diversification_and_range_size")
#setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)
library(phangorn)
library(diversitree)
library(secsse)
library(parallel)
library(popbio)




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


#detele tips with no range info
tipstodelete=setdiff(mamphylo$tip.label,bakedranges$species)
mamphylo=drop.tip(mamphylo,tipstodelete)

#categorize species by median
bakedranges$abovemedian=as.numeric(median(bakedranges$area)<bakedranges$area)

#make secsse all
tipstates=bakedranges$abovemedian+1

num_concealed_states=2
idparslist=cla_id_paramPos(tipstates,num_concealed_states)

lambdas_speltout=list()
lambdas_speltout[[1]]=matrix(0,4,4)
lambdas_speltout[[2]]=matrix(0,4,4)
lambdas_speltout[[3]]=matrix(0,4,4)
lambdas_speltout[[4]]=matrix(0,4,4)
idparslist$lambdas=lambdas_speltout

idparslist$lambdas[[1]][,]=0
idparslist$lambdas[[1]][1,1]=1

idparslist$lambdas[[2]][,]=0
idparslist$lambdas[[2]][2,2]=2
idparslist$lambdas[[2]][2,1]=3

idparslist$lambdas[[3]][,]=0
idparslist$lambdas[[3]][3,3]=4

idparslist$lambdas[[4]][,]=0
idparslist$lambdas[[4]][4,4]=5
idparslist$lambdas[[4]][4,3]=6

idparslist$mus[1]=7
idparslist$mus[2]=8
idparslist$mus[3]=9
idparslist$mus[4]=10

idparslist$Q=idparslist$Q
idparslist$Q[,]=0
diag(idparslist$Q[,])=NA

idparslist$Q [1,2]=11
idparslist$Q [2,1]=12
idparslist$Q [3,4]=13
idparslist$Q [4,3]=14

idparslist$Q [2,4]=15
idparslist$Q [4,2]=16
idparslist$Q [1,3]=15
idparslist$Q [3,1]=16


p=starting.point.classe(mamphylo,2)[c(1,6,5,1,6,5,7,8,7,8,9,10,9,10,10,10)]
fitsecsse=cla_secsse_ml(phy = mamphylo,tipstates, num_concealed_states = num_concealed_states,
                        cond = "proper_cond", root_state_weight = "proper_weights",
                        idparslist = idparslist,
                        idparsopt = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
                        initparsopt= p,
                        idparsfix = 0,
                        parsfix = 0,
                        sampling_fraction =c(1,1),
                        optimmethod = "subplex",
                        tol = c(1e-05, 1e-06, 1e-08))

paste("all together model finished")

#by order
table(bakedranges$order)
orders=c("CARNIVORA", "CETARTIODACTYLA", "CHIROPTERA","EULIPOTYPHLA", "PRIMATES", "RODENTIA")
geosseislandsbyorders=function(order){
  bakedrangesord=bakedranges[which(bakedranges$order==order),]
  mamphyloord=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesord$species))
  
  tipstatesord=bakedrangesord$abovemedian+1
  
  num_concealed_statesord=2
  idparslistord=cla_id_paramPos(tipstatesord,num_concealed_statesord)
  
  lambdas_speltoutord=list()
  lambdas_speltoutord[[1]]=matrix(0,4,4)
  lambdas_speltoutord[[2]]=matrix(0,4,4)
  lambdas_speltoutord[[3]]=matrix(0,4,4)
  lambdas_speltoutord[[4]]=matrix(0,4,4)
  idparslistord$lambdas=lambdas_speltoutord
  
  idparslistord$lambdas[[1]][,]=0
  idparslistord$lambdas[[1]][1,1]=1
  
  idparslistord$lambdas[[2]][,]=0
  idparslistord$lambdas[[2]][2,2]=2
  idparslistord$lambdas[[2]][2,1]=3
  
  idparslistord$lambdas[[3]][,]=0
  idparslistord$lambdas[[3]][3,3]=4
  
  idparslistord$lambdas[[4]][,]=0
  idparslistord$lambdas[[4]][4,4]=5
  idparslistord$lambdas[[4]][4,3]=6
  
  idparslistord$mus[1]=7
  idparslistord$mus[2]=8
  idparslistord$mus[3]=9
  idparslistord$mus[4]=10
  
  idparslistord$Q[,]=0
  diag(idparslistord$Q[,])=NA
  
  idparslistord$Q [1,2]=11
  idparslistord$Q [2,1]=12
  idparslistord$Q [3,4]=13
  idparslistord$Q [4,3]=14
  
  idparslistord$Q [2,4]=15
  idparslistord$Q [4,2]=16
  idparslistord$Q [1,3]=15
  idparslistord$Q [3,1]=16
  
  
  pord=starting.point.classe(mamphyloord,2)[c(1,6,5,1,6,5,7,8,7,8,9,10,9,10,10,10)]
  fitsecsseord=cla_secsse_ml(phy = mamphyloord,tipstatesord, num_concealed_states = num_concealed_statesord,
                             cond = "proper_cond", root_state_weight = "proper_weights",
                          idparslist = idparslistord,
                          idparsopt = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
                          initparsopt= pord,
                          idparsfix = 0,
                          parsfix = 0,
                          sampling_fraction =c(1,1),
                          optimmethod = "subplex",
                          tol = c(1e-05, 1e-06, 1e-08))
  paste(order, "model finished")
  return(fitsecsseord)
}


orderresult=mcmapply(geosseislandsbyorders, orders, SIMPLIFY=F,mc.cores=20)


save.image("medianSecSSEmodeliii.RData")

