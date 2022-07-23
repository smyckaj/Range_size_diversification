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

num_concealed_states=4
idparslist=id_paramPos(tipstates,num_concealed_states)

idparslist$lambdas[1]=1
idparslist$lambdas[2]=1
idparslist$lambdas[3]=2
idparslist$lambdas[4]=2
idparslist$lambdas[5]=3
idparslist$lambdas[6]=3
idparslist$lambdas[7]=4
idparslist$lambdas[8]=4

idparslist$mus[1]=5
idparslist$mus[2]=5
idparslist$mus[3]=6
idparslist$mus[4]=6
idparslist$mus[5]=7
idparslist$mus[6]=7
idparslist$mus[7]=8
idparslist$mus[8]=8

idparslist$Q=idparslist$Q
idparslist$Q[,]=0
diag(idparslist$Q[,])=NA

idparslist$Q [1,2]=9
idparslist$Q [2,1]=10
idparslist$Q [3,4]=9
idparslist$Q [4,3]=10
idparslist$Q [5,6]=9
idparslist$Q [6,5]=10
idparslist$Q [7,8]=9
idparslist$Q [8,7]=10

#AB
idparslist$Q [2,4]=11
idparslist$Q [1,3]=11
idparslist$Q [4,2]=12
idparslist$Q [3,1]=12

#AC
idparslist$Q [2,6]=13
idparslist$Q [1,5]=13
idparslist$Q [6,2]=14
idparslist$Q [5,1]=14

#AD
idparslist$Q [2,8]=15
idparslist$Q [1,7]=15
idparslist$Q [8,2]=16
idparslist$Q [7,1]=16

#BC
idparslist$Q [4,6]=17
idparslist$Q [3,5]=17
idparslist$Q [6,4]=18
idparslist$Q [5,3]=18

#BD
idparslist$Q [4,8]=19
idparslist$Q [3,7]=19
idparslist$Q [8,4]=20
idparslist$Q [7,3]=20

#CD
idparslist$Q [6,8]=21
idparslist$Q [5,7]=21
idparslist$Q [8,6]=22
idparslist$Q [7,5]=22

p=starting.point.musse(mamphylo,4)[c(1:20,20,20)]
fitsecsse=secsse_ml(phy = mamphylo,tipstates, num_concealed_states = num_concealed_states,
                        cond = "proper_cond", root_state_weight = "proper_weights",
                        idparslist = idparslist,
                        idparsopt = c(1:22),
                        initparsopt= p,
                        idparsfix = 0,
                        parsfix = 0,
                        sampling_fraction =c(1,1),
                        optimmethod = "subplex",
                        tol = c(1e-05, 1e-06, 1e-08))

paste("all together model finished")

#by order
table(bakedranges$order)
orders=c("CARNIVORA", "CETARTIODACTYLA","CHIROPTERA","EULIPOTYPHLA", "PRIMATES", "RODENTIA")
geosseislandsbyorders=function(order){
  bakedrangesord=bakedranges[which(bakedranges$order==order),]
  mamphyloord=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesord$species))
  
  tipstatesord=bakedrangesord$abovemedian+1
  
  num_concealed_statesord=4
  idparslistord=id_paramPos(tipstatesord,num_concealed_statesord)
  
  idparslistord$lambdas[1]=1
  idparslistord$lambdas[2]=1
  idparslistord$lambdas[3]=2
  idparslistord$lambdas[4]=2
  idparslistord$lambdas[5]=3
  idparslistord$lambdas[6]=3
  idparslistord$lambdas[7]=4
  idparslistord$lambdas[8]=4
  
  idparslistord$mus[1]=5
  idparslistord$mus[2]=5
  idparslistord$mus[3]=6
  idparslistord$mus[4]=6
  idparslistord$mus[5]=7
  idparslistord$mus[6]=7
  idparslistord$mus[7]=8
  idparslistord$mus[8]=8
  
  idparslistord$Q=idparslistord$Q
  idparslistord$Q[,]=0
  diag(idparslistord$Q[,])=NA
  
  idparslistord$Q [1,2]=9
  idparslistord$Q [2,1]=10
  idparslistord$Q [3,4]=9
  idparslistord$Q [4,3]=10
  idparslistord$Q [5,6]=9
  idparslistord$Q [6,5]=10
  idparslistord$Q [7,8]=9
  idparslistord$Q [8,7]=10
  
  #AB
  idparslistord$Q [2,4]=11
  idparslistord$Q [1,3]=11
  idparslistord$Q [4,2]=12
  idparslistord$Q [3,1]=12
  
  #AC
  idparslistord$Q [2,6]=13
  idparslistord$Q [1,5]=13
  idparslistord$Q [6,2]=14
  idparslistord$Q [5,1]=14
  
  #AD
  idparslistord$Q [2,8]=15
  idparslistord$Q [1,7]=15
  idparslistord$Q [8,2]=16
  idparslistord$Q [7,1]=16
  
  #BC
  idparslistord$Q [4,6]=17
  idparslistord$Q [3,5]=17
  idparslistord$Q [6,4]=18
  idparslistord$Q [5,3]=18
  
  #BD
  idparslistord$Q [4,8]=19
  idparslistord$Q [3,7]=19
  idparslistord$Q [8,4]=20
  idparslistord$Q [7,3]=20
  
  #CD
  idparslistord$Q [6,8]=21
  idparslistord$Q [5,7]=21
  idparslistord$Q [8,6]=22
  idparslistord$Q [7,5]=22
  
  pord=starting.point.musse(mamphylo,4)[c(1:20,20,20)]
  fitsecsseord=secsse_ml(phy = mamphyloord,tipstatesord, num_concealed_states = num_concealed_statesord,
                          cond = "proper_cond", root_state_weight = "proper_weights",
                          idparslist = idparslistord,
                          idparsopt = c(1:22),
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


save.image("medianSecSSEmodeliv.RData")

