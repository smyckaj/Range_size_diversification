#!/usr/bin/Rscript
setwd("~/Desktop/honza/Diversification_and_range_size")
#setwd("~/Dropbox/CTS postdoc/Diversification_and_range_size/")

library(ape)
library(phytools)
library(phangorn)
library(diversitree)
library(secsse)
library(parallel)
library(Rmpfr)

load("medianSecSSEmodeliii.RData")
source("./cla_secsse_loglikab.R")

cores=30

#carnivora
bakedrangescarn=bakedranges[which(bakedranges$order=="CARNIVORA"),]
mamphylocarn=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangescarn$species))
tipstatescarn=bakedrangescarn$abovemedian+1

#cetartiodactyla
bakedrangesceta=bakedranges[which(bakedranges$order=="CETARTIODACTYLA"),]
mamphyloceta=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesceta$species))
tipstatesceta=bakedrangesceta$abovemedian+1

#chiroptera
bakedrangeschiro=bakedranges[which(bakedranges$order=="CHIROPTERA"),]
mamphylochiro=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeschiro$species))
tipstateschiro=bakedrangeschiro$abovemedian+1

#eulipotyphla
bakedrangeseuli=bakedranges[which(bakedranges$order=="EULIPOTYPHLA"),]
mamphyloeuli=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeseuli$species))
tipstateseuli=bakedrangeseuli$abovemedian+1

#primates
bakedrangesprim=bakedranges[which(bakedranges$order=="PRIMATES"),]
mamphyloprim=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesprim$species))
tipstatesprim=bakedrangesprim$abovemedian+1

#rodentia
bakedrangesrode=bakedranges[which(bakedranges$order=="RODENTIA"),]
mamphylorode=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesrode$species))
tipstatesrode=bakedrangesrode$abovemedian+1


#################
#calculations
tipestimatescarn=cla_secsse_propab_all(orderresult[[1]]$MLpars,phy=mamphylocarn, tipstatescarn, mc.cores = cores)
paste("carnivora finished")
save(tipestimatescarn, file="tipestmedianSecSSEcarn.RData")

tipestimatesceta=cla_secsse_propab_all(orderresult[[2]]$MLpars,phy=mamphyloceta, tipstatesceta, mc.cores = cores)
paste("cetartiodactyla finished")
save(tipestimatesceta, file="tipestmedianSecSSEceta.RData")

tipestimateschiro=cla_secsse_propab_all(orderresult[[3]]$MLpars,phy=mamphylochiro, tipstateschiro, mc.cores = cores)
paste("chiroptera finished")
save(tipestimateschiro, file="tipestmedianSecSSEchiro.RData")

tipstimateseuli=cla_secsse_propab_all(orderresult[[4]]$MLpars,phy=mamphyloeuli, tipstateseuli, mc.cores = cores)
paste("eulipotyphla finished")
save(tipstimateseuli, file="tipestmedianSecSSEeuli.RData")

tipestimatesprim=cla_secsse_propab_all(orderresult[[5]]$MLpars,phy=mamphyloprim, tipstatesprim, mc.cores = cores, reltol=5e-12)
paste("primates finished")
save(tipestimatesprim, file="tipestmedianSecSSEprim.RData")

tipestimatesrode=cla_secsse_propab_all(orderresult[[6]]$MLpars,phy=mamphylorode, tipstatesrode, mc.cores = cores)
paste("rodentia finished")
save(tipestimatesrode, file="tipestmedianSecSSErode.RData")

tipestimatesall=cla_secsse_propab_all(fitsecsse$MLpars,phy=mamphylo, tipstates, mc.cores = cores)
paste("alltogether finished")
save(tipestimatesall, file="tipestmedianSecSSEall.RData")

