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

load("250SecSSEmodeliii.RData")
source("./cla_secsse_loglikab.R")

cores=2

#carnivora
bakedrangescarn=bakedranges[which(bakedranges$order=="CARNIVORA"),]
mamphylocarn=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangescarn$species))
tipstatescarn=bakedrangescarn$above250000+1

#cetartiodactyla
bakedrangesceta=bakedranges[which(bakedranges$order=="CETARTIODACTYLA"),]
mamphyloceta=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesceta$species))
tipstatesceta=bakedrangesceta$above250000+1

#chiroptera
bakedrangeschiro=bakedranges[which(bakedranges$order=="CHIROPTERA"),]
mamphylochiro=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeschiro$species))
tipstateschiro=bakedrangeschiro$above250000+1

#eulipotyphla
bakedrangeseuli=bakedranges[which(bakedranges$order=="EULIPOTYPHLA"),]
mamphyloeuli=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeseuli$species))
tipstateseuli=bakedrangeseuli$above250000+1

#primates
bakedrangesprim=bakedranges[which(bakedranges$order=="PRIMATES"),]
mamphyloprim=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesprim$species))
tipstatesprim=bakedrangesprim$above250000+1

#rodentia
bakedrangesrode=bakedranges[which(bakedranges$order=="RODENTIA"),]
mamphylorode=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesrode$species))
tipstatesrode=bakedrangesrode$above250000+1


#################
#calculations
tipestimatescarn=cla_secsse_propab_all(orderresult[[1]]$MLpars,phy=mamphylocarn, tipstatescarn, mc.cores = cores)
paste("carnivora finished")
save(tipestimatescarn, file="tipest250SecSSEcarn.RData")

tipestimatesceta=cla_secsse_propab_all(orderresult[[2]]$MLpars,phy=mamphyloceta, tipstatesceta, mc.cores = cores)
paste("cetartiodactyla finished")
save(tipestimatesceta, file="tipest250SecSSEceta.RData")

tipestimateschiro=cla_secsse_propab_all(orderresult[[3]]$MLpars,phy=mamphylochiro, tipstateschiro, mc.cores = cores)
paste("chiroptera finished")
save(tipestimateschiro, file="tipest250SecSSEchiro.RData")

tipstimateseuli=cla_secsse_propab_all(orderresult[[4]]$MLpars,phy=mamphyloeuli, tipstateseuli, mc.cores = cores)
paste("eulipotyphla finished")
save(tipstimateseuli, file="tipest250SecSSEeuli.RData")

tipestimatesprim=cla_secsse_propab_all(orderresult[[5]]$MLpars,phy=mamphyloprim, tipstatesprim, mc.cores = cores, reltol=2e-12)
paste("primates finished")
save(tipestimatesprim, file="tipest250SecSSEprim.RData")

tipestimatesrode=cla_secsse_propab_all(orderresult[[6]]$MLpars,phy=mamphylorode, tipstatesrode, mc.cores = cores)
paste("rodentia finished")
save(tipestimatesrode, file="tipest250SecSSErode.RData")

tipestimatesall=cla_secsse_propab_all(fitsecsse$MLpars,phy=mamphylo, tipstates, mc.cores = cores)
paste("alltogether finished")
save(tipestimatesall, file="tipest250SecSSEall.RData")

