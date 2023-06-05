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

load("varmedianSecSSEmodeliii.RData")
load("varmedianSecSSEmodeliiirode.RData")
source("./cla_secsse_loglikab.R")

cores=12


#rodentia
bakedrangesrode=bakedranges[which(bakedranges$order=="RODENTIA"),]
mamphylorode=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesrode$species))
tipstatesrode=as.numeric(median(bakedrangesrode$area)<bakedrangesrode$area)+1


#################
#calculations

tipestimatesrode=cla_secsse_propab_all(roderesult[[1]]$MLpars,phy=mamphylorode, tipstatesrode, mc.cores = cores)
paste("rodentia finished")
save(tipestimatesrode, file="tipestvarmedianSecSSErode.RData")

