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


#detele tips with no range info, island endemics and chiroptera
tipstodelete=setdiff(mamphylo$tip.label,bakedranges$species)
mamphylo=drop.tip(mamphylo,tipstodelete)

#categorize species by treshold
bakedranges$above250000=as.numeric(250000<bakedranges$area)

#make secsse all
tipstates=bakedranges$above250000+1

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


#by order
table(bakedranges$order)
orders=c("CARNIVORA", "CETARTIODACTYLA", "CHIROPTERA","EULIPOTYPHLA", "PRIMATES", "RODENTIA")
geosseislandsbyorders=function(order){
  bakedrangesord=bakedranges[which(bakedranges$order==order),]
  mamphyloord=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesord$species))
  
  tipstatesord=bakedrangesord$above250000+1
  
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


save.image("250SecSSEmodeliii.RData")

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

save.image("250SecSSEmodeiii.RData")


# load("medianSecSSEaf.RData")
#
# #ASR
# #carnivora
# bakedrangescarn=bakedranges[which(bakedranges$order=="CARNIVORA"),]
# mamphylocarn=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangescarn$species))
# tipstatescarn=bakedrangescarn$abovemedian+1
# 
# asrcarn=cla_secsse_loglik(orderresult[[1]]$MLpars, 
#                           phy=mamphylocarn, tipstatescarn, num_concealed_states=num_concealed_states,
#                             use_fortran = T, methode = "ode45", cond = "proper_cond",
#                             root_state_weight = "proper_weights", sampling_fraction=c(1,1),
#                             run_parallel = FALSE, setting_calculation = NULL,
#                             setting_parallel = NULL, see_ancestral_states = T,
#                             loglik_penalty = 0)
# 
# ls=cbind(asrcarn$ancestral_states[,1]+asrcarn$ancestral_states[,3],asrcarn$ancestral_states[,2]+asrcarn$ancestral_states[,4])
# ab=cbind(asrcarn$ancestral_states[,1]+asrcarn$ancestral_states[,2],asrcarn$ancestral_states[,3]+asrcarn$ancestral_states[,4])
# 
# plot(mamphylocarn)
# tiplabels(col=tipstatescarn+1,pch = 20)
# #large and small
# nodelabels(pie=ls,pch = 20,cex=0.5,piecol=c("red","green"))
# #first and second
# nodelabels(pie=ab,pch = 20,cex=0.5,piecol=c("blue","orange"))
# 
# ab[which(ab[,2]>0.7),]
# ex=which(ab[,2]>0.7)+length(mamphylocarn$tip.label)
# 
# carnexamples=lapply(ex, extract.clade, phy=mamphylocarn)
# 
# plot(carnexamples[[1]])
# plot(carnexamples[[2]])
# 
# 
# #ungulates
# bakedrangesceta=bakedranges[which(bakedranges$order=="CETARTIODACTYLA"),]
# mamphyloceta=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesceta$species))
# tipstatesceta=bakedrangesceta$abovemedian+1
# 
# asrceta=cla_secsse_loglik(orderresult[[2]]$MLpars, 
#                           phy=mamphyloceta, tipstatesceta, num_concealed_states=num_concealed_states,
#                           use_fortran = T, methode = "ode45", cond = "proper_cond",
#                           root_state_weight = "proper_weights", sampling_fraction=c(1,1),
#                           run_parallel = FALSE, setting_calculation = NULL,
#                           setting_parallel = NULL, see_ancestral_states = T,
#                           loglik_penalty = 0)
# 
# ls=cbind(asrceta$ancestral_states[,1]+asrceta$ancestral_states[,3],asrceta$ancestral_states[,2]+asrceta$ancestral_states[,4])
# ab=cbind(asrceta$ancestral_states[,1]+asrceta$ancestral_states[,2],asrceta$ancestral_states[,3]+asrceta$ancestral_states[,4])
# 
# plot(mamphyloceta)
# tiplabels(col=tipstatesceta+1,pch = 20)
# #large and small
# nodelabels(pie=ls,pch = 20,cex=0.5,piecol=c("red","green"))
# #first and second
# nodelabels(pie=ab,pch = 20,cex=0.5,piecol=c("blue","orange"))
# 
# ab[which(ab[,2]>0.7),]
# ex=which(ab[,2]>0.7)+length(mamphyloceta$tip.label)
# 
# cetaexamples=lapply(ex, extract.clade, phy=mamphyloceta)
# 
# plot(cetaexamples[[1]])
# plot(cetaexamples[[2]])
# plot(cetaexamples[[3]])
# plot(cetaexamples[[4]])
# plot(cetaexamples[[5]])
# 
# #eulipotyphla
# bakedrangeseuli=bakedranges[which(bakedranges$order=="EULIPOTYPHLA"),]
# mamphyloeuli=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangeseuli$species))
# tipstateseuli=bakedrangeseuli$abovemedian+1
# 
# asreuli=cla_secsse_loglik(orderresult[[3]]$MLpars, 
#                           phy=mamphyloeuli, tipstateseuli, num_concealed_states=num_concealed_states,
#                           use_fortran = T, methode = "ode45", cond = "proper_cond",
#                           root_state_weight = "proper_weights", sampling_fraction=c(1,1),
#                           run_parallel = FALSE, setting_calculation = NULL,
#                           setting_parallel = NULL, see_ancestral_states = T,
#                           loglik_penalty = 0)
# 
# ls=cbind(asreuli$ancestral_states[,1]+asreuli$ancestral_states[,3],asreuli$ancestral_states[,2]+asreuli$ancestral_states[,4])
# ab=cbind(asreuli$ancestral_states[,1]+asreuli$ancestral_states[,2],asreuli$ancestral_states[,3]+asreuli$ancestral_states[,4])
# 
# plot(mamphyloeuli)
# tiplabels(col=tipstateseuli+1,pch = 20)
# #large and small
# nodelabels(pie=ls,pch = 20,cex=0.5,piecol=c("red","green"))
# #first and second
# nodelabels(pie=ab,pch = 20,cex=0.5,piecol=c("blue","orange"))
# 
# ab[which(ab[,2]>0.7),]
# ex=which(ab[,2]>0.7)+length(mamphyloeuli$tip.label)
# 
# euliexamples=lapply(ex, extract.clade, phy=mamphyloeuli)
# 
# plot(euliexamples[[1]])
# plot(euliexamples[[2]])
# plot(euliexamples[[3]])
# plot(euliexamples[[4]])
# plot(euliexamples[[5]])
# plot(euliexamples[[6]])
# plot(euliexamples[[7]])
# plot(euliexamples[[8]])
# plot(euliexamples[[9]])
# plot(euliexamples[[10]])
# plot(euliexamples[[11]])
# plot(euliexamples[[12]])
# plot(euliexamples[[13]])
# plot(euliexamples[[14]])
# plot(euliexamples[[15]])
# plot(euliexamples[[16]])
# plot(euliexamples[[17]])
# plot(euliexamples[[18]])
# plot(euliexamples[[19]])
# plot(euliexamples[[20]])
# plot(euliexamples[[21]])
# 
# #primates
# bakedrangesprim=bakedranges[which(bakedranges$order=="PRIMATES"),]
# mamphyloprim=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesprim$species))
# tipstatesprim=bakedrangesprim$abovemedian+1
# 
# asrprim=cla_secsse_loglik(orderresult[[4]]$MLpars, 
#                           phy=mamphyloprim, tipstatesprim, num_concealed_states=num_concealed_states,
#                           use_fortran = T, methode = "ode45", cond = "proper_cond",
#                           root_state_weight = "proper_weights", sampling_fraction=c(1,1),
#                           run_parallel = FALSE, setting_calculation = NULL,
#                           setting_parallel = NULL, see_ancestral_states = T,
#                           loglik_penalty = 0)
# 
# ls=cbind(asrprim$ancestral_states[,1]+asrprim$ancestral_states[,3],asrprim$ancestral_states[,2]+asrprim$ancestral_states[,4])
# ab=cbind(asrprim$ancestral_states[,1]+asrprim$ancestral_states[,2],asrprim$ancestral_states[,3]+asrprim$ancestral_states[,4])
# 
# plot(mamphyloprim)
# tiplabels(col=tipstatesprim+1,pch = 20)
# #large and small
# nodelabels(pie=ls,pch = 20,cex=0.5,piecol=c("red","green"))
# #first and second
# nodelabels(pie=ab,pch = 20,cex=0.5,piecol=c("blue","orange"))
# 
# ab[which(ab[,2]>0.9),]
# ex=which(ab[,2]>0.9)+length(mamphyloprim$tip.label)
# 
# primexamples=lapply(ex, extract.clade, phy=mamphyloprim)
# 
# plot(primexamples[[1]])
# plot(primexamples[[2]])
# plot(primexamples[[3]])
# plot(primexamples[[4]])
# plot(primexamples[[5]])
# plot(primexamples[[6]])
# plot(primexamples[[7]])
# plot(primexamples[[8]])
# plot(primexamples[[9]])
# plot(primexamples[[10]])
# plot(primexamples[[11]])
# plot(primexamples[[12]])
# plot(primexamples[[13]])
# plot(primexamples[[14]])
# plot(primexamples[[15]])
# plot(primexamples[[16]])
# plot(primexamples[[17]])
# plot(primexamples[[18]])
# plot(primexamples[[19]])
# plot(primexamples[[20]])
# plot(primexamples[[21]])
# plot(primexamples[[22]])
# plot(primexamples[[23]])
# plot(primexamples[[24]])
# 
# #rodents
# bakedrangesrode=bakedranges[which(bakedranges$order=="RODENTIA"),]
# mamphylorode=drop.tip(mamphylo,setdiff(mamphylo$tip.label,bakedrangesrode$species))
# tipstatesrode=bakedrangesrode$abovemedian+1
# 
# asrrode=cla_secsse_loglik(orderresult[[5]]$MLpars, 
#                           phy=mamphylorode, tipstatesrode, num_concealed_states=num_concealed_states,
#                           use_fortran = T, methode = "ode45", cond = "proper_cond",
#                           root_state_weight = "proper_weights", sampling_fraction=c(1,1),
#                           run_parallel = FALSE, setting_calculation = NULL,
#                           setting_parallel = NULL, see_ancestral_states = T,
#                           loglik_penalty = 0)
# 
# ls=cbind(asrrode$ancestral_states[,1]+asrrode$ancestral_states[,3],asrrode$ancestral_states[,2]+asrrode$ancestral_states[,4])
# ab=cbind(asrrode$ancestral_states[,1]+asrrode$ancestral_states[,2],asrrode$ancestral_states[,3]+asrrode$ancestral_states[,4])
# 
# plot(mamphylorode)
# tiplabels(col=tipstatesrode+1,pch = 20)
# #large and small
# nodelabels(pie=ls,pch = 20,cex=0.5,piecol=c("red","green"))
# #first and second
# nodelabels(pie=ab,pch = 20,cex=0.5,piecol=c("blue","orange"))
# 
# ab[which(ab[,2]>0.47),]
# ex=which(ab[,2]>0.47)+length(mamphylorode$tip.label)
# 
# rodeexamples=lapply(ex, extract.clade, phy=mamphylorode)
# 
# plot(rodeexamples[[1]])
# plot(rodeexamples[[2]])
# plot(rodeexamples[[3]])
# plot(rodeexamples[[4]])
# 
# ls[which(ls[,1]>0.90),]
# ex=which(ls[,1]>0.90)+length(mamphylorode$tip.label)
# 
# rodeexamples=lapply(ex, extract.clade, phy=mamphylorode)
# 
# plot(rodeexamples[[1]])
# plot(rodeexamples[[2]])
# plot(rodeexamples[[3]])
# plot(rodeexamples[[4]])
# 
# # orderresultdf=data.frame(matrix(unlist(orderresult), nrow=length(orderresult), byrow=T))
# # orderresultdf=orderresultdf[,c(1,5:10)]
# # names(orderresultdf)=c("speciation small", "speciation clado", "speciation large", 
# #                "extinction small","extinction large", "expansion", "shrinkage")
# # 
# # row.names(orderresultdf)=orders
# # 
# # plot(1:5,orderresultdf$`speciation small`-orderresultdf$`extinction small`,ylab="net diversification of small ranged species",xaxt="n",pch=19,xlab="",xlim=c(0,7), ylim=c(-1,1))
# # abline(0,0, lty=2)
# # gplots::angleAxis(1,at=1:5,labels=orders)
# # 
# # plot(1:5,orderresultdf$`speciation large`-orderresultdf$`extinction large`,ylab="net diversification of large ranged species",xaxt="n",pch=19,xlab="",xlim=c(0,7), ylim=c(-1,1))
# # abline(0,0, lty=2)
# # gplots::angleAxis(1,at=1:5,labels=orders)
# # 
# # plot(1:5,orderresultdf$`speciation clado`,ylab="state change speciation",xaxt="n",pch=19,xlab="",xlim=c(0,7), ylim=c(-1,1))
# # abline(0,0, lty=2)
# # gplots::angleAxis(1,at=1:5,labels=orders)
# # 
# # 
# # 
# # orderresultdf$repvallarge=0
# # 
# # 
# # for (i in 1:length(orders)){
# #   fulp=c(orderresultdf$`speciation small`[i],0,0,0,orderresultdf$`speciation clado`[i],orderresultdf$`speciation large`[i],
# #          orderresultdf$`extinction small`[i], orderresultdf$`extinction large`[i],
# #          orderresultdf$expansion[i],orderresultdf$shrinkage[i])
# #   tmat=diversitree:::projection.matrix.classe(fulp,2)
# #   orderresultdf$repvallarge[i]=reproductive.value(tmat)[2]
# # }
# # 
# # plot(1:5,log(orderresultdf$repvallarge),ylab="log(reproductive value of large ranged species)",xaxt="n",pch=19,xlab="",xlim=c(0,7), ylim=c(-5,20))
# # abline(0,0, lty=2)
# # gplots::angleAxis(1,at=1:5,labels=orders)
# # 
# # write.csv(orderresultdf,file="medianSSEresults.csv")
# # 
# # 
# # 
# # 
