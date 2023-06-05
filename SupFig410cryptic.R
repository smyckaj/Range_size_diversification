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
library(scales)
library(vioplot)
library(rphylopic)

#load secsse estimates
load("cryptic10medianSecSSEmodeliii.RData")

#load tip state predictions
load("tipestcryptic10medianSecSSEall.RData")
load("tipestcryptic10medianSecSSEcarn.RData")
load("tipestcryptic10medianSecSSEceta.RData")
load("tipestcryptic10medianSecSSEchiro.RData")
load("tipestcryptic10medianSecSSEeuli.RData")
load("tipestcryptic10medianSecSSEprim.RData")
load("tipestcryptic10medianSecSSErode.RData")


#get the tip estimates
##########################################

ndall=c(fitsecsse$MLpars[[1]][[1]][1,1]-fitsecsse$MLpars[[2]][1],
        fitsecsse$MLpars[[1]][[2]][2,1]+fitsecsse$MLpars[[1]][[2]][2,2]-fitsecsse$MLpars[[2]][2],
        fitsecsse$MLpars[[1]][[3]][3,3]-fitsecsse$MLpars[[2]][3],
        fitsecsse$MLpars[[1]][[4]][4,3]+fitsecsse$MLpars[[1]][[4]][4,4]-fitsecsse$MLpars[[2]][4])
        

tipndall=rowSums(t(t(tipestimatesall) * ndall))


ndcarn=c(orderresult$CARNIVORA$MLpars[[1]][[1]][1,1]-orderresult$CARNIVORA$MLpars[[2]][1],
         orderresult$CARNIVORA$MLpars[[1]][[2]][2,1]+orderresult$CARNIVORA$MLpars[[1]][[2]][2,2]-orderresult$CARNIVORA$MLpars[[2]][2],
         orderresult$CARNIVORA$MLpars[[1]][[3]][3,3]-orderresult$CARNIVORA$MLpars[[2]][3],
         orderresult$CARNIVORA$MLpars[[1]][[4]][4,3]+orderresult$CARNIVORA$MLpars[[1]][[4]][4,4]-orderresult$CARNIVORA$MLpars[[2]][4])


tipndcarn=rowSums(t(t(tipestimatescarn) * ndcarn))

ndceta=c(orderresult$CETARTIODACTYLA$MLpars[[1]][[1]][1,1]-orderresult$CETARTIODACTYLA$MLpars[[2]][1],
         orderresult$CETARTIODACTYLA$MLpars[[1]][[2]][2,1]+orderresult$CETARTIODACTYLA$MLpars[[1]][[2]][2,2]-orderresult$CETARTIODACTYLA$MLpars[[2]][2],
         orderresult$CETARTIODACTYLA$MLpars[[1]][[3]][3,3]-orderresult$CETARTIODACTYLA$MLpars[[2]][3],
         orderresult$CETARTIODACTYLA$MLpars[[1]][[4]][4,3]+orderresult$CETARTIODACTYLA$MLpars[[1]][[4]][4,4]-orderresult$CETARTIODACTYLA$MLpars[[2]][4])


tipndceta=rowSums(t(t(tipestimatesceta) * ndceta))

ndchiro=c(orderresult$CHIROPTERA$MLpars[[1]][[1]][1,1]-orderresult$CHIROPTERA$MLpars[[2]][1],
         orderresult$CHIROPTERA$MLpars[[1]][[2]][2,1]+orderresult$CHIROPTERA$MLpars[[1]][[2]][2,2]-orderresult$CHIROPTERA$MLpars[[2]][2],
         orderresult$CHIROPTERA$MLpars[[1]][[3]][3,3]-orderresult$CHIROPTERA$MLpars[[2]][3],
         orderresult$CHIROPTERA$MLpars[[1]][[4]][4,3]+orderresult$CHIROPTERA$MLpars[[1]][[4]][4,4]-orderresult$CHIROPTERA$MLpars[[2]][4])


tipndchiro=rowSums(t(t(tipestimateschiro) * ndchiro))

ndeuli=c(orderresult$EULIPOTYPHLA$MLpars[[1]][[1]][1,1]-orderresult$EULIPOTYPHLA$MLpars[[2]][1],
          orderresult$EULIPOTYPHLA$MLpars[[1]][[2]][2,1]+orderresult$EULIPOTYPHLA$MLpars[[1]][[2]][2,2]-orderresult$EULIPOTYPHLA$MLpars[[2]][2],
          orderresult$EULIPOTYPHLA$MLpars[[1]][[3]][3,3]-orderresult$EULIPOTYPHLA$MLpars[[2]][3],
          orderresult$EULIPOTYPHLA$MLpars[[1]][[4]][4,3]+orderresult$EULIPOTYPHLA$MLpars[[1]][[4]][4,4]-orderresult$EULIPOTYPHLA$MLpars[[2]][4])


tipndeuli=rowSums(t(t(tipstimateseuli) * ndeuli))

ndprim=c(orderresult$PRIMATES$MLpars[[1]][[1]][1,1]-orderresult$PRIMATES$MLpars[[2]][1],
         orderresult$PRIMATES$MLpars[[1]][[2]][2,1]+orderresult$PRIMATES$MLpars[[1]][[2]][2,2]-orderresult$PRIMATES$MLpars[[2]][2],
         orderresult$PRIMATES$MLpars[[1]][[3]][3,3]-orderresult$PRIMATES$MLpars[[2]][3],
         orderresult$PRIMATES$MLpars[[1]][[4]][4,3]+orderresult$PRIMATES$MLpars[[1]][[4]][4,4]-orderresult$PRIMATES$MLpars[[2]][4])


tipndprim=rowSums(t(t(tipestimatesprim) * ndprim))

ndrode=c(orderresult$RODENTIA$MLpars[[1]][[1]][1,1]-orderresult$RODENTIA$MLpars[[2]][1],
         orderresult$RODENTIA$MLpars[[1]][[2]][2,1]+orderresult$RODENTIA$MLpars[[1]][[2]][2,2]-orderresult$RODENTIA$MLpars[[2]][2],
         orderresult$RODENTIA$MLpars[[1]][[3]][3,3]-orderresult$RODENTIA$MLpars[[2]][3],
         orderresult$RODENTIA$MLpars[[1]][[4]][4,3]+orderresult$RODENTIA$MLpars[[1]][[4]][4,4]-orderresult$RODENTIA$MLpars[[2]][4])


tipndrode=rowSums(t(t(tipestimatesrode) * ndrode))

#get the transition estimates
##########################################
#construction of transition matrices
transmatsecsse=function(fitsecsse){
  alltransmat=matrix(0,4,4)
  
  #1A
  alltransmat[1,1]=fitsecsse$MLpars[[1]][[1]][1,1]+ #sympatric
    -fitsecsse$MLpars[[2]][1]+ #extinction
    -fitsecsse$MLpars[[3]][1,2]+ #emigration to large
    -fitsecsse$MLpars[[3]][1,3] #emigration to rate shift
  
  alltransmat[1,2]=fitsecsse$MLpars[[1]][[2]][2,1]+ #range shift speciation
    fitsecsse$MLpars[[3]][2,1] #imigration from large
  
  alltransmat[1,3]=fitsecsse$MLpars[[3]][3,1] #imigration from rate shift
  
  #2A
  alltransmat[2,1]=fitsecsse$MLpars[[3]][1,2] #imigration from small
  
  alltransmat[2,2]=fitsecsse$MLpars[[1]][[2]][2,2]+ #sympatric
    -fitsecsse$MLpars[[2]][2]+ #extinction
    -fitsecsse$MLpars[[3]][2,1]+ #emigration to small
    -fitsecsse$MLpars[[3]][2,4] #emigration to rate shift
  
  alltransmat[2,4]=fitsecsse$MLpars[[3]][4,2] #imigration from rate shift
  
  #1B
  alltransmat[3,1]=fitsecsse$MLpars[[3]][1,3] #imigration from rate shift
  
  alltransmat[3,3]=fitsecsse$MLpars[[1]][[3]][3,3]+ #sympatric
    -fitsecsse$MLpars[[2]][3]+ #extinction
    -fitsecsse$MLpars[[3]][3,4]+ #emigration to large
    -fitsecsse$MLpars[[3]][3,1] #emigration to rate shift
  
  alltransmat[3,4]=fitsecsse$MLpars[[1]][[4]][4,3]+ #range shift speciation
    fitsecsse$MLpars[[3]][4,3] #imigration from large
  
  #2B
  alltransmat[4,2]=fitsecsse$MLpars[[3]][2,4] #imigration from rate shift
  
  alltransmat[4,3]=fitsecsse$MLpars[[3]][3,4] #imigration from small
  
  alltransmat[4,4]=fitsecsse$MLpars[[1]][[4]][4,4]+ #sympatric
    -fitsecsse$MLpars[[2]][4]+ #extinction
    -fitsecsse$MLpars[[3]][4,3]+ #emigration to small
    -fitsecsse$MLpars[[3]][4,2] #emigration to rate shift
  return(alltransmat)
}

transmatorders=lapply(orderresult,transmatsecsse)
ssorders=lapply(transmatorders,stable.stage)

transmatall=transmatsecsse(fitsecsse)
ssall=stable.stage(transmatall)



#plots
##########################################
# par(mfrow=c(2,1))
# 
# par(mar=c(0,4,4,4))

#png("sf2.png",3000,1500,res=300)
svg("sfcryptic10.svg", 18/2.54, 9/2.54, pointsize=8)

plot(NULL, xlim=c(0,21),ylim=c(-1.2,1.2),xaxt="n",xlab="",ylab="net diversification rate")

vioplot(tipndall[tipestimatesall[,1]>0],at=1,col=2,add=T)
vioplot(tipndall[tipestimatesall[,2]>0],at=2,col=4,add=T)

vioplot(tipndcarn[tipestimatescarn[,1]>0],at=4,col=2,add=T)
vioplot(tipndcarn[tipestimatescarn[,2]>0],at=5,col=4,add=T)

vioplot(tipndceta[tipestimatesceta[,1]>0],at=7,col=2,add=T)
vioplot(tipndceta[tipestimatesceta[,2]>0],at=8,col=4,add=T)

vioplot(tipndchiro[tipestimateschiro[,1]>0],at=10,col=2,add=T)
vioplot(tipndchiro[tipestimateschiro[,2]>0],at=11,col=4,add=T)        

vioplot(tipndeuli[tipstimateseuli[,1]>0],at=13,col=2,add=T)
vioplot(tipndeuli[tipstimateseuli[,2]>0],at=14,col=4,add=T)

vioplot(tipndprim[tipestimatesprim[,1]>0],at=16,col=2,add=T)
vioplot(tipndprim[tipestimatesprim[,2]>0],at=17,col=4,add=T)

vioplot(tipndrode[tipestimatesrode[,1]>0],at=19,col=2,add=T)
vioplot(tipndrode[tipestimatesrode[,2]>0],at=20,col=4,add=T)

abline(0,0, lty=2)
abline(v=3)

carnsil <- image_data("e2015ba3-4f7e-4950-9bde-005e8678d77b", size = "1024")
cetasil <- image_data("dd52da0c-2d10-4f57-8013-6ffefbedcc40", size = "1024")
chirosil <- image_data("18bfd2fc-f184-4c3a-b511-796aafcc70f6", size = "1024")
eulisil <- image_data("823fa387-793a-41b0-8c91-73f4519504ce", size = "1024")
primsil <- image_data("d6cfb28f-136e-4a20-a5ac-8eb353c7fc4a", size = "1024")
rodesil <- image_data("be8670c2-a5bd-4b44-88e8-92f8b0c7f4c6", size = "1024")


add_phylopic_base2=function (img, x = NULL, y = NULL, ysize = NULL, alpha = 0.2, 
                             color = NULL,rxy=1) 
{
  img <- rphylopic:::recolor_phylopic(img, alpha, color)
  dims <- dim(img)[1:2]
  AR <- dims[1]/dims[2]
  xsize <- AR * ysize*rxy
  graphics::rasterImage(img, x - xsize/2, y - ysize/2, x + 
                          xsize/2, y + ysize/2, interpolate = F)
}

add_phylopic_base2(carnsil[[1]],4.5,1,0.3,alpha=1,rxy=12)
add_phylopic_base2(cetasil[[1]],7.5,1,0.3,alpha=1,rxy=3)
add_phylopic_base2(chirosil[[1]],10.5,1,0.25,alpha=1,rxy=10)
add_phylopic_base2(eulisil[[1]],13.5,1,0.2,alpha=1,rxy=15)
add_phylopic_base2(primsil[[1]],16.5,1,0.3,alpha=1,rxy=2)
add_phylopic_base2(rodesil[[1]],19.5,1,0.2,alpha=1,rxy=15)


# par(mar=c(6,4,0.5,4))
# 
# plot(NULL, xlim=c(0,21),ylim=c(-1,1),xaxt="n",xlab="",ylab="net diversification rate")
# points(rep(1.5,4),
#        ndall,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssall*10)
# points(rep(4.5,4),
#        ndcarn,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssorders$CARNIVORA*10)
# points(rep(7.5,4),
#        ndceta,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssorders$CETARTIODACTYLA*10)
# points(rep(10.5,4),
#        ndchiro,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssorders$CHIROPTERA*10)
# points(rep(13.5,4),
#        ndeuli,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssorders$EULIPOTYPHLA*10)
# points(rep(16.5,4),
#        ndprim,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssorders$PRIMATES*10)
# points(rep(19.5,4),
#        ndrode,
#        pch=20,col=alpha(c(2,3,2,3),0.7),cex=ssorders$RODENTIA*10)
# 
# abline(0,0, lty=2)

legend("bottomright",c("small-ranged","large-ranged"), col=c("red","blue"),pch=19,inset=0.0)
gplots::angleAxis(1,at=c(2,5,8,11,14,17,20),labels=c("all together", "Carnivora", "Cetartiodactyla", "Chiroptera", "Eulipotyphla","Primates", "Rodentia"))

text(1,0.6,paste(length(tipndall[tipestimatesall[,1]>0])),cex=0.8)
text(2,0.6,paste(length(tipndall[tipestimatesall[,2]>0])),cex=0.8)

text(4,0.6,paste(length(tipndcarn[tipestimatescarn[,1]>0])),cex=0.8)
text(5,0.6,paste(length(tipndcarn[tipestimatescarn[,2]>0])),cex=0.8)

text(7,0.6,paste(length(tipndceta[tipestimatesceta[,1]>0])),cex=0.8)
text(8,0.6,paste(length(tipndceta[tipestimatesceta[,2]>0])),cex=0.8)

text(10,0.6,paste(length(tipndchiro[tipestimateschiro[,1]>0])),cex=0.8)
text(11,0.6,paste(length(tipndchiro[tipestimateschiro[,2]>0])),cex=0.8)

text(13,0.6,paste(length(tipndeuli[tipstimateseuli[,1]>0])),cex=0.8)
text(14,0.6,paste(length(tipndeuli[tipstimateseuli[,2]>0])),cex=0.8)

text(16,0.6,paste(length(tipndprim[tipestimatesprim[,1]>0])),cex=0.8)
text(17,0.6,paste(length(tipndprim[tipestimatesprim[,2]>0])),cex=0.8)

text(19,0.6,paste(length(tipndrode[tipestimatesrode[,1]>0])),cex=0.8)
text(20,0.6,paste(length(tipndrode[tipestimatesrode[,2]>0])),cex=0.8)

dev.off()





