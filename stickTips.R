stickTips <- 
  function(tree, tab , prun=T, subtree.type = "polytomy", b=0.01, d=0, p=0.3, comments=TRUE){
    
    
    if(!require(ape)) {stop("library ape is needed!\n")}
    if(!require(adephylo)) {stop("library adephylo is needed!\n")}
    library(ape)
    library(adephylo)
    
    tab[,1] <- as.character(tab[,1])
    tab[,2] <- as.character(tab[,2])
    
    nochanges <- which(tab[,1]==tab[,2])
    if( sum(duplicated(c(unique(tab[-nochanges,1]), unique(tab[-nochanges,2])))) != 0 ) {stop("two tips have the same name!\n")}   # this will avoid the creation of a weird object
    
    #select table names that are not in the tips
    outers.tab <- unique(tab[!(tab[,2] %in% tree$tip.label),2])
    #warning for these names and removing of the matching lines
    if (length(outers.tab)!=0) { warning(paste(outers.tab,"is/are not in the tips of your tree !\n"),call.=FALSE)
      tab <- tab[- which((tab[,2] %in% outers.tab)),]
    }
    
    #select tips names that are not in the table
    outers.tree <- unique(tree$tip.label[!(tree$tip.label %in% tab[,2])])
    #warning if the tree contains more names than the table
    if (length(outers.tree)!=0) { if(comments==TRUE) warning(paste(outers.tree,"is/are not in your tab !\n"),call.=FALSE)    }
    
    # sticking
    for(tip in as.character(unique(tab[,2]))){ #each tip
      
      toStick <- tab[tab[,2]==tip,1]  #new tips list
      if(comments==TRUE) print(paste("sticking to : ",as.character(tip)))
      #if there is no element, we do nothing (case only if prun=F)
      #substitution if only 1 element
      if(length(toStick)==1) tree$tip.label[tree$tip.label==tip] <- as.character(toStick[1])
      #if several elements
      if(length(toStick)>1) {
        # depth of the branch:
        depth <- tree$edge.length[which.edge(tree,which(tree$tip.label==tip))]
        
        if(subtree.type == "polytomy"){
          #shorten the receiver branch
          tree$edge.length[which.edge(tree,which(tree$tip.label==tip))] = depth/2
          #make the polytomic tree to stick
          polytomy <-  as.phylo(~newTip, data=data.frame(newTip=toStick))
          subtree = compute.brlen(polytomy,depth/2)
        }else {
          if(subtree.type=="brownian"){
            if(!require(geiger)) {stop("library geiger is needed!\n")}
            library(geiger)
            # subtree
            subtree <- treedata(birthdeath.tree(b, d, taxa.stop=(length(toStick) + 1)), data.frame(rnorm(length(toStick))), warnings=FALSE)$phy
            
            # create tmp to get root length
            tmp <- treedata(birthdeath.tree(b, d, taxa.stop=4), data.frame(rnorm(3)), warnings=FALSE)$phy
            root.tmp <- which(tmp$edge[,2] == 5)
            # shorten the receiver branch
            # old: receiv.brlen <- mean(subtree$edge.length) * depth/max(distRoot(subtree))
            receiv.brlen <- tmp$edge.length[root.tmp] * depth/(max(distRoot(subtree)) + tmp$edge.length[root.tmp])
            
          }else{
            if(!require(apTreeshape)) {stop("library apTreeshape is needed!\n")}
            library(apTreeshape)
            #make the tree to stick
            subtree <- as.phylo(rtreeshape(1, tip.number = length(toStick), model=subtree.type, p=p)[[1]])
            
            #shorten the receiver branch
            receiv.brlen <- depth/max(distRoot(subtree))
            
          }
          #shorten the receiver branch
          tree$edge.length[which.edge(tree,which(tree$tip.label==tip))]  <- receiv.brlen
          #rescale new branch length
          subtree$edge.length <- subtree$edge.length*(depth-receiv.brlen)/max(distRoot(subtree))
          # tips names
          subtree$tip.label <- toStick
        }
        #stick the subtree
        tree <- bind.tree(tree,subtree, which(tree$tip.label==tip))
      }
    }
    #remove tips that are not in the table
    if (prun == TRUE & length(outers.tree)!=0) tree <- drop.tip(tree,outers.tree)
    return(tree)
  }
