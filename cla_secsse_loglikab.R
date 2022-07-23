######################################
#This R script provides a set of functions that allow to determine marginal probabilities of species being in different concealed regimes of SecSSE diversification model. The current version of the script works for claSecSSE model with tho states, specifically the model iii as defined in Smycka et al. manuscript on range size~diversification dynamics in mammals.
#The workflow of the script is following: The functions cla_secsse_loglik a or b calculate a log-likelihood of SecSSE model parameters, given the phylogeny with mapped non-concealed states and known value of concealed state a or b for the focal species. In other words, these functions calculate a "log-probability" of focal species being in concealed state a or b and a phylogeny with mapped non-concealed states, given the SecSSe model parameters. The function cla_secsse_propab calculates a probability of focal species being in a or b, given the phylogeny with mapped non-concealed states, i.e. it exponentializes the cla_secsse_loglik a and b values, and marginalizes them. The cla_secsse_propab_all applies the cla_secsse_propab to all species in the given phylogeny. This function can be run in parallel on UNIX machines.
#The script also contains cla_calThruNodes2 and cla_secsse_loglik2 adopted from the SecSSE package, but allowing to change absolute tolerances of the ODE solver. This feature was added because for some of the species in Smycka et. al. mammals dataset, the likelihoods were below the limits of default ODE solver settings.
######################################


######################################
#Likelihood of focal species in concealed state a and phylogeny with tip traits, given the SecSSE parameters. The code is based on likelihood function in SecSSE, modified sections are highlighted with comments. The argument focalspecies is a string defining the species of interest from phy$tip.label, other arguments are the same as cla_secsse_loglik in SecSSE.

cla_secsse_loglika=function (parameter, phy, traits, focalspecies, num_concealed_states, use_fortran = TRUE, 
          methode = "ode45", cond = "proper_cond", root_state_weight = "proper_weights", 
          sampling_fraction, run_parallel = FALSE, setting_calculation = NULL, 
          setting_parallel = NULL, see_ancestral_states = FALSE, loglik_penalty = 0, reltol=1e-12,abstol=1e-16,hmax=NULL) 
{
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  if (run_parallel == TRUE) {
    if (is.null(setting_calculation)) {
      secsse:::check_input(traits, phy, sampling_fraction, root_state_weight)
      setting_calculation <- secsse:::build_initStates_time_bigtree(phy, 
                                                           traits, num_concealed_states, sampling_fraction)
    }
    
    #force the focal species to be in the concealed state a
    states <- setting_calculation$states
    focalspeciesorder=which(phy$tip.label==focalspecies)
    if(states[focalspeciesorder,8]==1 & states[focalspeciesorder,6]==1){
      states[focalspeciesorder,8]=0
      states[focalspeciesorder,6]=1
    }
    
    if(states[focalspeciesorder,7]==1 & states[focalspeciesorder,5]==1){
      states[focalspeciesorder,7]=0
      states[focalspeciesorder,5]=1
    }
    #end of modified code
    
    
    forTime <- setting_calculation$forTime
    ancesSub1 <- setting_calculation$ancesSub1
    ancesSub2 <- setting_calculation$ancesSub2
    ancesRest <- setting_calculation$ancesRest
    if (num_concealed_states != round(num_concealed_states)) {
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 
                                                      1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    take_ancesSub <- list(ancesSub1, ancesSub2)
    if (is.null(setting_parallel)) {
      cl <- parallel::makeCluster(2)
      doParallel::registerDoParallel(cl)
    }
    statesNEW <- cla_doParalThing(take_ancesSub, states, 
                                  loglik, forTime, parameter, use_fortran, methode, 
                                  phy)
    comingfromSub1 <- statesNEW[[1]][[1]]
    comingfromSub2 <- statesNEW[[2]][[1]]
    loglik <- statesNEW[[1]][[2]] + statesNEW[[2]][[2]]
    thoseCalculated <- which(comingfromSub2[, ncol(comingfromSub2)] > 
                               0 & comingfromSub2[, ncol(comingfromSub2)] < 1 & 
                               (is.na(comingfromSub2[, ncol(comingfromSub2)]) == 
                                  FALSE))
    comingfromSub1[thoseCalculated, ] <- comingfromSub2[thoseCalculated, 
                                                        ]
    states <- comingfromSub1
    for (i in 1:length(ancesRest)) {
      calcul <- cla_calThruNodes2(ancesRest[i], states, 
                                 loglik, forTime, parameter, use_fortran = use_fortran, 
                                 methode = methode, phy = phy,reltol=reltol,abstol=abstol,hmax=hmax)
      states <- calcul$states
      loglik <- calcul$loglik
    }
  }
  else {
    if (is.null(setting_calculation)) {
      secsse:::check_input(traits, phy, sampling_fraction, root_state_weight)
      setting_calculation <- secsse:::build_initStates_time(phy, 
                                                   traits, num_concealed_states, sampling_fraction)
    }
    states <- setting_calculation$states
    focalspeciesorder=which(phy$tip.label==focalspecies)
    if(states[focalspeciesorder,8]==1 & states[focalspeciesorder,6]==1){
      states[focalspeciesorder,8]=0
      states[focalspeciesorder,6]=1
    }
    
    if(states[focalspeciesorder,7]==1 & states[focalspeciesorder,5]==1){
      states[focalspeciesorder,7]=0
      states[focalspeciesorder,5]=1
    }
    

    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances
    if (num_concealed_states != round(num_concealed_states)) {
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 
                                                      1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    for (i in 1:length(ances)) {
      calcul <- cla_calThruNodes2(ances[i], states, loglik, 
                                 forTime, parameter, use_fortran = use_fortran, 
                                 methode = methode, phy = phy,reltol=reltol,abstol=abstol,hmax=hmax)
      states <- calcul$states
      loglik <- calcul$loglik
      nodeN <- calcul$nodeN
    }
  }
  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  mergeBranch2 <- (mergeBranch)
  if (class(root_state_weight) == "numeric") {
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights, num_concealed_states)
  }
  else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    if (root_state_weight == "proper_weights") {
      numerator <- NULL
      for (j in 1:length(mergeBranch2)) {
        numerator <- c(numerator, (mergeBranch2[j]/(sum(lambdas[[j]] * 
                                                          (1 - nodeM[1:d][j])^2))))
      }
      denomin <- NULL
      for (j in 1:length(mergeBranch2)) {
        denomin <- c(denomin, (mergeBranch2[j]/(sum(lambdas[[j]] * 
                                                      (1 - nodeM[1:d][j])^2))))
      }
      weightStates <- numerator/sum(denomin)
    }
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1/length(mergeBranch2), length(mergeBranch2))
    }
  }
  if (cond == "maddison_cond") {
    preCond <- NULL
    for (j in 1:length(weightStates)) {
      preCond <- c(preCond, sum(weightStates[j] * lambdas[[j]] * 
                                  (1 - nodeM[1:d][j])^2))
    }
    mergeBranch2 <- mergeBranch2/(sum(preCond))
  }
  if (cond == "proper_cond") {
    preCond <- NULL
    for (j in 1:length(mergeBranch2)) {
      preCond <- c(preCond, sum((lambdas[[j]] * (1 - nodeM[1:d][j])^2)))
    }
    mergeBranch2 <- mergeBranch2/preCond
  }
  atRoot <- ((mergeBranch2) * (weightStates))
  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - secsse:::penalty(pars = parameter, 
                                          loglik_penalty = loglik_penalty)
  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states), 
                               ]
    ancestral_states <- ancestral_states[, -(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, states=states, LL = LL))
  }
  else {
    return(LL)
  }
}


######################################
#Likelihood of focal species in concealed state b and phylogeny with tip traits, given the SecSSE parameters. The code is based on likelihood function in SecSSE, modified sections are highlighted with comments.  The argument focalspecies is a string defining the species of interest from phy$tip.label, other arguments are the same as cla_secsse_loglik in SecSSE.

cla_secsse_loglikb=function (parameter, phy, traits, focalspecies, num_concealed_states, use_fortran = TRUE, 
                             methode = "ode45", cond = "proper_cond", root_state_weight = "proper_weights", 
                             sampling_fraction, run_parallel = FALSE, setting_calculation = NULL, 
                             setting_parallel = NULL, see_ancestral_states = FALSE, loglik_penalty = 0, reltol=1e-12,abstol=1e-16,hmax=NULL) 
{
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  if (run_parallel == TRUE) {
    if (is.null(setting_calculation)) {
      secsse:::check_input(traits, phy, sampling_fraction, root_state_weight)
      setting_calculation <- secsse:::build_initStates_time_bigtree(phy, 
                                                                    traits, num_concealed_states, sampling_fraction)
    }
    
    #force the focal species to be in the concealed state b
    states <- setting_calculation$states
    focalspeciesorder=which(phy$tip.label==focalspecies)
    if(states[focalspeciesorder,8]==1 & states[focalspeciesorder,6]==1){
      states[focalspeciesorder,8]=1
      states[focalspeciesorder,6]=0
    }
    
    if(states[focalspeciesorder,7]==1 & states[focalspeciesorder,5]==1){
      states[focalspeciesorder,7]=1
      states[focalspeciesorder,5]=0
    }
    #end of modified code
    
    
    forTime <- setting_calculation$forTime
    ancesSub1 <- setting_calculation$ancesSub1
    ancesSub2 <- setting_calculation$ancesSub2
    ancesRest <- setting_calculation$ancesRest
    if (num_concealed_states != round(num_concealed_states)) {
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 
                                                      1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    take_ancesSub <- list(ancesSub1, ancesSub2)
    if (is.null(setting_parallel)) {
      cl <- parallel::makeCluster(2)
      doParallel::registerDoParallel(cl)
    }
    statesNEW <- cla_doParalThing(take_ancesSub, states, 
                                  loglik, forTime, parameter, use_fortran, methode, 
                                  phy)
    comingfromSub1 <- statesNEW[[1]][[1]]
    comingfromSub2 <- statesNEW[[2]][[1]]
    loglik <- statesNEW[[1]][[2]] + statesNEW[[2]][[2]]
    thoseCalculated <- which(comingfromSub2[, ncol(comingfromSub2)] > 
                               0 & comingfromSub2[, ncol(comingfromSub2)] < 1 & 
                               (is.na(comingfromSub2[, ncol(comingfromSub2)]) == 
                                  FALSE))
    comingfromSub1[thoseCalculated, ] <- comingfromSub2[thoseCalculated, 
                                                        ]
    states <- comingfromSub1
    for (i in 1:length(ancesRest)) {
      calcul <- cla_calThruNodes2(ancesRest[i], states, 
                                          loglik, forTime, parameter, use_fortran = use_fortran, 
                                          methode = methode, phy = phy,reltol=reltol,abstol=abstol,hmax=hmax)
      states <- calcul$states
      loglik <- calcul$loglik
    }
  }
  else {
    if (is.null(setting_calculation)) {
      secsse:::check_input(traits, phy, sampling_fraction, root_state_weight)
      setting_calculation <- secsse:::build_initStates_time(phy, 
                                                            traits, num_concealed_states, sampling_fraction)
    }
    states <- setting_calculation$states
    focalspeciesorder=which(phy$tip.label==focalspecies)
    if(states[focalspeciesorder,8]==1 & states[focalspeciesorder,6]==1){
      states[focalspeciesorder,8]=1
      states[focalspeciesorder,6]=0
    }
    
    if(states[focalspeciesorder,7]==1 & states[focalspeciesorder,5]==1){
      states[focalspeciesorder,7]=1
      states[focalspeciesorder,5]=0
    }
    
    
    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances
    if (num_concealed_states != round(num_concealed_states)) {
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 
                                                      1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    for (i in 1:length(ances)) {
      calcul <- cla_calThruNodes2(ances[i], states, loglik, 
                                          forTime, parameter, use_fortran = use_fortran, 
                                          methode = methode, phy = phy,reltol=reltol,abstol=abstol,hmax=hmax)
      states <- calcul$states
      loglik <- calcul$loglik
      nodeN <- calcul$nodeN
    }
  }
  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  mergeBranch2 <- (mergeBranch)
  if (class(root_state_weight) == "numeric") {
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights, num_concealed_states)
  }
  else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    if (root_state_weight == "proper_weights") {
      numerator <- NULL
      for (j in 1:length(mergeBranch2)) {
        numerator <- c(numerator, (mergeBranch2[j]/(sum(lambdas[[j]] * 
                                                          (1 - nodeM[1:d][j])^2))))
      }
      denomin <- NULL
      for (j in 1:length(mergeBranch2)) {
        denomin <- c(denomin, (mergeBranch2[j]/(sum(lambdas[[j]] * 
                                                      (1 - nodeM[1:d][j])^2))))
      }
      weightStates <- numerator/sum(denomin)
    }
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1/length(mergeBranch2), length(mergeBranch2))
    }
  }
  if (cond == "maddison_cond") {
    preCond <- NULL
    for (j in 1:length(weightStates)) {
      preCond <- c(preCond, sum(weightStates[j] * lambdas[[j]] * 
                                  (1 - nodeM[1:d][j])^2))
    }
    mergeBranch2 <- mergeBranch2/(sum(preCond))
  }
  if (cond == "proper_cond") {
    preCond <- NULL
    for (j in 1:length(mergeBranch2)) {
      preCond <- c(preCond, sum((lambdas[[j]] * (1 - nodeM[1:d][j])^2)))
    }
    mergeBranch2 <- mergeBranch2/preCond
  }
  atRoot <- ((mergeBranch2) * (weightStates))
  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - secsse:::penalty(pars = parameter, 
                                                   loglik_penalty = loglik_penalty)
  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states), 
                               ]
    ancestral_states <- ancestral_states[, -(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, states=states, LL = LL))
  }
  else {
    return(LL)
  }
}


######################################
#Function turning the cla_secsse_loglik a and b results into marginal probabilities of species being in a and b. Output is a vector of length 2 summing to 1.  The argument focalspecies is a string defining the species of interest from phy$tip.label, other arguments are the same as cla_secsse_loglik in SecSSE.

cla_secsse_propab=function (parameter, phy, traits, focalspecies, num_concealed_states=2, use_fortran = TRUE, 
                            methode = "ode45", cond = "proper_cond", root_state_weight = "proper_weights", 
                            sampling_fraction=c(1,1), run_parallel = FALSE, setting_calculation = NULL, 
                            setting_parallel = NULL, reltol=1e-12,abstol=1e-16,hmax=NULL, expprecision=64){
  
  asra=cla_secsse_loglika(parameter=parameter,
                              phy=phy, traits = traits,focalspecies = focalspecies, num_concealed_states=num_concealed_states,
                              use_fortran = use_fortran, methode = methode, cond = cond,
                              root_state_weight = root_state_weight, sampling_fraction=sampling_fraction,
                              run_parallel = run_parallel, setting_calculation = setting_calculation,
                              setting_parallel = setting_parallel, see_ancestral_states = F,
                              loglik_penalty = 0, reltol=reltol,abstol=abstol,hmax=hmax)
  
  asrb=cla_secsse_loglikb(parameter=parameter,
                              phy=phy, traits = traits,focalspecies = focalspecies, num_concealed_states=num_concealed_states,
                              use_fortran = use_fortran, methode = methode, cond = cond,
                              root_state_weight = root_state_weight, sampling_fraction=sampling_fraction,
                              run_parallel = run_parallel, setting_calculation = setting_calculation,
                              setting_parallel = setting_parallel, see_ancestral_states = F,
                              loglik_penalty = 0, reltol=reltol,abstol=abstol,hmax=hmax)
  
  proplha=as.numeric(exp(mpfr(asra,expprecision))/(exp(mpfr(asra,expprecision))+exp(mpfr(asrb,expprecision))))
  proplhb=as.numeric(exp(mpfr(asrb,expprecision))/(exp(mpfr(asra,expprecision))+exp(mpfr(asrb,expprecision))))
  
  return(c(proplha,proplhb))
} 


######################################
#Function applying the cla_secsse_propab to all the species in the phylogeny. Output is a data.frame with columns representing marginal probabilities of a and b, and rows representing species in the same order as phy$tip.labels. The argument mc.cores is a handle for using mclapply procedure, other arguments are the same as cla_secsse_loglik in SecSSE. 

cla_secsse_propab_all=function (parameter, phy, traits, num_concealed_states=2, use_fortran = TRUE, 
                            methode = "ode45", cond = "proper_cond", root_state_weight = "proper_weights", 
                            sampling_fraction=c(1,1), reltol=1e-12,abstol=1e-16,hmax=NULL,mc.cores=1){
  
  species=phy$tip.label
  
  fn=function(focalspecies){
    cla_secsse_propab(parameter, phy, traits, focalspecies,num_concealed_states=2, use_fortran = TRUE, 
                      methode = "ode45", cond = "proper_cond", root_state_weight = "proper_weights", 
                      sampling_fraction=c(1,1), reltol=reltol,abstol=abstol,hmax=hmax)
  }
  
  r=t(simplify2array(mclapply(species, fn, mc.cores=mc.cores)))
  
  res=matrix(0, length(species),4)
  
  for (i in 1:length(species)){
    res[i,traits[i]]=r[i,1]
    res[i,traits[i]+2]=r[i,2]
    
  }
  return(res)
  
}


######################################
#These are modified cla_calThruNodes and cla_secsse_loglik functions from SecSSE so that the internal ODE solver can use custom absolute tolerance.

cla_calThruNodes2=function (ances, states, loglik, forTime, parameter, use_fortran, 
                            methode, phy,  reltol=1e-12,abstol=1e-16,hmax=NULL) 
{
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  nb_node <- phy$Nnode
  ly <- ncol(states)
  d <- ncol(states)/2
  focal <- ances
  desRows <- which(phy$edge[, 1] == focal)
  desNodes <- phy$edge[desRows, 2]
  nodeM <- numeric()
  nodeN <- numeric()
  for (desIndex in 1:2) {
    y <- states[desNodes[desIndex], ]
    timeInte <- forTime[which(forTime[, 2] == desNodes[desIndex]), 
                        3]
    if (use_fortran == FALSE) {
      nodeMN <- deSolve::ode(y = y, func = secsse:::cla_secsse_loglik_rhs, 
                             times = c(0, timeInte), parms = parameter, rtol = reltol, 
                             atol = abstol, hmax = NULL, method = methode)
    }
    else {
      nodeMN <- secsse:::ode_FORTRAN(y = y, func = "cla_secsse_runmod", 
                                     times = c(0, timeInte), parms = parameter, rtol = reltol, 
                                     atol = abstol, method = methode)
    }
    if (desIndex == 1) {
      nodeN <- nodeMN
    }
    if (desIndex == 2) {
      nodeM <- nodeMN
    }
  }
  nodeM <- as.numeric(nodeM[2, -1])
  nodeN <- as.numeric(nodeN[2, -1])
  ff <- secsse:::normalize_loglik(nodeM[(1:d) + d], loglik)
  nodeM[(1:d) + d] <- ff$probs
  loglik <- ff$loglik
  ff <- secsse:::normalize_loglik(nodeN[(1:d) + d], loglik)
  nodeN[(1:d) + d] <- ff$probs
  loglik <- ff$loglik
  all_states <- cbind(nodeM[(d + 1):length(nodeM)], nodeN[(d + 
                                                             1):length(nodeN)])
  a <- cbind(all_states[, 2], all_states[, 1])
  b <- t(all_states)
  cross_M_N <- a %*% b
  mergeBranch <- 0.5 * (unlist(lapply(lapply(lambdas, "*", 
                                             cross_M_N), sum)))
  ff <- secsse:::normalize_loglik(mergeBranch, loglik)
  mergeBranch <- ff$probs
  loglik <- ff$loglik
  newstate <- nodeM[1:d]
  newstate <- c(newstate, mergeBranch)
  states[focal, ] <- newstate
  return(list(states = states, loglik = loglik, mergeBranch = mergeBranch, 
              nodeM = nodeM))
}


cla_secsse_loglik2=function (parameter, phy, traits, num_concealed_states, use_fortran = TRUE, 
                             methode = "ode45", cond = "proper_cond", root_state_weight = "proper_weights", 
                             sampling_fraction, run_parallel = FALSE, setting_calculation = NULL, 
                             setting_parallel = NULL, see_ancestral_states = FALSE, loglik_penalty = 0,reltol=1e-12,abstol=1e-16,hmax=NULL) 
{
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  Q <- parameter[[3]]
  if (run_parallel == TRUE) {
    if (is.null(setting_calculation)) {
      secsse:::check_input(traits, phy, sampling_fraction, root_state_weight)
      setting_calculation <- secsse:::build_initStates_time_bigtree(phy, 
                                                                    traits, num_concealed_states, sampling_fraction)
    }
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ancesSub1 <- setting_calculation$ancesSub1
    ancesSub2 <- setting_calculation$ancesSub2
    ancesRest <- setting_calculation$ancesRest
    if (num_concealed_states != round(num_concealed_states)) {
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 
                                                      1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    take_ancesSub <- list(ancesSub1, ancesSub2)
    if (is.null(setting_parallel)) {
      cl <- parallel::makeCluster(2)
      doParallel::registerDoParallel(cl)
    }
    statesNEW <- cla_doParalThing(take_ancesSub, states, 
                                  loglik, forTime, parameter, use_fortran, methode, 
                                  phy)
    comingfromSub1 <- statesNEW[[1]][[1]]
    comingfromSub2 <- statesNEW[[2]][[1]]
    loglik <- statesNEW[[1]][[2]] + statesNEW[[2]][[2]]
    thoseCalculated <- which(comingfromSub2[, ncol(comingfromSub2)] > 
                               0 & comingfromSub2[, ncol(comingfromSub2)] < 1 & 
                               (is.na(comingfromSub2[, ncol(comingfromSub2)]) == 
                                  FALSE))
    comingfromSub1[thoseCalculated, ] <- comingfromSub2[thoseCalculated, 
    ]
    states <- comingfromSub1
    for (i in 1:length(ancesRest)) {
      calcul <- cla_calThruNodes2(ancesRest[i], states, 
                                  loglik, forTime, parameter, use_fortran = use_fortran, 
                                  methode = methode, phy = phy,reltol=reltol,abstol=abstol,hmax=hmax)
      states <- calcul$states
      loglik <- calcul$loglik
    }
  }
  else {
    if (is.null(setting_calculation)) {
      secsse:::check_input(traits, phy, sampling_fraction, root_state_weight)
      setting_calculation <- secsse:::build_initStates_time(phy, 
                                                            traits, num_concealed_states, sampling_fraction)
    }
    states <- setting_calculation$states
    forTime <- setting_calculation$forTime
    ances <- setting_calculation$ances
    if (num_concealed_states != round(num_concealed_states)) {
      d <- ncol(states)/2
      new_states <- states[, c(1:sqrt(d), (d + 1):((d + 
                                                      1) + sqrt(d) - 1))]
      new_states <- states[, c(1, 2, 3, 10, 11, 12)]
      states <- new_states
    }
    loglik <- 0
    ly <- ncol(states)
    d <- ncol(states)/2
    for (i in 1:length(ances)) {
      calcul <- cla_calThruNodes2(ances[i], states, loglik, 
                                  forTime, parameter, use_fortran = use_fortran, 
                                  methode = methode, phy = phy,reltol=reltol,abstol=abstol,hmax=hmax)
      states <- calcul$states
      loglik <- calcul$loglik
      nodeN <- calcul$nodeN
    }
  }
  mergeBranch <- calcul$mergeBranch
  nodeM <- calcul$nodeM
  mergeBranch2 <- (mergeBranch)
  if (class(root_state_weight) == "numeric") {
    giveWeights <- root_state_weight/num_concealed_states
    weightStates <- rep(giveWeights, num_concealed_states)
  }
  else {
    if (root_state_weight == "maddison_weights") {
      weightStates <- (mergeBranch2)/sum((mergeBranch2))
    }
    if (root_state_weight == "proper_weights") {
      numerator <- NULL
      for (j in 1:length(mergeBranch2)) {
        numerator <- c(numerator, (mergeBranch2[j]/(sum(lambdas[[j]] * 
                                                          (1 - nodeM[1:d][j])^2))))
      }
      denomin <- NULL
      for (j in 1:length(mergeBranch2)) {
        denomin <- c(denomin, (mergeBranch2[j]/(sum(lambdas[[j]] * 
                                                      (1 - nodeM[1:d][j])^2))))
      }
      weightStates <- numerator/sum(denomin)
    }
    if (root_state_weight == "equal_weights") {
      weightStates <- rep(1/length(mergeBranch2), length(mergeBranch2))
    }
  }
  if (cond == "maddison_cond") {
    preCond <- NULL
    for (j in 1:length(weightStates)) {
      preCond <- c(preCond, sum(weightStates[j] * lambdas[[j]] * 
                                  (1 - nodeM[1:d][j])^2))
    }
    mergeBranch2 <- mergeBranch2/(sum(preCond))
  }
  if (cond == "proper_cond") {
    preCond <- NULL
    for (j in 1:length(mergeBranch2)) {
      preCond <- c(preCond, sum((lambdas[[j]] * (1 - nodeM[1:d][j])^2)))
    }
    mergeBranch2 <- mergeBranch2/preCond
  }
  atRoot <- ((mergeBranch2) * (weightStates))
  wholeLike <- sum(atRoot)
  LL <- log(wholeLike) + loglik - secsse:::penalty(pars = parameter, 
                                                   loglik_penalty = loglik_penalty)
  if (see_ancestral_states == TRUE) {
    num_tips <- ape::Ntip(phy)
    ancestral_states <- states[(num_tips + 1):nrow(states), 
    ]
    ancestral_states <- ancestral_states[, -(1:(ncol(ancestral_states)/2))]
    rownames(ancestral_states) <- ances
    return(list(ancestral_states = ancestral_states, states=states, LL = LL))
  }
  else {
    return(LL)
  }
}

