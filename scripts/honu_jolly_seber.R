# Initial multi-state jolly-seber model for honu 
# This model only includes nesting data collected from Tern Island in 2023
# Created by: Shelbie Ishimaru
# Created on: 2024-09-09
# Last edited: 2024-09-19
################################################################################
# Preface ----------------------------------------------------------------------
# This model is a learning tool for me to begin applying Bayesian analysis to real honu data
# This model incorporates Tern Island night (nesting) survey data from the 2023 field season
# This model is highly simplified and assumes only 3 states (arrival (not yet nested), nesting, and internesting interval)
# This state assumption cause finished with nesting to be lumped into internesting, so it does not account for number of clutches

# Initialization ---------------------------------------------------------------
library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators

# Constructing JAGS Model ------------------------------------------------------
jags.js.ms.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #-----------------------------------
  #Parameters:
  # phi: Renesting (persistence) probability, idea from Kendall et al. 2019 
  # gamma: removal entry probability
  # p: capture probability
  #-----------------------------------
  #States (S):
  # 1 Arrived to Lalo (not yet entered)
  # 2 Nesting on Tern
  # 3 Internesting Interval (not nesting)
  #Observations (O):
  # 1 seen
  # 2 not seen
  #-----------------------------------
  
  # Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi[t] <- mean.phi
    gamma[t] ~ dunif(0, 1) #Prior for entry probabilities
    p[t] <- mean.p
  }
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.p ~ dunif(0, 1) #Prior for mean capture
  
  #Define state-transition and observation matrices
  for (i in 1:M){
    #Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
      ps[1,i,t,1] <- 1-gamma[t] #prob a turtle that arrived to Lalo will not nest on Tern
      ps[1,i,t,2] <- gamma[t] #prob a turtle that arrived at Lalo will nest on Tern
      ps[1,i,t,3] <- 0 #prob a turtle that arrived at Lalo will have nested on Tern (zero b/c they didn't nest yet)
      ps[2,i,t,1] <- 0 #prob a turtle that nested at Tern will arrive at Lalo (0 b/c they already arrived and nested)
      ps[2,i,t,2] <- 0 #prob a turtle that nested at Tern will nest at Tern (0 b/c we know they take a break between clutches)
      ps[2,i,t,3] <- 1 #prob a turtle that nested at Tern will begin an internesting interval (1 bc we know after they nest they will take a break)
      ps[3,i,t,1] <- 0 #prob a turtle that completed an internesting interval will arrive at Lalo (0 b/c they already arrived and nested)
      ps[3,i,t,2] <- phi[t] #prob a turtle that completed an internesting interval will nest on Tern (equal to renesting persistence)
      ps[3,i,t,3] <- 1-phi[t] #prob a turtle that completed an internesting interval will not nest (equal to the inverse of the renesting persistance)
      #I feel like the last line is a bit simplified... to know if a turtle is done nesting for the season you likely will need to account for number of nests. But for the sake of learning this will do
      
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0 #prob a turtle that arrived to Lalo is seen (0 b/c it has not nested so no opportunity to see it on survey)
      po[1,i,t,2] <- 1 #prob a turtle that arrived to Lalo is seen (1 b/c it has not nested so no opportunity to see it on survey)
      po[2,i,t,1] <- p[t] #prob a turtle that is nesting on Tern is seen (equal to capture/recapture prob)
      po[2,i,t,2] <- 1-p[t] #prob a turtle that is nesting on Tern is not seen (equal to the inverse of the capture/recapture prob)
      po[3,i,t,1] <- 0 #prob a turtle is seen within an internesting interval (0 bc it will not be seen on a night survey)
      po[3,i,t,2] <- 1 #prob a turtle is not seen within an internesting interval (1 bc it will not be seen on a night survey)
    } 
  } 
  #Likelihood
  for (i in 1:M){
    #Define latent state at first occasion
    z[i,1] <- 1 #Make sure that all M individuals are in state 1 at t=1
    for (t in 2:n.occasions){
      #State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
    } 
  } 
  #Calculate derived population parameters
  for (t in 1:(n.occasions-1)){
    qgamma[t] <- 1-gamma[t]
  }
  cprob[1] <- gamma[1]
  for (t in 2:(n.occasions-1)){
    cprob[t] <- gamma[t] * prod(qgamma[1:(t-1)])
  } 
  psi <- sum(cprob[]) #Inclusion probability
  for (t in 1:(n.occasions-1)){
    b[t] <- cprob[t] / psi #Entry probability
  } 
  for (i in 1:M){
    for (t in 2:n.occasions){
      al[i,t-1] <- equals(z[i,t], 2)
    } 
    for (t in 1:(n.occasions-1)){
      d[i,t] <- equals(z[i,t]-al[i,t],0)
    } 
    alive[i] <- sum(al[i,])
  } 
  for (t in 1:(n.occasions-1)){
    N[t] <- sum(al[,t]) #Actual population size
    B[t] <- sum(d[,t]) #Number of entries
  } 
  for (i in 1:M){
    w[i] <- 1-equals(alive[i],0)
  } 
  Nsuper <- sum(w[]) #Superpopulation size
}

# Model Analysis ---------------------------------------------------------------
#Load data
honu <- as.matrix(read.table(here("outputs", "weekly_detection.txt"), sep= "", header= F))
honu <- honu %>% select(!TURTLEID)
nz <- 299
CH.aug <- rbind(honu, matrix(0, ncol= dim(honu)[2], nrow= nz))

#Add dummy occasion
CH.du <- cbind(rep(0, dim(CH.aug)[1]), CH)

# Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1 START
my.z.init <- CH.du

first.one <- apply(my.z.init[,1:ncol(CH.du)], 1, function(x) min(which(x == 1)))
last.one  <- apply(my.z.init[,1:ncol(CH.du)], 1, function(x) max(which(x == 1)))

for(i in 1:nrow(my.z.init)) {
  my.z.init[i,     first.one[i]  : last.one[i]        ] = 2
  if(first.one[i] > 1)               my.z.init[i,                1  : (first.one[i] - 1) ] = 1
  if(last.one[i]  < ncol(my.z.init)) my.z.init[i, (last.one[i] + 1) : ncol(my.z.init)    ] = 3
}
# Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1 END

#Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms[CH.ms==0] <- 2 #Not seen = 2, seen = 1

#Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1
my.z.init.ms <- rbind(my.z.init, matrix(0, ncol = dim(my.z.init)[2], nrow = nz))
my.z.init.ms[my.z.init.ms==0] <- 1

#Bundle data
jags.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1])

inits <- function(){list(mean.phi = runif(1, 0, 1), #Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1
                         mean.p = runif(1, 0, 1),
                         z = cbind(rep(NA, dim(my.z.init.ms)[1]), my.z.init.ms[,-1]))}

#Parameters monitored
params <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B")

# MCMC settings
ni <- 20000
nt <- 3
nb <- 5000
nc <- 3

#Call JAGS from R
js.ms <- jags(data  = jags.data,
              inits = inits,
              parameters.to.save = params,
              model.file = jags.js.ms.txt,
              n.chains = nc,
              n.thin= nt,
              n.iter = ni,
              n.burnin = nb)

print(js.ms, digits = 3)

k<-mcmcplots::as.mcmc.rjags(js.ms)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan
