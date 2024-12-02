# Initial multi-state jolly-seber model for honu 
# This model only includes nesting data collected from Tern Island in 2023
# Created by: Shelbie Ishimaru
# Created on: 2024-09-09
# Last edited: 2024-11-29
################################################################################
# Preface ----------------------------------------------------------------------
# This model is a learning tool for me to begin applying Bayesian analysis to real honu data
# This model incorporates Tern Island night (nesting) survey data from the 2023 field season
# This model is highly simplified and assumes 4 states (Pre Lalo (not yet nested), nesting, internesting interval, and post Lalo)

# Initialization ---------------------------------------------------------------
library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators
library(here) #to create unbreakable file paths

# Constructing JAGS Model ------------------------------------------------------
jags.js.ms.txt <- function(){  #CHANGED FROM BOOK SINK FUNCTION
  
  #-----------------------------------
  #Parameters:
  # phi: survival probability
  # gamma: entry probability
  # p: capture probability
  #-----------------------------------
  #States (S):
  # 1 Pre Lalo
  # 2 Nesting
  # 3 Internesting Interval 
  # 4 Post Lalo 
  #Observations (O):
  # 1 seen
  # 2 not seen
  #-----------------------------------
  
  # Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi[t] <- mean.phi #survival
    gamma[t] ~ dunif(0, 1) #Prior for entry probabilities
    p[t] <- mean.p #capture
  }
  mean.phi ~ dunif(0, 1) #Prior for mean survival
  mean.p ~ dunif(0, 1) #Prior for mean capture
  
  #Define state-transition and observation matrices
  for (i in 1:M){
    #Define probabilities of state S(t+1) given S(t)
    for (t in 1:(n.occasions-1)){
      ps[1,i,t,1] <- 1-gamma[t] #prob a turtle that arrived to Lalo will not nest 
      ps[1,i,t,2] <- gamma[t] #prob a turtle that arrived at Lalo will nest
      ps[1,i,t,3] <- 0 #prob a turtle that arrived at Lalo will internest (zero b/c they didn't nest yet)
      ps[1,i,t,4] <- 0 #prob a turtle that arrived at Lalo will leave (zero b/c they didn't nest yet)
      
      ps[2,i,t,1] <- 0 #prob a turtle that nested will arrive at Lalo (0 b/c they already arrived and nested)
      ps[2,i,t,2] <- 0 #prob a turtle that nested will nest (0 b/c we know they take a break between clutches)
      ps[2,i,t,3] <- phi[t] #prob a turtle that nested will begin an internesting interval (if they survive they will have an internesting interval then nest again)
      ps[2,i,t,4] <- 1-phi[t] #prob a turtle that nested will leave Lalo (if they "die" then they leave)
      
      ps[3,i,t,1] <- 0 #prob a turtle that completed an internesting interval will arrive at Lalo (0 b/c they already arrived and nested)
      ps[3,i,t,2] <- 1 #prob a turtle that completed an internesting interval will nest (1 b/c that will always happen)
      ps[3,i,t,3] <- 0 #prob a turtle that completed an internesting interval will have another internesting interval (0 b/c they just took a break)
      ps[3,i,t,4] <- 0 #prob at turtle that completed an internesting interval will leave Lalo (0 b/c they can only leave after nesting)
      
      ps[4,i,t,1] <- 0 #prob a turtle that left lalo will arrive at Lalo (0 b/c they already arrived, nested, and left)
      ps[4,i,t,2] <- 0 #prob a turtle that left lalo will nest (0 b/c they left)
      ps[4,i,t,3] <- 0 #prob a turtle that left lalo will have another internesting interval (0 b/c they left)
      ps[4,i,t,4] <- 1 #prob at turtle that left lalo will leave Lalo (1 b/c they left and are gone for the resst of season)
      
      #Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- 0 #prob a turtle that arrived to Lalo is seen (0 b/c it has not nested so no opportunity to see it on survey)
      po[1,i,t,2] <- 1 #prob a turtle that arrived to Lalo is not seen (1 b/c it has not nested so no opportunity to see it on survey)
      
      po[2,i,t,1] <- p[t] #prob a turtle that is nesting is seen (equal to capture probability)
      po[2,i,t,2] <- 1-p[t] #prob a turtle that is nesting is not seen (equal to the inverse of the capture probability)
      
      po[3,i,t,1] <- 0 #prob a turtle is seen within an internesting interval (0 b/c it will not be seen on a night survey)
      po[3,i,t,2] <- 1 #prob a turtle is not seen within an internesting interval (1 b/c it will not be seen on a night survey)
      
      po[4,i,t,1] <- 0 #prob a turtle is seen within after leaving (0 b/c it will not be seen on a night survey)
      po[4,i,t,2] <- 1 #prob a turtle is not seen within an leaving  (1 b/c it will not be seen on a night survey)
    } 
  } 
  #Likelihood
  for (i in 1:M){
    #Define latent state at first occasion
    z[i,1] <- 1 #Make sure that all M individuals are in state 1 at t=1
    for (t in 2:n.occasions){
      #State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,]) #calculate categorical likelihood
      #Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,]) #calculate categorical likelihood
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
honu <- read.table(here("outputs", "weekly_detection.txt")) #pull in detection txt
honu <- honu[, colnames(honu) != "TURTLEID"] #remove turtle ID column
honu <- unname(honu) #remove header
rownames(honu) <- NULL
honu <- as.matrix(honu)

nz <- 500 #number of rows we want our observation matrix to be augmented

#Add dummy occasion
CH.du <- cbind(rep(0, dim(honu)[1]), honu) # add extra col at start so that honu are not seen on first day

#Augment the data
CH.du <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

CH.du #look at edited df

#Recode CH matrix: a 0 is not allowed in WinBUGS!
CH.ms <- CH.du
CH.ms[CH.ms==0] <- 2 #Not seen = 2, seen = 1
CH.ms  #final capture history matrix to be passed to data

known.state.ms <- CH.ms #create final known state z matrix to pass to data
known.state.ms[known.state.ms==2] <- NA #correctly identify all other instances
known.state.ms

my.z.init.ms <- CH.ms #create initial z matrix to pass to init
my.z.init.ms[my.z.init.ms==1] <- NA #correctly identify all other instances
my.z.init.ms

#NOT USING RN: START -------------------------------------------------------------
# Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1 START
#my.z.init <- CH.du #create z matrix out of new df with added 0 at the start

#first.one <- apply(my.z.init[,1:ncol(CH.du)], 1, function(x) min(which(x == 1))) #identify the first nesting occurrence for each honu
#last.one  <- apply(my.z.init[,1:ncol(CH.du)], 1, function(x) max(which(x == 1))) #identify the last nesting occurrence for each honu

#for(i in 1:nrow(my.z.init)) { #create function to edit each row to reflect true start and end nesting periods for each honu
  #if(first.one[i] > 1)               my.z.init[i,                1  : (first.one[i] - 1) ] = 5 #make all values before the first nest= 4 (to be changed later...)
  #if(last.one[i]  < ncol(my.z.init)) my.z.init[i, (last.one[i] + 1) : ncol(my.z.init)    ] = 4 #make all values after the last nest (correctly identify all post Lalo instances)
#} #at this point: 5= Pre Lalo, 1= Nesting, 0= Internesting, and 4= Post Lalo

#my.z.init.ms <- rbind(my.z.init, matrix(0, ncol = dim(my.z.init)[2], nrow = nz)) #create final z matrix to pass to init
#my.z.init.ms[my.z.init.ms==1] <- NA #correctly identify all nesting instances 
#my.z.init.ms[my.z.init.ms==0] <- 3 #correctly identify all internesting instances
#my.z.init.ms[my.z.init.ms==5] <- 1 #correctly identify all pre Lalo instances
#my.z.init.ms #look at z matrix. Now: 1= Pre Lalo, NA= Nesting, 3= Internesting, 4= Post Lalo
# Code to fix error from: https://groups.google.com/g/hmecology/c/S4HO-tnzep8?pli=1 END

#Augment the data
#CH.du <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz))

#known.state.ms <- CH.du #create final z matrix to pass to data
#known.state.ms[known.state.ms==1] <- 2 #correctly identify all nesting instances 
#known.state.ms[known.state.ms==0] <- NA #correctly identify all other instances
#known.state.ms

#CH.ms <- rbind(CH.du, matrix(0, ncol = dim(CH.du)[2], nrow = nz)) #create observation matrix

#Recode CH matrix: a 0 is not allowed in WinBUGS!
#CH.ms <- CH.du
#CH.ms[CH.ms==0] <- 2 #Not seen = 2, seen = 1
#CH.ms  #final capture history matrix to be passed to data
#NOT USING RN: END -------------------------------------------------------------

#Bundle data
jags.data <- list(y = CH.ms, n.occasions = dim(CH.ms)[2], M = dim(CH.ms)[1], z= known.state.ms) 

inits <- function(){list(mean.phi = runif(1, 0, 1), 
                         mean.p = runif(1, 0, 1),
                         z= cbind(rep(NA, dim(my.z.init.ms)[1]), my.z.init.ms[,-1]))} 

#Parameters monitored
params <- c("mean.p", "mean.phi", "b", "Nsuper", "N", "B")

# MCMC settings
ni <- 20000 #number of interactions
nt <- 3 #thinning rate
nb <- 5000 #burn-in length
nc <- 3 #number of chains, to check for convergence

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
