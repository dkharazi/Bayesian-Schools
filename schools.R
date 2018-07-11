## Import libraries

library(rjags)
library(ggplot2)
library(reshape2)
library(scales)

## Set the seed
set.seed(482)

# Read in the student data
studentData <- read.table("student.txt",header=T)
attach(studentData)

# Number of schools
nSchools <- length(unique(School)) 

# Number of students
nStudents <- dim(studentData)[1]    

# Number of student covariates 
# (plus 1 for the intercept) 
pStudents <- 5     

# -------------------------------------------
# MODEL 1:  Mixed-Effects Model with Student Covariates
# -------------------------------------------

# -------------
# JAGS Set-up
# -------------

# Create a data list
dataList <- list( "Y"=Y,
	                "LRT" = LRT,
	                "VR1" = VR1, 
                  "VR2"=VR2,
	                "Gender" = Gender,
	                "nStudents" = nStudents,
	                "nSchools" = nSchools,
	                "pStudents" = pStudents,
	                "School" = School)
	
# List of parameters to be monitored  
parameters <- c("beta",
                "tau2",
                "alpha",
                "tau2_alpha")
	
# Set initial values
initsValues <- list("beta" = rep(0, pStudents),
                    "tau2" = .1,
                    "alpha" = rep(0, nSchools),
                    "tau2_alpha" = .1)

# Number of iteration for "tuning" 
adaptSteps <- 5000 

# Number of iterations for "burn-in" 
burnInSteps <- 5000   

# Number of chains to run
nChains <- 2          

# Total number of iterations to save
numSavedSteps <- 10000           

# "Thinning" (1 = keep every interation)
thinSteps <- 1                  

# Iterations per chain
ITER <- ceiling((numSavedSteps * thinSteps )/ nChains) 

# -------------
# Run JAGS
# -------------

# Create, initialize, and adapt the model
jagsModel <- jags.model("schools-model1.txt", 
                        data = dataList, 
                        inits = initsValues, 
                        n.chains = nChains, 
                        n.adapt = adaptSteps)

# Burn-in the algorithm
update(jagsModel, 
       n.iter = burnInSteps)

# Run algorithm to get interations for inference
codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameters, 
                            n.iter = ITER, 
                            thin = thinSteps)

# -------------
# Look at posterior samples
# -------------

# Make a dataframe with the posterior samples
mcmcChainDF <- data.frame(as.matrix(codaSamples, 
                                    iters = T, 
                                    chains = T))

# Create a vector with the variable names
varNames <- names(mcmcChainDF)[3:(dim(mcmcChainDF)[2])]

# Number of variables
nVars <- length(varNames)

mcmcChainDF$CHAIN <- as.factor(mcmcChainDF$CHAIN)

# Construct trace plots
par(ask = T)
for( k in 1:nVars )
{
  plot_frame <- mcmcChainDF
  plot_frame$dep_var <- mcmcChainDF[ , varNames[k]]
  print(ggplot( plot_frame, 
                aes( x = ITER, 
                     y = dep_var))  +
    geom_line( aes( color = CHAIN ) ) + 
    labs( y = varNames[k] ))
  flush.console()
}

# Examine the posterior 
# Distribution of the fixed effects
par(ask = F)
varBeta <- sapply(1:5, function(i) paste0("beta.",i,"."))
postDFreshape <- melt( mcmcChainDF, 
                       id.vars = "ITER",
                       measure.vars = varBeta)
                         
ggplot( postDFreshape, 
        aes(x = variable, y = value )) +
        geom_boxplot() +
  ylab( "posterior" ) +
  xlab( "" ) 

# Examine the posterior 
# Distribution of the random effects
varAlpha <- sapply(1:nSchools, function(i) paste0("alpha.",i,"."))
postDFreshape <- melt( mcmcChainDF, 
                       id.vars = "ITER",
                       measure.vars = varAlpha)

ggplot( postDFreshape, 
        aes(x = variable, y = value )) +
  geom_boxplot() +
  scale_x_discrete( labels = 1:nSchools) +
  ylab( "posterior" ) +
  xlab( "alpha" ) 

# Rank the schools based
# on the alphas
alphaRank <- matrix(NA, numSavedSteps, nSchools)
for(k in 1:numSavedSteps){
  alphaRank[k,] <- apply(mcmcChainDF[k , varAlpha], 2, rank)
}

alphaRank <- apply(mcmcChainDF[, varAlpha], 1, rank)
avgRank <- rank(apply(alphaRank, 1, mean))

alphaRange <- apply(mcmcChainDF[, varAlpha], 
                    2, 
                    function(i) quantile(i, c(0.025,0.975)))
alphaMean <- apply(mcmcChainDF[, varAlpha], 
                    2, 
                    mean)

plot(c(-1.1,1), c(1,38), type='n', axes=F, xlab="", ylab="", main="School Rankings")
axis(1, at=seq(-1, 1, length=5))
for(j in 1:nSchools)
  {
    lines(alphaRange[ , j], rep(avgRank[j],2))
    points(alphaMean[j], avgRank[j], pch=19)
    text(-1.1, avgRank[j], j, cex=.5)
  }

# -------------------------------------------
# MODEL 2:  Mixed-Effects Model with Student and
# School Covariates
# -------------------------------------------

# -------------
# Data set-up
# -------------

# Read in the school data
schoolData <- read.table("school.txt",header=T)

attach(schoolData)

# Number of student covariates 
# (NO plus 1 -- no intercept needed) 
pSchools <- 5      

# -------------
# JAGS Set-up
# -------------

# Create a data list
dataList <- list("Y" = Y, 
                 "LRT" = LRT, 
                 "VR1" = VR1,
                 "VR2" = VR2,
                 "Gender" = Gender,
                 "nStudents" = nStudents,
                 "nSchools" = nSchools,
                 "pStudents" = pStudents,
                 "School" = School,
                 "pSchools" = pSchools,
                 "CE" = CE,
                 "RC" = RC,
                 "Other" = Other,
                 "Girls" = Girls,
                 "Boys" = Boys)

# List of parameters to be monitored  
parameters <- c("beta_st",
                "beta_sc",
                "tau2",
                "alpha",
                "tau2_alpha")

# Set initial values
initsValues <- list("beta_st" = rep(0, pStudents),
                    "beta_sc" = rep(0, pSchools),
                    "tau2"= .1, 
                    "alpha" = rep(0, nSchools), 
                    "tau2_alpha" = .1)

# Number of iteration for "tuning" 
adaptSteps <- 5000 

# Number of iterations for "burn-in" 
burnInSteps <- 5000   

# Number of chains to run
nChains <- 2          

# Total number of iterations to save
numSavedSteps <- 10000           

# "Thinning" (1 = keep every interation)
thinSteps <- 1                  

# Iterations per chain
ITER <- ceiling((numSavedSteps * thinSteps )/ nChains) 

# -------------
# Run JAGS
# -------------

# Create, initialize, and adapt the model
jagsModel <- jags.model("schools-model2.txt", 
                        data = dataList, 
                        inits = initsValues, 
                        n.chains = nChains, 
                        n.adapt = adaptSteps)

# Burn-in the algorithm
update(jagsModel, 
       n.iter = burnInSteps)

# Run algorithm to get interations for inference
codaSamples <- coda.samples(jagsModel, 
                            variable.names = parameters, 
                            n.iter = ITER, 
                            thin = thinSteps)

# -------------
# Look at posterior samples
# -------------

# Make a dataframe with the posterior samples
mcmcChainDF <- data.frame(as.matrix(codaSamples, 
                                    iters = T, 
                                    chains = T))

# Create a vector with the variable names
varNames <- names(mcmcChainDF)[3:(dim(mcmcChainDF)[2])]

# Number of variables
nVars <- length(varNames)

mcmcChainDF$CHAIN <- as.factor(mcmcChainDF$CHAIN)

# Construct trace plots
par(ask = F)
for( k in 1:nVars )
{
  plot_frame <- mcmcChainDF
  plot_frame$dep_var <- mcmcChainDF[ , varNames[k]]
  print(ggplot( plot_frame, 
                aes( x = ITER, 
                     y = dep_var))  +
          geom_line( aes( color = CHAIN ) ) + 
          labs( y = varNames[k] ))
  flush.console()
}

# Examine the posterior 
# Distribution of the student fixed effects
varBetaST <- sapply(1:5, function(i) paste0("beta_st.",i,"."))
postDFreshape <- melt( mcmcChainDF, 
                       id.vars = "ITER",
                       measure.vars = varBetaST)

par(ask = F)
ggplot( postDFreshape, 
        aes(x = variable, y = value )) +
  geom_boxplot() +
  ylab( "posterior" ) +
  xlab( "" ) 

# Examine the posterior 
# Distribution of the student school effects
varBetaSC <- sapply(1:5, function(i) paste0("beta_sc.",i,"."))
postDFreshape <- melt( mcmcChainDF, 
                       id.vars = "ITER",
                       measure.vars = varBetaSC)

ggplot( postDFreshape,  
        aes(x = variable, y = value )) +
  geom_boxplot() +
  ylab( "posterior" ) +
  xlab( "" ) 

# Examine the posterior 
# Distribution of the random effects
varAlpha <- sapply(1:nSchools, function(i) paste0("alpha.",i,"."))
postDFreshape <- melt( mcmcChainDF, 
                       id.vars = "ITER",
                       measure.vars = varAlpha)

ggplot( postDFreshape, 
        aes(x = variable, y = value )) +
  geom_boxplot() +
  scale_x_discrete( labels = 1:nSchools) +
  ylab( "posterior" ) +
  xlab( "alpha" ) 

