# Install these packages and libraries when you open R or Rstudio 
install.packages("rjags")
install.packages("coda")
install.packages("R2jags")
install.packages("hdi")
install.packages("loo")
library(R2jags)
library(loo)

# Load data - be sure to set you working directory to the location of this file
dat <- read.csv("ex1.csv", fill = FALSE, header = TRUE) 

m1.dat<-list(c=dat$tlh,v=dat$mass,subj=dat$sample,sal=dat$salinity)

# Fitting serial dilution data to lognormal-binomial 
# The model to be input into JAGS
cat(
  "model{
  for (i in 1:38) {
  c[i] ~ dbin(p[i],3)
  p[i] <- 1-exp(-lambda[i]*v[i])
  lambda[i] ~ dlnorm(mu[i],tau_all)
  mu[i] <- b0[subj[i]]+b1*sal[i]
  }
  for (s in 1:7) {
  b0[s] ~ dnorm(mu_logb0,tau_logb0)
  }
  mu_logb0 ~ dnorm(0,0.1)
  tau_logb0 ~ dgamma(0.1,0.1)
  tau_all ~ dgamma(0.1,0.1)
  b1 ~ dnorm(0,0.1)
  }",
  file="m1.jag"
)

# Initial parameters values to feed into JAGS
m1.inits <- list(list("b0"=c(0.1,0,0.1,0,0.1,0,0.1),
                      "mu_logb0"=0,"tau_logb0"=0.1,"tau_all"=0.3,"b1"=0),
                 list("b0"=c(0,0.01,0,0.01,0,0.01,0),
                      "mu_logb0"=0,"tau_logb0"=0.2,"tau_all"=0.1,"b1"=0),
                 list("b0"=c(0.1,0.01,0.1,0.01,0.1,0.01,0.1),
                      "mu_logb0"=0,"tau_logb0"=0.3,"tau_all"=0.2,"b1"=0))

# Parameters we ask JAGS to save
parameters <- c("b1","mu_logb0","tau_logb0","tau_all","b1")

m1 <- jags(data = m1.dat,
           inits = m1.inits,
           parameters.to.save = parameters,
           model.file = "m1.jag",
           n.chains = 3,
           n.iter = 10000,
           n.burnin = 1000,
           n.thin = 3)  
