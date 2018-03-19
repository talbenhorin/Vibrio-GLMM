setwd("~/R/MPN")
rm(list=ls(all=TRUE))

install.packages("rjags")
install.packages("coda")
install.packages("R2jags")
install.packages("hdi")
install.packages("loo")

library(R2jags)
library(loo)

dat <- read.csv("ex1.csv", fill = FALSE, header = TRUE) 
m1.dat<-list(c=dat$tlh,v=dat$mass,subj=dat$sample,sal=dat$salinity)

# Fitting MPN data to lognormal-binomial 

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
  b0[s] ~ dnorm(mu_b0,tau_b0)
  }
  mu_b0 ~ dnorm(0,0.1)
  tau_b0 ~ dgamma(0.1,0.1)
  tau_all ~ dgamma(0.1,0.1)
  b1 ~ dnorm(0,0.1)
  }",
  file="m1.jag"
)

# Initial parameters values to feed into JAGS
m1.inits <- list(list("b0"=c(0,0,0,0,0,0,0),
                      "mu_b0"=0,"tau_b0"=0.1,"tau_all"=0.3,"b1"=0),
                 list("b0"=c(0,0,0,0,0,0,0),
                      "mu_b0"=0,"tau_b0"=0.2,"tau_all"=0.1,"b1"=0),
                 list("b0"=c(0,0,0,0,0,0,0),
                      "mu_b0"=0,"tau_b0"=0.3,"tau_all"=0.2,"b1"=0))

# Parameters we ask JAGS to save
parameters <- c("b1")

m1 <- jags(data = m1.dat,
           inits = m1.inits,
           parameters.to.save = parameters,
           model.file = "m1.jag",
           n.chains = 3,
           n.iter = 200000,
           n.burnin = 5000,
           n.thin = 3)  