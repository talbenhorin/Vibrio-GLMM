rm(list=ls(all=TRUE))

# Download JAGS at https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
# Load these packages the first time your run R
#install.packages("RTools")
#install.packages("rjags")
#install.packages("coda")
#install.packages("R2jags")
#install.packages("hdi")
#install.packages("MCMCvis")

library(R2jags)
library(loo)
library(HDInterval)
library(MCMCvis)

dat <- read.csv("seagrantvibrio.csv", fill = FALSE, header = TRUE) 
m1.dat<-list(c=dat$path,v=dat$mass,samp=dat$samp)

# Mixed-effects model
# Only intercepts are random, but slopes are identical for all groups 

# Fitting MPN data to lognormal-binomial 
cat(
  "model{
    for (i in 1:5320) {
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-lambda[samp[i]]*v[i])
    }
    for (s in 1:996) {
      lambda[s] ~ dlnorm(mu[s],tau_all)
      mu[s] ~ dnorm(0,0.1)
    }
    tau_all ~ dgamma(0.1,0.1)
  }",
  file="m1.jag"
)

# Initial params BOTH YEARS
m1.inits <- list(list("mu"=numeric(996),"tau_all"=0.1),
                 list("mu"=numeric(996),"tau_all"=1),
                 list("mu"=numeric(996),"tau_all"=1))

parameters <- c("lambda")

m1 <- jags(data = m1.dat,
           inits = m1.inits,
           parameters.to.save = parameters,
           model.file = "m1.jag",
           n.chains = 3,
           n.iter = 10000,
           n.burnin = 2000,
           n.thin = 4)

# Posterior Median and Highest Posterior Density Intervals
out<-MCMCpstr(m1,
              params = parameters,
              func = median,
              type = 'summary')
write.csv(out,file = "output.csv", row.names = FALSE)
