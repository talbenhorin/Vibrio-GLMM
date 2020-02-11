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
m1.dat<-list(c=dat$path,v=dat$mass,tlh=dat$oysterMPN)

# Mixed-effects model
# Only intercepts are random, but slopes are identical for all groups 

# Fitting MPN data to lognormal-binomial 
cat(
  "model{
    for (i in 1:5320) {
      c[i] ~ dbin(p[i],3)
      p[i] <- 1-exp(-tlh[i]*rho*v[i])
    }
    logit(rho) <- b0
    b0 ~ dnorm(0,0.1)
  }",
  file="m1.jag"
)

# Initial params BOTH YEARS
m1.inits <- list(list("b0"=1),
                 list("b0"=-1),
                 list("b0"=0.5))

parameters <- c("b0")

m1 <- jags(data = m1.dat,
           inits = m1.inits,
           parameters.to.save = parameters,
           model.file = "m1.jag",
           n.chains = 3,
           n.iter = 5000,
           n.burnin = 2000,
           n.thin = 3)  

# Posterior Median and Highest Posterior Density Intervals
out<-MCMCpstr(m1,
              params = parameters,
              func = median,
              type = 'summary')
out95<-hdi(list(m1$BUGSoutput$sims.list$b0,m1$BUGSoutput$sims.list$b1,m1$BUGSoutput$sims.list$b2,m1$BUGSoutput$sims.list$b3,m1$BUGSoutput$sims.list$b4,m1$BUGSoutput$sims.list$b5))
# output file
med <- rbind(out[3],out[4],out[5],out[6],out[7],out[8])
lower <- rbind(out95[[1]][1,],out95[[2]][1,],out95[[3]][1,],out95[[4]][1,],out95[[5]][1,],out95[[6]][1,])
upper <- rbind(out95[[1]][2,],out95[[2]][2,],out95[[3]][2,],out95[[4]][2,],out95[[5]][2,],out95[[6]][2,])
Iwant <-data.frame(as.numeric(med), lower, upper)
write.csv(Iwant,file = "output.csv", row.names = FALSE)
