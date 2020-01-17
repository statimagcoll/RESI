#devtools::install('./')

# tests that simulations from paper run.
library(testthat)
library(pracma)
library(sandwich)
library(lmtest)
set.seed(666)
nsim=2 # low number of sims for this example
alpha=0.05
# controls skewness of gamma
shapes = c(0.5, 10)
ns = c(25, 50, 100, 250, 500, 1000)
SIs = c(0, 0.1, 0.25, 0.4, 0.6)
m1s = c(1,3,5)
m0s = c(2,5)
rhosqs = c(0, .6)
hetero = 1 # x variable that induces heterogeneity. 0 for none
out = expand.grid(shape=shapes, n=ns, m1=m1s, m0=m0s, S=SIs, rhosq=rhosqs)
params = names(out)
out[, c('bias', 'coverage.central', 'coverage.sr', 'width.central', 'width.sr')] = NA
# testing SI=SIs[1]; n=ns[1]; shape=shapes[1]; m1=m1s[1]; m0=m0s[1]; rhosq=rhosqs[2]; hetero=1
# function that collapses two compiled simulation functions.
sim.func = function(shape, n, m1, m0, S, rhosq, hetero){
  simparameters = simSetup(S=S, m1=m1, m0=m0, rhosq=rhosq, hetero=hetero)
  simFunc(simparameters, shape=shape, n=n, alpha=0.05)
}
for( rhosq in rhosqs){
  for(SI in SIs){
    for(m1 in m1s){
      for(m0 in m0s){
        for(n in ns){
          for(shape in shapes){
            cat(paste(shape, n, m1, m0, SI, rhosq, collapse=','), '\t')
            temp = t(replicate(nsim, sim.func(shape, n, m1, m0, SI, rhosq, hetero)))
            out[which( out$shape==shape & out$n==n & out$m1==m1 & out$m0==m0 & out$S==SI & out$rhosq==rhosq), c('variance', 'bias', 'coverage.central', 'coverage.sr', 'width.central', 'width.sr')] = c(var(temp[,1]), colMeans(temp))
          }
        }
      }
    }
  }
}
