# Assembling data on risk measure estimation for Simulation (AR(1)-GARCH(1,1) with skewed t innovations)
# This is to be used for tests and other related analyses


install.packages("rugarch")
install.packages("fGarch")
install.packages("MASS")
install.packages("expectreg")
install.packages("ismev")
install.packages("lmom")
install.packages("QRM")
install.packages("skewt")

library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")

source("Rfns.R")

n=5000 # out-of-sample size

# Levels for risk measures
avec=c(.90, .95, .99) # vector of alpha levels for VaR
tvec=c(.96561, .98761, .99855) # vector of tau levels for expectile
nvec=c(.754, .875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

# Parameters
nu=5; xi=1.5
mst = mean.st(shape=nu,skew=xi) # mean and sd for optimal forecast computations
sst = sqrt(var.st(shape=nu,skew=xi))

# recovering the mu[t] and sig[t] used in the data generation
# keeping values relevant for the optimal forecast computations
load("simdat_ms.RDATA")
mut=tail(ms$mut,n)
sigt=tail(ms$sigt,n)

# --------------------------------------------------------
# Verifying observations on which to assess forecasts
load("simdat.RDATA")
y=tail(simdat,n)

# --------------------------------------------------------
# Forecasts

#=====================================================
# VaR_alpha
#=====================================================

VaRout1 = getVaR(k=1, levels=avec,mut=mut,sigt=sigt)
VaRout2 = getVaR(k=2, levels=avec,mut=mut,sigt=sigt)
VaRout3 = getVaR(k=3, levels=avec,mut=mut,sigt=sigt)

#=====================================================
# tau-Expectile
#=====================================================

EXPout1 = getEXP(k=1, levels=tvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
EXPout2 = getEXP(k=2, levels=tvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
EXPout3 = getEXP(k=3, levels=tvec,mut=mut,sigt=sigt,mst=mst,sst=sst)

#=====================================================
# (VaR_nu, ES_nu)
#=====================================================

out=getVaRES(k=1, levels=nvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
VaRout1b = out$VaR; ESout1 = out$ES

out=getVaRES(k=2, levels=nvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
VaRout2b = out$VaR; ESout2 = out$ES

out=getVaRES(k=3, levels=nvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
VaRout3b = out$VaR; ESout3 = out$ES

#=====================================================
# --------------------------------------------------------
# Volatility estimates

load("Sim3norm3.RDATA")
sigt.norm = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

load("Sim3std3.RDATA")
sigt.std = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

load("Sim3sstd3.RDATA")
sigt.sstd = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

sigt.mat = rbind(sigt.norm,sigt.std,sigt.sstd,sigt)