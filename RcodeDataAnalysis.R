# Illustration of backtesting and method comparison methodologies 
# using the NASDAQ Composite index

library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library(plotrix) #color2D.matplot

source("Rfns.R")


# Data
dat = rev(read.csv("nasdaq2016.csv")$Close) # NASDAQ Composite index (Feb 8,1971 - May 18,2016)
dates = as.Date(rev(read.csv("nasdaq2016.csv")$Date))[-1] # dates for returns
x = - log(dat[-1]/dat[-length(dat)])*100 # negated percentage log returns
N=length(x) # sample size
w=500 # moving window size (same as in simulations)
(n=N-w) # out-of-sample size = 10920

# time series plot
plot(dates,x,type="l")

# risk measure levels
avec=.99 # vector of alpha levels for VaR
tvec=0.99855 # vector of tau levels for expectile
nvec=.975 # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own and in pair with ES
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

# =======================================================
# MODEL ESTIMATION AND FORECAST COMPUTATION
# =======================================================
spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1)
esmodel = NULL # expected shortfall for the assumed model/distribution
for(i in inu)
  esmodel = c(esmodel, integrate(function(x) x*dnorm(x), qmodel[i], Inf)$value/(1-VaR.levels[i]))
emodel = enorm(asy=tvec)

VaR <- matrix(nrow=n,ncol=length(VaR.levels))
VaRfhs <- matrix(nrow=n,ncol=length(VaR.levels))
VaRevt <- matrix(nrow=n,ncol=length(VaR.levels))

ES <- matrix(nrow=n,ncol=length(nvec))
ESfhs <- matrix(nrow=n,ncol=length(nvec))
ESevt <- matrix(nrow=n,ncol=length(nvec))

EXP <- matrix(nrow=n,ncol=length(tvec))
EXPfhs <- matrix(nrow=n,ncol=length(tvec))
EXPevt <- matrix(nrow=n,ncol=length(tvec))

# estimated parameters
xi <- vector(mode="numeric", length=n)
beta  <- vector(mode="numeric", length=n)
fit.par <- matrix(nrow=n, ncol=5) # 5 model parameters
mut  <- vector(mode="numeric", length=n)
sigt  <- vector(mode="numeric", length=n)


# =======================================================
# standard normal innovations
# =======================================================

for(i in 1:n)
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i)
  
  VaR[i,]=foc$VaRmodel
  VaRfhs[i,] = foc$VaRfhs
  VaRevt[i,] = foc$VaRevt
  
  EXP[i,]=foc$EXPmodel
  EXPfhs[i,] = foc$EXPfhs
  EXPevt[i,] = foc$EXPevt
  
  ES[i,]=foc$ESmodel
  ESfhs[i,] = foc$ESfhs
  ESevt[i,] = foc$ESevt
  xi[i] = foc$evt.shape
  beta[i] = foc$evt.scale
  mut[i] = foc$mut
  sigt[i] = foc$sigt
  
  if(i %% 5 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)
save(out, file="NASDAQnorm.RDATA")

# =======================================================
# skewed Student t innovations
# The version of Fernandez & Steel (1998)
# =======================================================


spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")

VaR <- matrix(nrow=n,ncol=length(VaR.levels))
VaRfhs <- matrix(nrow=n,ncol=length(VaR.levels))
VaRevt <- matrix(nrow=n,ncol=length(VaR.levels))

ES <- matrix(nrow=n,ncol=length(nvec))
ESfhs <- matrix(nrow=n,ncol=length(nvec))
ESevt <- matrix(nrow=n,ncol=length(nvec))

EXP <- matrix(nrow=n,ncol=length(tvec))
EXPfhs <- matrix(nrow=n,ncol=length(tvec))
EXPevt <- matrix(nrow=n,ncol=length(tvec))

xi <- vector(mode="numeric", length=n)
beta  <- vector(mode="numeric", length=n)
fit.par <- matrix(nrow=n, ncol=7) # 7 model parameters
mut  <- vector(mode="numeric", length=n)
sigt  <- vector(mode="numeric", length=n)

for(i in 1:n)
  
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  
  # mean and std dev of skewed t rv
  nu=coef(fit)["shape"]; ga=coef(fit)["skew"]
  
  m = mean.st(shape=nu,skew=ga)
  s = sqrt(var.st(shape=nu,skew=ga))
  
  qmodel = qskt(p=VaR.levels,df=nu,gamma=ga) # for skew t variable in standard form
  
  esmodel = NULL # expected shortfall for the assumed model/distribution
  for(j in inu)
    esmodel = c(esmodel, integrate(function(x) x*dskt(x, df=nu, gamma=ga), qmodel[j], Inf)$value/(1-VaR.levels[j]))
  esmodel=(esmodel-m)/s
  
  qmodel = (qmodel-m)/s # adjustment to have mean zero, variance one
  
  emodel = (est(asy=tvec,shape=nu,skew=ga)-m)/s
  
  foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i)
  
  VaR[i,]=foc$VaRmodel
  VaRfhs[i,] = foc$VaRfhs
  VaRevt[i,] = foc$VaRevt
  
  EXP[i,]=foc$EXPmodel
  EXPfhs[i,] = foc$EXPfhs
  EXPevt[i,] = foc$EXPevt
  
  ES[i,]=foc$ESmodel
  ESfhs[i,] = foc$ESfhs
  ESevt[i,] = foc$ESevt
  
  xi[i] = foc$evt.shape
  beta[i] = foc$evt.scale
  mut[i] = foc$mut
  sigt[i] = foc$sigt
  
  if(i %% 5 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)
save(out, file="NASDAQsstd.RDATA")

# =======================================================
# TRADITIONAL AND COMPARATIVE BACKTESTING
# =======================================================

# loading the data
load("NASDAQnorm.RDATA")
VaRout = rbind(out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
EXPout = rbind(out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
VaRout2 = rbind(out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
ESout = rbind(out$ES[,1],out$ESfhs[,1],out$ESevt[,1])
sigt.norm = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

load("NASDAQsstd.RDATA")
VaRout = rbind(VaRout,out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
EXPout = rbind(EXPout,out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
VaRout2 = rbind(VaRout2,out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
ESout = rbind(ESout,out$ES[,1],out$ESfhs[,1],out$ESevt[,1])
sigt.sstd = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

sigt.mat = rbind(sigt.norm,sigt.sstd)

# Average forecasts and method ranks
nm=dim(VaRout)[1]; n=dim(VaRout)[2]
y = tail(x,n)
plot(y, type="l")

# Matrices to store score values
smatVaR1 <-smatVaR0 <- matrix(nrow=nm,ncol=n)
smatEXP1<-smatEXP0 <- matrix(nrow=nm,ncol=n)
smatVaRES1 <- smatVaRES0 <- matrix(nrow=nm,ncol=n)


# Matrices to store p-values of two-sided calibration tests
hac=FALSE
cctVaR <- cctEXP <- cct <- matrix(nrow=nm, ncol=2)

for(i in 1:nm)
{
  smatVaR1[i,] <- sfVaR(r=VaRout[i,],x=y,a=avec,h=1)
  smatVaR0[i,] <- sfVaR(r=VaRout[i,],x=y,a=avec,h=0)
  tmp=cct.2s.VaR(x=y, r=VaRout[i,],lev=avec[1],hac=hac)
  cctVaR[i,] = c(tmp$pv.avg, tmp$pv.cond)	
  }

barVaR =  apply(VaRout, 1, mean)
ymat = matrix(rep(y, nm), nrow = nm, byrow=TRUE)
hmat = ymat > VaRout
viol = apply(hmat, 1, "sum")/n*100 # % percent violations for VaR
sbarVaR1 = apply(smatVaR1, 1, mean)
sbarVaR0 = apply(smatVaR0, 1, mean)

barEXP =  apply(EXPout, 1, mean)
sbarEXP1 = apply(smatEXP1, 1, mean)
sbarEXP0 = apply(smatEXP0, 1, mean)

barES =  apply(ESout, 1, mean)
sbarVaRES1 = apply(smatVaRES1, 1, mean)
sbarVaRES0 = apply(smatVaRES0, 1, mean)


### Traffic light matrices

tlmVaR1 <- TLMfn(r=VaRout, y=y, rm=1, lev=avec[1],h=1)
tlmVaR0 <- TLMfn(r=VaRout, y=y, rm=1, lev=avec[1],h=0)

tlmEXP1 <- TLMfn(r=EXPout, y=y, rm=2, lev=tvec[1],h=1)
tlmEXP0 <- TLMfn(r=EXPout, y=y, rm=2, lev=tvec[1],h=0)

tlmVaRES1 <- TLMfn2(r=VaRout2, r2 = ESout, y=y, lev=nvec[1],h=1)
tlmVaRES0 <- TLMfn2(r=VaRout2, r2 = ESout, y=y, lev=nvec[1],h=0)

postscript(file="plot_NASDAQ_TLMs.eps", width=6.5, height=10, horizontal=FALSE)

par(mfrow=c(3,2))

plotTLM(tlmVaR1$TLM05, rm=1, lev=avec[1])
plotTLM(tlmVaR0$TLM05, rm=1, lev=avec[1])

plotTLM(tlmEXP1$TLM05, rm=2, lev=tvec[1])
plotTLM(tlmEXP0$TLM05, rm=2, lev=tvec[1])

plotTLM(tlmVaRES1$TLM05, rm=3, lev=nvec[1])
plotTLM(tlmVaRES0$TLM05, rm=3, lev=nvec[1])

dev.off()