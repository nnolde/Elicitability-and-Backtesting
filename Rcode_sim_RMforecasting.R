# Computation of risk measure forecasts in the simulation study
# Risk measures considered: VaR_alpha, tau-expectile and (VaR_nu, ES_nu)
# Levels of risk measures are chosen to have roughly the same magnitude under the standard normal distribution

library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")

source("Rfns.R")

# =======================================================
# SIMULATION SET-UP
# Synthetic dataset generated from an AR(1)-GARCH(1,1) process with skewed t innovations
# AR-GARCH filter parameters:
# mu=-.05; ar1 = .3 # AR(1) part
# omega=.01; al=.1; be=.85 # GARCH(1,1) parameters
# Innovation distribution parameters:
# nu=5 # shape parameter
# ga=1.5 # skewness parameter
# Burn-in period of 1000 points was used
# The simulated series is saved as "simdat.RDATA"
# =======================================================

load("simdat.RDATA")

avec=c(.90, .95, .99) # vector of alpha levels for VaR
tvec=c(.96561, .98761, .99855) # vector of tau levels for expectile
nvec=c(.754, .875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

n=5000 # out-of-sample size to evaluate forecasts

w=500 # moving window size
x=simdat[-(1:(1000-w))] # 1000 is the in-sample to accomodate moving estimation windows of up to size 1000

x=tail(simdat,n+w) # time series to be used for fitting and forecasting

# =======================================================
# Normal innovations
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

# ----------------------------------------------------

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
save(out, file="Sim3norm3.RDATA")

# =======================================================
# Student t innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="std")

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
fit.par <- matrix(nrow=n, ncol=6) # 6 model parameters
mut  <- vector(mode="numeric", length=n)
sigt  <- vector(mode="numeric", length=n)


# ----------------------------------------------------
for(i in 1:n)
{
	fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
 	fit.par[i,] = coef(fit)

	qmodel = qdist("std",p=VaR.levels,mu=0,sigma=1,shape=nu)

	esmodel = NULL # expected shortfall for the assumed model/distribution
	for(j in inu)
	esmodel = c(esmodel, integrate(function(x) x*dt(x, df=nu), qmodel[j], Inf)$value/(1-VaR.levels[j]))
	esmodel=sqrt((nu-2)/nu)*esmodel # assuming nu>2

	emodel = sqrt((nu-2)/nu)*et(asy=tvec, df=nu) # assuming nu>2

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

save(out, file="Sim3std3.RDATA")


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

save(out, file="Sim3sstd3.RDATA")
