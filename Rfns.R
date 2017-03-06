# FUNCTIONS

# ===========================================================
# Scoring functions
# ===========================================================

# Inputs: 
# r: forecast
# x: verifying observation
# a: risk measure confidence level 
# h: parameter to identify which homogeneous scoring function to use. 
## There are two choices: h=1 (b-homog. for a certain b>0), h=0 (zero-homog.)
 
# VaR-alpha

# h=1: "standard" choice, 1-homogeneous
sfVaR <- function(r,x,a,h=1)
{
  s <- vector(mode="numeric",length=length(r))
  ind1 = (x>r) & (r>0)
  ind2 = (x<r) & (r>0)
  if(h==1){
    s[ind1] <- -a*r[ind1] + x[ind1]
    s[ind2] <- (1-a)*r[ind2]
    }
	if(h==0){
	  s[ind1] <- -a*log(r[ind1]) + log(x[ind1])
	  s[ind2] <- (1-a)*log(r[ind2])
	  }
	return(s)
}

# tau-expectile

# h=1: "standard" choice, 2-homogeneous

sfEXP <- function(r,x,a,h=1)
{
  s <- vector(mode="numeric",length=length(r))
  ind1 = (x>r) & (r>0)
  ind2 = (x<r) & (r>0)
	if(h==1){
	  s[ind1] <- (-(1-2*a)*(x-r)^2 + (1-a)*r*(r-2*x))[ind1]
	  s[ind2] <- ((1-a)*r*(r-2*x))[ind2]
	  }
	if(h==0){
	  s[ind1] <- (1-2*a)*(log(x[ind1]/r[ind1]) + 1 - (x[ind1]/r[ind1])) + (1-a)*(log(r[ind1]) -1 +(x/r)[ind1])
	  s[ind2] <- (1-a)*(log(r[ind2]) -1 +(x/r)[ind2])
	  }
	return(s)
}

# VaR_nu, ES_nu
# r1: VaR forecast; r2: ES forecast

# h=1: 1/2 - homogeneous
sfVaRES <- function(r1,r2,x,a,h=1)
{
  s <- vector(mode="numeric",length=length(x))
  ind = (r2>0)
	if(h==1){s[ind]=((x[ind]>r1[ind])*(x[ind]-r1[ind]) + (1-a)*(r1[ind]+r2[ind]))/(2*sqrt(r2[ind]))}
	if(h==0){s[ind]=((x>r1)*(x-r1)/r2)[ind] + (1-a)*((r1/r2)[ind]-1+log(r2[ind]))}
	return(s)	
}



# ===========================================================
# HAC estimator for the asymptotic variance of the average relative scores
# HAC: heteroskedasticity and autocorrelation-consistent variance estimator
# using Parzen window as the lag window (set of weights)
# M: trancation point, set to approx. 2sqrt(n)	

HACestimator <- function(x)
{
	n = length(x)
	m=ceiling(2*sqrt(n))
	gam = acf(x,lag.max=m,type="covariance",plot=F)$acf
	k1 = 1:ceiling(m/2)
	k2 = (ceiling(m/2)+1):m
	
	lam = c(1, 2*(1-6*(k1/m)^2+6*(k1/m)^3),2*2*(1-k2/m)^3)
	(gam %*% lam)
}

## Computation of an HAC estimator (with Parzen window) of covariance matrix of q-vectors hV
## The procedure is time consuming, hence we do computation once and save the output
## for future use
## x: hV matrix (2 x n)

HACestimator2 <- function(x)
{
  n = dim(x)[2]; q = dim(x)[1]
  
  # Parzen weights
  m=ceiling(2*sqrt(n))
  k1 = 1:ceiling(m/2)
  k2 = (ceiling(m/2)+1):m
  lam = c(1, 2*(1-6*(k1/m)^2+6*(k1/m)^3),2*2*(1-k2/m)^3) # Parzen weights
  
  Omega0 = matrix(0, nrow=q,ncol=q)   # covariance matrix
  for(t in 1:n) Omega0 = Omega0 + (x[,t] %*% t(x[,t]))	# q x q
  
  tmp=0 # contribution due to (auto-)covariance terms
  for (t in 1:(n-m-1))
  {for (s in (t+1):(t+m)) tmp = tmp + lam[(s-t)]* x[,t]%*% t(x[,s])}
  
  Omega = Omega0 + 2*tmp	# q x q	
  Omega0 = Omega0/n; Omega=Omega/n
  return(list(Omega0=Omega0, Omega=Omega))
}

# ===========================================================
# Asymptotic variance estimates 2 pi f(0), where f: spectral density
# using Parzen lag window
# x: times series
# M: truncation point (e.g. 2 sqrt(N))
# ===========================================================

var.asymp <- function(x,M=0)
{
	if(M==0) M=2*sqrt(length(x))
	if(M %% 2!=0) M=M+1 # even M
	# autocovariance estimates
	ck = acf(x, type="covariance",plot=F, lag.max=M)$acf
	# weights based on the Parzen lag window
	kM1 = (0:(M/2))/M; kM2=((M/2+1):M)/M
	lk1 = 1-6*kM1^2+6*kM1^3; lk2 = 2*(1-kM2)^3
	lk=c(lk1,lk2)
	f = lk[1]*ck[1] + 2*sum(lk[-1]*ck[-1]) # 2*pi*f
	return(f)
}

# ===========================================================
# Tests of the fit using likelihood ratios for VaR violations
# ===========================================================

# Test of unconditional coverage
# h: vector of indicator variables H[t]=ind(r[t]<-VaR[t])
LRuc <- function(lam, h)
{
	n1 = sum(h); lamh=n1/length(h)
	mllh = sum(log(lamh*h+(1-lamh)*(1-h))) # maximum log-likelihood
	llh0 = sum(log(lam*h+(1-lam)*(1-h))) # log-likelihood under null
	LR = 2*(mllh - llh0)
	pval = pchisq(q=LR, df=1, lower.tail=F)
	return(list(LR=LR,pvalue=pval))
}

# ===========================================================
# Model-based expectiles
# ===========================================================

# Standard normal

Gfn.norm <- function(z, tau) (dnorm(z)+z*pnorm(z))/(2*dnorm(z)+z*(2*pnorm(z)-1))-tau

#(eX <- uniroot(Gfn.norm, interval=c(-4,4), tau=.99)$root)

# Student's t distribution with nu>1 degrees of freedom
# Add check nu>1, else expectile does not exist

Gfn.t <- function(z, nu, tau)
((nu+z^2)/(nu-1)*dt(z,nu)+z*pt(z,nu))/(2*(nu+z^2)/(nu-1)*dt(z,nu)+z*(2*pt(z,nu)-1))-tau

#(eX <- uniroot(Gfn.t, interval=c(-25,25), nu=4, tau=.99)$root)

# Skewed Student's t distribution with nu>1 degrees of freedom
# Add check nu>1, else expectile does not exist
# nu: shape; ga=skewness parameter

#cdf
pst <- function(x, shape=4, skew=1.5)
{
	nu=shape; ga=skew
	c=2/(ga+1/ga)
	c*ifelse(x<0,pt(ga*x,df=nu)/ga,pt(0,df=nu)*(1/ga-ga)+ga*pt(x/ga,df=nu))
}
#pdf
dst <- function(x,shape=4,skew=1.5)
{
	nu=shape; ga=skew
	c=2/(ga+1/ga)
	c*ifelse(x<0,dt(ga*x,df=nu),dt(x/ga,df=nu))
}

Gfn.st <- function(z, shape=4, skew=1)
{
	nu=shape; ga=skew
	
	mu=mean.st(shape=nu, skew=ga)
	fz=dst(z,shape=nu,skew=ga)
	Fz=pst(z,shape=nu,skew=ga)

	k=dst(0,shape=nu,skew=ga)
	Gz=ifelse(z<0,-(z^2+nu/ga^2)/(nu-1)*fz,-(z^2+nu*ga^2)/(nu-1)*fz+k*nu/(nu-1)*(ga^2-1/ga^2))
	
	a=Gz-z*Fz
	a/(2*a+z-mu)
}


est <- function(asy, shape=4, skew=1)
{
	zz = 0 * asy
    for (k in 1:length(asy)) {
    	root=function(z) Gfn.st(z,shape,skew)-asy[k]
        z = uniroot(root, interval = c(-1, 20), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}

# Variance of a skewed t random variable with location zero and scale one
mean.st <- function(shape=4, skew=1)
{
	nu=shape; ga=skew
	k=gamma((nu+1)/2)/sqrt(pi*nu)/gamma(nu/2)
	m=2*k*nu/(nu-1)*(1-2/(1+ga^2))*(ga+1/ga)
	return(m)
}

var.st <- function(shape=4, skew=1)
{
	nu=shape; ga=skew
	
	# m=mean.st(shape=nu, skew=ga)
	# v=(ga+1/ga)^2*nu/(nu-2)*(1-3*ga^2/(1+ga^2)^2) - m^2
	
	k=gamma((nu+1)/2)/sqrt(pi*nu)/gamma(nu/2)
	v=nu/(nu-2)*(1-3*ga^2/(1+ga^2)^2)-4*k^2*(nu/(nu-1))^2*(1-2/(1+ga^2))^2
	v=v*(ga+1/ga)^2
	return(v)
}


# Simulation of the asymmetric skewed t random variable
# Ref.: Zhu & Galbraith (2010)
# n: sample size
# a: skewness parameter in (0,1) (a=.5 => no skewness provided nu1=nu2)
# nu1, nu2: lower/upper tail shape parameters
rast <- function(n, a, nu1, nu2, mean=0, sd=1)
{
	u=runif(n); t1=rt(n,nu1); t2=rt(n,nu2)
	
	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)
	y=a.star*abs(t1)*(sign(u-a)-1) + (1-a.star)*abs(t2)*(sign(u-a)+1)
	(y-mean.ast(a,nu1,nu2))/sqrt(var.ast(a,nu1,nu2))*sd + mean
}

# Mean of the asymmetric skewed t random variable
# Ref: Eq.(14) in Zhu & Galbraith (2010)
mean.ast<-function(a, nu1, nu2)
{
	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)
	b=a*k1+(1-a)*k2
	4*b*(-a.star^2*nu1/(nu1-1)+(1-a.star)^2*nu2/(nu2-1))
}
# Variance of the asymmetric skewed t random variable
# Ref: Eq.(15) in Zhu & Galbraith (2010)
var.ast<-function(a, nu1, nu2)
{
	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)
	b=a*k1+(1-a)*k2
	4*(a*a.star^2*nu1/(nu1-2)+(1-a)*(1-a.star)^2*nu2/(nu2-2)) - 16*b^2*(-a.star^2*nu1/(nu1-1)+(1-a.star)^2*nu2/(nu2-1))^2
}

past<-function(x, a, nu1, nu2)
{
	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)
	2*a*ifelse(x<=0,pt(x/(2*a.star),df=nu1),pt(0,df=nu1)+(1-a)/a*(pt(x/(2*(1-a.star)),df=nu2)-pt(0,df=nu2)))	
}

qast<-function(p, a, nu1, nu2)
{
	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)
	
	2*ifelse(p<=a,a.star*qt((p/2/a),df=nu1),(1-a.star)*qt(p/2/(1-a)-a/(1-a)*pt(0,df=nu1)+pt(0,df=nu2),df=nu2))
}

dast<-function(x, a, nu1, nu2)
{
	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)
	ifelse(x<=0,a/a.star*dt(x/(2*a.star),df=nu1),(1-a)/(1-a.star)*dt(x/(2*(1-a.star)),df=nu2))	
}

Gfn.ast <- function(z, a, nu1, nu2)
{	
	mu=mean.ast(a, nu1, nu2)
	fz=dast(z,a, nu1, nu2)
	Fz=past(z,a, nu1, nu2)

	k1=dt(0,df=nu1); k2=dt(0,df=nu2)
	a.star=a*k1/(a*k1+(1-a)*k2)

	b=a*k1+(1-a)*k2
	k=dast(0,a, nu1, nu2)
	
	Gz=-4*ifelse(z<=0,a.star^2*nu1/(nu1-1)*(1+(z/2/a.star)^2/nu1)*fz,(1-a.star)^2*nu2/(nu2-1)*(1+(z/2/(1-a.star))^2/nu2)*fz+b*(a.star^2*nu1/(nu1-1)-(1-a.star)^2*nu2/(nu2-1)))
	
	a=Gz-z*Fz
	a/(2*a+z-mu)
}

east <- function(asy, a, nu1, nu2)
{
	zz = 0 * asy
    for (k in 1:length(asy)) {
    	root=function(z) Gfn.ast(z,a, nu1, nu2)-asy[k]
        z = uniroot(root, interval = c(-10, 10), tol = 1e-06)
        zz[k] = z$root
    }
    return(zz)
}


# ===========================================================
# EVT-based expectiles
# ===========================================================
# tau: expectile level
# Interest is in the upper tail (i.e. tau >.5)
# k: the order of the upper order statistic used as a threshold
# n: sample size
# u: threshold used to fit GPD to excesses
# beta, xi: fitted values of the GP distribution
# par: list with values of k, n, u, beta, xi
Gfn.evt <- function(z, tau, dat, param)
{	
	k=param$k; n=param$n; u=param$u; beta=param$beta; xi=param$xi
		
	a=1+xi/beta*(z-u)
	num = k/n*beta/(1-xi)*a^(-1/xi+1) # numerator of the omega ratio
	
	zbaru = sum(dat[dat<=u])/n # denominator of the omega ratio
	c = zbaru +k/n*(u+beta/(1-xi))
	den=z + k/n*a^(-1/xi)*(xi/(1-xi)*(z-u)+beta/(1-xi)-u) - c
	
	omega = num/den
	G = 1/(1+omega) - tau
	
	return(G)	
}

expectile.evt <- function(x,tau,k,u,beta,xi)
{
	eEVT=expectile(as.vector(x,"numeric"), probs=tau) #empirical estimates

#	eEVT=expectile(x, probs=tau) #empirical estimates
	e.emp=eEVT

	# list of parameters for the EVT estimator
	param=list(k=k,n=length(x),u=u,beta=beta,xi=xi)

	# lowest tau level at which the EVT estimator works; otherwise, the empirical estimator can be used
	tau.thresh = Gfn.evt(u, tau=.9, dat=x, param=param) + .9

	indEVT=(tau>tau.thresh) # indicator whether an EVT estimator is used 

	if(xi>=0){
		for (i in which(indEVT))
	{
		ubnd<- eEVT[i]*4 # an upper bound on expectile values 
		eEVT[i]<-uniroot(Gfn.evt, interval=c(u,ubnd), tau=tau[i], dat=x, param)$root
	}} else	{
		ubx<- 0.99999*(u-beta/xi) # an upper bound on expectile values 
		for (i in which(indEVT)){
			ubnd<- min(eEVT[i]*2,ubx)
			eEVT[i]<-uniroot(Gfn.evt, interval=c(u,ubnd), tau=tau[i], dat=x, param)$root	
		} 
	}
	
	return(list(expectile=eEVT,ind=indEVT, emp=e.emp))		
}

# ===========================================================
# RISK MEASURE FORECASTS
# Risk measures considered: VaR_alpha, expectile_tau and (VaR_nu,Expected Shortfall_nu)
# Forecasting methods: 
# (1) model-based forecasts (using a functional of the assumed error distribution)
# (2) Filtered historical simulation:
# Nonparametric, sampling with replacement is used (uses "seed" for reproducibility)
# returns 1-period ahead risk measure forecasts based on the empirical distribution of the simulated returns
# (3) EVT (extreme value theory based forecasts)

# INPUTS:
# - fit: ugarchfit object "fit" (rugarch package) (after fitting negated log-returns = "log-losses")
# - alpha: VaR level, can be vectorized
# - tau: expectile level, >1/2 (upper tail), can be vectorized
# - nu: level for pair (VaR_nu, ES_nu), can be vectorized
# - qmodel: model-based quantiles corresponding to alpha level(s)
# - emodel: model-based expectile(s) corresponding to tau level(s)
# - esmodel: model-based expected shortfall corresponding to nu level(s)
# - seed: used for non-parameteric sampling for FHS
# - n: sample size used in FHS

# OUTPUTS: 
# - list of 1-period-ahead risk measure forecasts of negated log-return series under the considered estimation methods
# - Parameter estimates when fitting the GPD distribution to threshold excesses
# - Forecasts of conditional mean and volatility mu[t+1] and sigma[t+1]
# ===========================================================

RM.forecasts3 <- function(fit,alpha,tau,nu,qmodel,emodel,esmodel,n=10^5,seed)
{
  ia = 1:length(alpha) #index set for alpha levels of VaR
  inu = (length(alpha)+1):(length(alpha)+length(nu)) #index set for nu levels of VaR
  a = c(alpha,nu) # combined levels for VaR computation
  
  resid=residuals(fit) # residuals = returns - fitted values
  sigma.t=sigma(fit) # conditional volatility estimates
  z=resid/sigma.t # standardized residuals (representing iid errors)
  z=as.numeric(z)
  
  frcst = ugarchforecast(fit,n.ahead=1)
  sigt1 = sigma(frcst)
  mut1 = fitted(frcst)
  
  # sampling with replacement from the standardized residuals
  set.seed(seed)
  zt1=sample(z, size=n, replace=TRUE)
  
  qFHS = quantile(zt1,a) # VaR_alpha and VaR_nu
  qFHSnu = quantile(zt1,nu) # VaR_nu, to be paired with ES_nu
  esFHS = NULL
  for (i in 1:length(nu)) esFHS = c(esFHS, mean(zt1[zt1>= qFHSnu[i]]))
  
  
  # quantile based on EVT-POT method with threshold X[(k+1)] for k=100
  # See McNeil & Frey(2000) for discussion of threshold choice
  k=.12*length(z) # upper order statistic used as a threshold
  u=sort(z,decreasing=T)[k+1] # threshold
  # Fit GPD to threshold excesses
  GPDfit = gpd.fit(z,threshold=u,show=F)
  beta=GPDfit$mle[1]; xi=GPDfit$mle[2] # MLE's of scale and shape parameters
  
  qPOT = u + beta/xi*(((1-a)/(k/length(z)))^(-xi)-1)
  esPOT = qPOT[inu]*(1/(1-xi) + (beta - xi*u)/(1-xi)/qPOT[inu])
  
  out.evt = expectile.evt(x=z,tau=tau,k=k,u=u,beta=beta,xi=xi)
  eEVT=out.evt$expectile
  
  # alpha- and nu- VaR forecasts for the next period
  VaRmodel=mut1 + sigt1*qmodel
  VaRfhs=mut1 + sigt1*qFHS
  VaRevt=mut1 + sigt1*qPOT
  
  # tau-expectile forecasts for the next period
  e.model= mut1 + sigt1*emodel
  e.fhs  = mut1 + sigt1*expectile(zt1,probs=tau)
  e.evt  = mut1 + sigt1*eEVT
  
  # nu- ES forecasts for the next period
  es.model= mut1 + sigt1*esmodel
  es.fhs  = mut1 + sigt1*esFHS
  es.evt  = mut1 + sigt1*esPOT
  
  return(list(VaRmodel=VaRmodel,VaRfhs=VaRfhs,VaRevt=VaRevt,
              EXPmodel=e.model,EXPfhs=e.fhs,EXPevt=e.evt,ESmodel=es.model, ESfhs=es.fhs, ESevt=es.evt,evt.shape=xi, evt.scale=beta, mut=mut1, sigt=sigt1))	
}



# ===========================================================
# CONDITIONAL CALIBRATION TESTS (CCT)
# ===========================================================
# Two-sided Calibration Tests
# Inputs:
# x: verifying observations
# r: risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# lev: risk measure level ("alpha" for VaR, "tau" for expectile, "nu" for pair (VaR,ES))

# Output:
# p-values for two tests: average and conditional two-sided calibration tests

cct.2s.VaR <- function(x, r, lev=0.95, hac=FALSE)
{
	
	n = length(x) # out-of-sample size
	 
	# setting up the identification function
	v = (r>=x) - lev

	# average calibration test
	
	if(hac) Omega = HACestimator(v) #an HAC estimator of the asymptotic variance should be used for the average calibration tests (see GW2006, p.1557)
	else{ Omega = (v %*% v)/n  # zero truncation lag
	}
	
	tn.avg <- sqrt(n)*mean(v)/sqrt(Omega) 	#test statistic
	pv.avg <- 2*pnorm(-abs(tn.avg)) # p-value	
	
	# This test function leads to results similar to the average calib. test for NASQAQ data	
	N=n-1
	vt=v[-1] #v[t], t=2,...,n, identification function used in test statistic computation
	h=rbind(rep(1,N),v[-n])

	q=dim(h)[1]
		
	zbar=(h %*% (vt))/N # q x 1
	Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
	for (t in 1:N)
	{			
		 Omega = Omega + vt[t]^2*(h[,t]%*%t(h[,t]))	
	  
	}
	Omega = Omega/N
	Omega.inv = solve(Omega)
	
	tn.cond <- N*t(zbar)%*% Omega.inv %*% zbar 	#test statistic
	pv.cond <- pchisq(tn.cond, df=q, lower.tail=FALSE) # p-value	

	return(list(tn.avg = tn.avg, pv.avg = pv.avg, tn.cond = tn.cond, pv.cond = pv.cond))
	
}

cct.2s.EXP <- function(x, r, lev=0.95, sigt, hac=FALSE)
{
	
	n = length(x) # out-of-sample size
	 
	# setting up the identification function
	v = abs((r>=x)-lev)*(r-x)

	# average calibration test
	if(hac) Omega = HACestimator(v) #an HAC estimator of the asymptotic variance should be used for the average calibration tests (see GW2006, p.1557)
	else{ Omega = (v %*% v)/n  # zero truncation lag
	}
	
	tn.avg <- sqrt(n)*mean(v)/sqrt(Omega) 	#test statistic
	pv.avg <- 2*pnorm(-abs(tn.avg)) # p-value	
	
	# conditional calibration tests	

	h=rbind(rep(1,n),sigt^(-1))
	q=dim(h)[1]
		
	zbar=(h %*% (v))/n # q x 1
	Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
	for (t in 1:(n-1))
	{			
		Omega = Omega + v[t]^2*(h[,t]%*%t(h[,t]))	
	}
	Omega = Omega/(n)
	Omega.inv = solve(Omega)
	
	tn.cond <- (n)*t(zbar)%*% Omega.inv %*% zbar 	#test statistic
	pv.cond <- pchisq(tn.cond, df=q, lower.tail=FALSE) # p-value	

	return(list(pv.avg = pv.avg, pv.cond = pv.cond))
	
}

# Two-sided Conditional Calibration Tests (CCT) for k-variate risk measure with k>1
# The code is specific for (VaR,ES) pair
# Inputs:
# x: verifying observations
# r: risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# lev: risk measure level ("nu" for pair (VaR,ES))
# htype: type of the test function matrix h 

# Test functions
# 1: ht = (0,1)
# 2: ht = cbind(sigt^(-1)*c((r[2,t]-r[1,t])/(1-nu),1))


# Output:
# p-values for the cct tests corresponding to the above mentioned choices of the test functions

cct.twosided2 <- function(x, r, lev=0.95, sigt, hac=FALSE)
{
	n = length(x) # out-of-sample size
	 
	# setting up appropriate identification functions for the risk measure
	ind = (x>r[,1])
	v1 = 1 - lev - ind
	v2 = r[,1] - r[,2] - 1/(1-lev)*ind*(r[,1]-x)
	v = rbind(v1,v2) # identification function matrix (k x n)

	# average calibration test using test statistic tau3 (with h=identity or h=(0,1))
	  vbar = apply(v,1,"mean")
	  Omega = matrix(0, nrow=2,ncol=2)   # covariance matrix
	  for(t in 1:n) Omega = Omega + (v[,t] %*% t(v[,t]))	# k x k
	  Omega = Omega/n
	  Omega.inv=solve(Omega)
	  tn = n*t(vbar)%*% Omega.inv %*% vbar 	#test statistic
	  pv.avg <- pchisq(tn, df=2, lower.tail=FALSE) # p-value
	  
	# conditional calibration tests	
		q=1; N=n
		zsum = matrix(0, nrow=q,ncol=1) # storage variable
		Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix

		for (t in 1:N)
		{
			ht = sigt[t]^(-1)*c((r[t,2]-r[t,1])/(1-lev),1) # q x k (i.e., 1 x 2)
			zt = ht %*% v[,t]	# q x 1
			zsum = zsum + zt
			Omega = Omega + (zt %*% t(zt))	# q x q
		}
		zbar = zsum/(N)
		Omega = Omega/(N)
		Omega.inv = solve(Omega)
		
		tn <- N*t(zbar)%*% Omega.inv %*% zbar 	#test statistic
		pv.cond <- pchisq(tn, df=q, lower.tail=FALSE) # p-value
	
	return(list(pv.avg = pv.avg, pv.cond=pv.cond))
}

# ===================================================================
# ONE-SIDED TEST of super-calibration
# test statistic tau2 for average tests and tau4 for conditional tests
# ===================================================================

cct.1s.VaR <- function(x, r, lev=0.95, hac=FALSE) # super-calbration
{
	n = length(x) # out-of-sample size
	 
	# setting up appropriate identification functions for the risk measure
	v = (r>=x) - lev

	# average calibration test
	if(hac) Omega = HACestimator(v)
	else{ Omega = (v %*% v)/n  # zero truncation lag
	}

  tn.avg <- sqrt(n)*mean(v)/sqrt(Omega) 	#test statistic
  pv.avg <- pnorm(tn.avg) # p-value	
  	
	# conditional calibration tests	

	## same test function as for the two-sided test
	N=n
	vt=v
	h=rbind(rep(1,N),abs(r))
	
	q=dim(h)[1]

	zbar=(h %*% (vt))/N # q x 1
	Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
	for (t in 1:N)
	{		
		Omega = Omega + vt[t]^2*(h[,t]%*%t(h[,t]))	# V is a scalar for k=1
	}
	Omega = Omega/N
	
	# Hommel's (1983) procedure
	tn <- sqrt(N) * diag(Omega)^(-1/2) * zbar 	#test statistic
	pi <- sort(pnorm(tn))
	cq <- sum(1/(1:q))
	pv.cond <- min(q*cq*min(pi/(1:q)),1)
	
	pv.condB <- min(q*pi[1],1) # q*min(p_i) Bonferroni multiple test correction
	
	return(list(pv.avg = pv.avg, pv.cond = pv.cond, pv.condB = pv.condB))
}

cct.1s.EXP <- function(x, r, lev=0.95, sigt, hac=FALSE) # super-calbration
{
	n = length(x) # out-of-sample size
	 
	# setting up appropriate identification functions for the risk measure
	v = abs((r>=x)-lev)*(r-x)

	# average calibration test
	if(hac) Omega = HACestimator(v)
	else{ Omega = (v %*% v)/n  # zero truncation lag
	}
  

  tn.avg <- sqrt(n)*mean(v)/sqrt(Omega) 	#test statistic
  pv.avg <- pnorm(tn.avg) # p-value	
  
	# conditional calibration tests	
	
	h=sigt^(-1)
	q=1

	zbar=(h %*% (v))/(n) # q x 1
	Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
	for (t in 1:(n))
	{		
		Omega = Omega + v[t]^2*(h[t]%*%t(h[t]))	# V is a scalar for k=1
	}
	Omega = Omega/(n)

	# Hommel's (1983) procedure
	tn <- sqrt(n) * diag(Omega)^(-1/2) * zbar 	#test statistic
	pi <- sort(pnorm(tn))
	cq <- sum(1/(1:q))
	pv.cond <- min(q*cq*min(pi/(1:q)),1)
	
	return(list(pv.avg = pv.avg, pv.cond = pv.cond))
}



# One-sided Conditional Calibration Tests for k-variate risk measure with k>1
# The code is specific for (VaR,ES) pair
# Inputs:
# x: verifying observations
# r: risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# lev: risk measure level ("nu" for pair (VaR,ES))
# htype: type of the test function matrix h 

# Test functions
# 1: ht=diag(2) # identity matrix for the unconditional test
# 2: ht = cbind(c(1, abs(V1[t-1]), r1[t],0,0),c(0,0,0,1,sigt^(-1)) )
# 3: ht = cbind(c(1, abs(V1[t-1]), r1[t],0),c(0,0,0,sigt^(-1)) )

# Output:
# pv: p-value (minimum of the adj. p-values) for the cct test based on the Hommel's (1983) procedure
#     corresponding to one of the above mentioned choices of the test functions

cct.onesided2 <- function(x, r, lev=0.95, sigt, hac=FALSE) # sub-calibration
{
	n = length(x) # out-of-sample size
	 
	# setting up appropriate identification functions for the risk measure
	ind = (x>r[,1])
	v1 = 1 - lev - ind
	v2 = r[,1] - r[,2] - 1/(1-lev)*ind*(r[,1]-x)
	v = rbind(v1,v2) # identification function matrix (k x n)
	
	# average calibration test

	# average calibration test using test statistic tau2
	vbar = apply(v,1,"mean")
	
	if(hac) Omega = HACestimator2(v)$Omega
	else{
	  Omega = matrix(0, nrow=2,ncol=2)   # covariance matrix
	  for(t in 1:n) Omega = Omega + (v[,t] %*% t(v[,t]))	# k x k
	  Omega = Omega/n  
	}
	
	
	# Hommel's (1983) procedure
	tn <- sqrt(n) * diag(Omega)^(-1/2) * vbar 	#test statistic
	pi <- sort(1-pnorm(tn))
	cq <- sum(1/(1:2))
	pv.avg <- min(2*cq*min(pi/(1:2)),1) # q=2=k, dim of risk measure
	pv.avgB <- min(2*pi[1],1) # q*min(p_i) Bonferroni multiple test correction
	
	
	# conditional calibration tests	
	# ht = cbind(c(1, abs(r1[t]),0,0),c(0,0,1,sigt^(-1)) )
		q=4; N=n
		zsum = matrix(0, nrow=q,ncol=1) # storage variable
		Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
		
		for (t in 1:N)
		{
			ht = cbind(c(1, abs(r[t,1]),0,0),c(0,0,1,sigt[t]^(-1)) ) # q x k with q=q1+q2
			zt = ht %*% v[,t]	# q x 1 
			zsum = zsum + zt
			Omega = Omega + (zt %*% t(zt))	# q x q	
		}
	zbar = zsum/N
	Omega = Omega/N
	
	# Hommel's (1983) procedure
	
	tn <- sqrt(N) * diag(Omega)^(-1/2) * zbar 	#test statistic
	pi <- sort(1-pnorm(tn))
	cq <- sum(1/(1:q))
	pv.cond <- min(q*cq*min(pi/(1:q)),1)
	pv.condB <- min(q*pi[1],1) # q*min(p_i) Bonferroni multiple test correction
	
	return(list(pv.avg = pv.avg, pv.cond=pv.cond, pv.avgB = pv.avgB, pv.condB=pv.condB))
}


# ===================================================================
# Traffic light matrix 
# ===================================================================
# The function computes TLM's for three significance levels (1%, 5%, 10%)
# Inputs: 
# r: (nm x n) matrix of forecasts (n = out-of-sample size)
# y: (n) verifying observations
# rm: risk measure (rm=1: VaR, rm=2: expectile)
# lev: risk measure level (alpha for VaR, tau for expectile, nu for (VaR, ES))
# h: indicates which score function to use. 
## There are two choices available: h=1 ("standard" or b-homogeneous) and h=0 (0-homog.)

TLMfn <- function(r, y, rm, lev, h)
{
	m=dim(r)[1] # number of forecasting methods
	pvmat = matrix(nrow=m, ncol=m) # p-values for H0^+ hypothesis
	pvmatL = matrix(nrow=m, ncol=m) # p-values for H0^- hypothesis 
	
	if(rm==1)
	{
		for(j in 1:m)
		{
			i0 <- j
			xs <- r[,i0]
			tmpS <- sfVaR(r=r[i0,],x=y,a=lev,h=h)
	

			for(k in (1:m)[-(i0)])
			{
				tmp <- sfVaR(r=r[k,],x=y,a=lev,h=h)
				test=efp.test(tmp,tmpS,type="one-sided-ge") # H0^+
				pvmat[j,k] <- test$pvalue
				test=efp.test(tmp,tmpS,type="one-sided-le") # H0^-
				pvmatL[j,k] <- test$pvalue
			}}}
				
	if(rm==2)
	{
		for(j in 1:m)
		{
			i0 <- j
			xs <- r[,i0]
			tmpS <- sfEXP(r=r[i0,],x=y,a=lev,h=h)
			
			for(k in (1:m)[-(i0)])
			{
				tmp <- sfEXP(r=r[k,],x=y,a=lev,h=h)
				test=efp.test(tmp,tmpS,type="one-sided-ge") # H0^+
				pvmat[j,k] <- test$pvalue
				test=efp.test(tmp,tmpS,type="one-sided-le") # H0^-
				pvmatL[j,k] <- test$pvalue
			}}}

	TLM01 <- (pvmat <= 0.01) + (pvmatL > 0.01); diag(TLM01) <- 3 
	TLM05 <- (pvmat <= 0.05) + (pvmatL > 0.05); diag(TLM05) <- 3  
	TLM10 <- (pvmat <= 0.10) + (pvmatL > 0.10); diag(TLM10) <- 3  
	return(list(TLM01=TLM01, TLM05=TLM05, TLM10=TLM10, pvmat=pvmat, pvmatL=pvmatL))

}

# TLM for pair (VaR, ES)
# r1: nm x n matrix of VaR forecasts
# r2: nm x n matrix of ES forecasts

TLMfn2 <- function(r1, r2, y, lev=0.95,h=1)
{
	m=dim(r1)[1] # number of forecasting methods
	pvmat = matrix(nrow=m, ncol=m) # p-values for H0^+ hypothesis
	pvmatL = matrix(nrow=m, ncol=m) # p-values for H0^- hypothesis 

	for(j in 1:m)
	{
		i0 <- j
		tmpS <- sfVaRES(r1=r1[i0,],r2=r2[i0,],x=y,a=lev,h=h)
			
		for(k in (1:m)[-(i0)])
		{
			tmp <- sfVaRES(r1=r1[k,],r2=r2[k,],x=y,a=lev,h=h)
			test=efp.test(tmp,tmpS,type="one-sided-ge") # H0^+
			pvmat[j,k] <- test$pvalue
			test=efp.test(tmp,tmpS,type="one-sided-le") # H0^-
			pvmatL[j,k] <- test$pvalue
	}}
	
	TLM01 <- (pvmat <= 0.01) + (pvmatL > 0.01); diag(TLM01) <- 3 
	TLM05 <- (pvmat <= 0.05) + (pvmatL > 0.05); diag(TLM05) <- 3  
	TLM10 <- (pvmat <= 0.10) + (pvmatL > 0.10); diag(TLM10) <- 3  
	return(list(TLM01=TLM01, TLM05=TLM05, TLM10=TLM10, pvmat=pvmat, pvmatL=pvmatL))
}


# Plotting of TLM 
# Input: x is a (nm by nm) TLM with rows corresponding to the standard model and column to the internal model
# The diagonal of TLM is filled with dummy value "3"

plotTLM <- function(x, rm, lev)
{
	
	m = dim(x)[1] # number of methods considered
	cellcol<-matrix(nrow=m)
	cellcol[x==3]<-"white"
	cellcol[x==0]<-"red"
	cellcol[x==1]<-"yellow"
	cellcol[x==2]<-"green"

	if(rm==1) xlb = bquote(alpha == .(lev))
	if(rm==2) xlb = bquote(tau == .(lev))
	if(rm==3) xlb = bquote(nu == .(lev))

	
	par(mar=c(4,6,6,2)+0.1)
	color2D.matplot(x,cellcolors=cellcol, axes=FALSE, xlab="", ylab="")
	axis(side=3, at=(1:m)-0.5,labels=method, cex.axis=0.9, las=2)
	axis(side=2, at=(1:m)-0.5,labels=rev(method), cex.axis=0.9, las=2)
	mtext("internal model", side=3, line=4, cex=0.9)
	mtext("standard model", side=2, line=4, cex=0.9)
	mtext(xlb, side=1, line=1, cex=0.8)	
}

# ===================================================================
# Extracting and combining RM forecasts under different estimation methods for a given confidence level
# Input: index "k" for the RM level in "levels"
# See Rcode_sim_RMdataprep.R for details
# ===================================================================

getVaR <- function(k, levels, mut, sigt, nu=5, xi=1.5)
{
VaRa=qsstd(p=levels[k],nu=nu,xi=xi)  # a-VaR for the innovation distribution
VaRopt=mut + sigt*VaRa
	
load("Sim3norm3.RDATA")
VaRout = rbind(out$VaR[,k],out$VaRfhs[,k],out$VaRevt[,k])

load("Sim3std3.RDATA")
VaRout = rbind(VaRout,out$VaR[,k],out$VaRfhs[,k],out$VaRevt[,k])

load("Sim3sstd3.RDATA")
VaRout = rbind(VaRout,out$VaR[,k],out$VaRfhs[,k],out$VaRevt[,k])

VaRout=rbind(VaRout,VaRopt)

return(VaRout)
}


getEXP <- function(k, levels,mut,sigt,mst,sst,nu=5,xi=1.5)
{
EXPa = (est(asy=levels[k],shape=nu,skew=xi)-mst)/sst # a-expectile for the innovation distribution
EXPopt=mut + sigt*EXPa
	
load("Sim3norm3.RDATA")
EXPout = rbind(out$EXP[,k],out$EXPfhs[,k],out$EXPevt[,k])

load("Sim3std3.RDATA")
EXPout = rbind(EXPout,out$EXP[,k],out$EXPfhs[,k],out$EXPevt[,k])

load("Sim3sstd3.RDATA")
EXPout = rbind(EXPout,out$EXP[,k],out$EXPfhs[,k],out$EXPevt[,k])

EXPout=rbind(EXPout,EXPopt)

return(EXPout)
}

getVaRES <- function(k,levels,mut,sigt,mst,sst,nu=5,xi=1.5)
{
VaRa = qskt(p=levels[k],df=nu,gamma=xi) # for skew t variable in standard form of the density
ESa=integrate(function(x) x*dskt(x, df=nu, gamma=xi), VaRa, Inf)$value/(1-levels[k])

VaRa = (VaRa-mst)/sst
ESa=(ESa-mst)/sst # adjustment for skew t distribution with mean zero and sd=1
VaRopt=mut + sigt*VaRa
ESopt = mut + sigt*ESa

	
load("Sim3norm3.RDATA")
VaRout = rbind(out$VaR[,(3+k)],out$VaRfhs[,(3+k)],out$VaRevt[,(3+k)])
ESout = rbind(out$ES[,k],out$ESfhs[,k],out$ESevt[,k])

load("Sim3std3.RDATA")
VaRout = rbind(VaRout,out$VaR[,(3+k)],out$VaRfhs[,(3+k)],out$VaRevt[,(3+k)])
ESout = rbind(ESout,out$ES[,k],out$ESfhs[,k],out$ESevt[,k])

load("Sim3sstd3.RDATA")
VaRout = rbind(VaRout,out$VaR[,(3+k)],out$VaRfhs[,(3+k)],out$VaRevt[,(3+k)])
ESout = rbind(ESout,out$ES[,k],out$ESfhs[,k],out$ESevt[,k])

VaRout=rbind(VaRout,VaRopt)
ESout=rbind(ESout,ESopt)

return(list(VaR=VaRout,ES=ESout))
}


# ===================================================================
# Summary tables of risk measure forecasts and method comparisons based on the sample average of consistent scores.
# Inputs: tables with risk measure estimates; 
# y: vector of verifying observations
# aa,tt,nn: rm levels (alpha, tau and nu, resp)
# ===================================================================

rmf_comp <- function(VaRout, EXPout,VaRoutb, ESout, y, aa, tt, nn)
{
  nm=dim(VaRout)[1]; n=dim(VaRout)[2]
  
  # for LaTex output
  method=c("N-N   ", "N-FHS ","N-EVT ","t-t   ","t-FHS ","t-EVT ","ST-ST ","ST-FHS","ST-EVT","opt   ")
  lbr = rep("(",10); rbr = rep(")",10)
  and = rep("&",10); and3 = rep(and,3)
  
  # Matrices to store score values
  smatVaR1 <-smatVaR0 <- matrix(nrow=nm,ncol=n)
  smatEXP1<-smatEXP0 <- matrix(nrow=nm,ncol=n)
  smatVaRES1 <- smatVaRES0 <- matrix(nrow=nm,ncol=n)
  
  for(i in 1:nm)
  {
    smatVaR1[i,] <- sfVaR(r=VaRout[i,],x=y,a=aa,h=1)
    smatVaR0[i,] <- sfVaR(r=VaRout[i,],x=y,a=aa,h=0)
    
    smatEXP1[i,] <- sfEXP(r=EXPout[i,],x=y,a=tt,h=1)
    smatEXP0[i,] <- sfEXP(r=EXPout[i,],x=y,a=tt,h=0)
    
    smatVaRES1[i,] <- sfVaRES(r1=VaRoutb[i,],r2=ESout[i,], x=y,a=nn,h=1)
    smatVaRES0[i,] <- sfVaRES(r1=VaRoutb[i,],r2=ESout[i,], x=y,a=nn,h=0)
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
  
  out = cbind(method,and,round(barVaR,3),and,round(viol,1),and,round(sbarVaR1/(1-aa),4),lbr,rank(sbarVaR1),rbr,and,round(sbarVaR0/(1-aa),4),
              lbr,rank(sbarVaR0),rbr,and,
              round(barEXP,3),and,round(sbarEXP1/(1-tt),4),lbr,rank(sbarEXP1),rbr,and,round(sbarEXP0/(1-tt),4),lbr,rank(sbarEXP0),rbr,and,
              round(barES,3),and,round(sbarVaRES1/(1-nn),4),lbr,rank(sbarVaRES1),rbr,and,round(sbarVaRES0/(1-nn),4),lbr,rank(sbarVaRES0),rbr)

   return(out)
}


# ===================================================================
# Simulation of the data generating process used in the simulation studies
# Model: AR(1)-GARCH(1,1) with skew-t innovations
# Input:
# - length of the series (default = 500+250=750)
# - model parameters (default = as in the paper setting)
## mu=-.05; ar1 = .3 # AR(1) part
## omega=.01; al=.1; be=.85 # GARCH(1,1) parameters
## Innovation distribution: skewed t with nu=5, xi=1.5
#
# Output = list of the following items:
# - simulated time series of a given size
# - list with conditional mean (mu[t]) and volatility (sig[t])
# - seed
# Notes:
# - burn-in period: 1000 observations
# ===================================================================

simDGP <- function(N=750,mu=-.05,ar1=.3,omega=.01,al=.1,be=.85,nu=5,xi=1.5,seed){

burn=1000  
n=N+burn # desired length series plus the burn-in period

# Series of innovations
set.seed(seed=seed)
x=rsstd(n=n, nu=nu, xi=xi)
# plot(density(x))

# AR(1)-GARCH(1,1) filter
r=vector(mode="logical", length=n) # process, representing negated log-returns
mu.t=vector(mode="logical", length=n)
sig.t=vector(mode="logical", length=n)

# starting values
r[1]<-0; mu.t[1]<-0; sig.t[1]<-0

for (i in 2:n)
{
  mu.t[i]=mu+ar1*r[i-1]
  sig.t[i]=sqrt(omega+al*(r[i-1]-mu.t[i-1])^2+be*sig.t[i-1]^2)
  r[i]=mu.t[i]+sig.t[i]*x[i]
}

simdat = r[-(1:burn)]
mut = mu.t[-(1:burn)]
sigt = sig.t[-(1:burn)]
# plot(simdat, type="l", ylab="")

return(list(xt=simdat,mut=mut,sigt=sigt))
}

