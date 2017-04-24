## Investigation of the impact of a small out-of-sample size on the results
## of forecasting performance comparisons based on consistent scoring rules

install.packages("rugarch")
install.packages("fGarch")
install.packages("MASS")
install.packages("expectreg")
install.packages("ismev")
install.packages("lmom")
install.packages("QRM")
install.packages("skewt")
install.packages("reshape2")
install.packages("plotrix")


library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library("reshape2")
library("plotrix")

source("Rfns.R")

method=c("n-FP  ", "n-FHS ","n-EVT ","t-FP  ","t-FHS ","t-EVT ","st-FP ","st-FHS","st-EVT","opt   ")

## =======================================================================
## Simulation study
w=500 # estimation window
n=250 # out-of-sample size for assessment and testing

n=1000

nu=5; xi=1.5 #parameters of skewed t distribution of innovations
# mean and standard deviation for expectile under skew t for optimal forecast computations
mst = mean.st(shape=nu,skew=xi)
sst = sqrt(var.st(shape=nu,skew=xi))



## =======================================================================
## Confidence levels for risk measures
avec=c(.90,.99) # vector of alpha levels for VaR
tvec=c(.96561,.99855) # vector of tau levels for expectile
nvec=c(.754,.975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:2) and in pair with ES (3:4)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR
K=length(avec) # number of confidence levels considered 
## =======================================================================
## matrices to store forecasts and some parameters (temporarily within the loop)
## .n, .t and .st correspond to likelihood used to fit the AR-GARCH filter
VaR.n <- VaR.t <-VaR.st <-matrix(nrow=n,ncol=length(VaR.levels))
VaRfhs.n <- VaRfhs.t <- VaRfhs.st <- matrix(nrow=n,ncol=length(VaR.levels))
VaRevt.n <-VaRevt.t <-VaRevt.st <- matrix(nrow=n,ncol=length(VaR.levels))

ES.n <- ES.t <-ES.st <-matrix(nrow=n,ncol=length(nvec))
ESfhs.n <-ESfhs.t <-ESfhs.st <- matrix(nrow=n,ncol=length(nvec))
ESevt.n <- ESevt.t <- ESevt.st <- matrix(nrow=n,ncol=length(nvec))

EXP.n <- EXP.t <- EXP.st <- matrix(nrow=n,ncol=length(tvec))
EXPfhs.n <- EXPfhs.t <- EXPfhs.st <- matrix(nrow=n,ncol=length(tvec))
EXPevt.n <- EXPevt.t <- EXPevt.st <- matrix(nrow=n,ncol=length(tvec))

mut.n  <- mut.t  <- mut.st  <- vector(mode="numeric", length=n)
sigt.n  <- sigt.t  <- sigt.st  <- vector(mode="numeric", length=n)

flag.t <- 0 #flag to see if nu.hat <=2

## =======================================================================

B=1000 # number of samples to simulate and analyse

## =======================================================================
## Output matrices to scote average score values for each simulated sample
## Columns correspond to methods
sbarVaR1a <- sbarVaR1b <- NULL ## scoring function (h=1), confidence levels a & b (1 & 2)
sbarVaR0a <- sbarVaR0b <- NULL ## scoring function (h=0), confidence levels a & b (1 & 2)

sbarEXP1a <- sbarEXP1b <- NULL ## scoring function (h=1), confidence levels a & b (1 & 2)
sbarEXP0a <- sbarEXP0b <- NULL ## scoring function (h=0), confidence levels a & b (1 & 2)

sVaRESmat1a <- sVaRESmat1b <- matrix(nrow=n,ncol=10) ## scoring function (h=1), confidence levels a & b (1 & 2)
sVaRESmat0a <- sVaRESmat0b <- matrix(nrow=n,ncol=10) ## scoring function (h=0), confidence levels a & b (1 & 2)

sbarVaRES1a <- sbarVaRES1b <- NULL ## scoring function (h=1), confidence levels a & b (1 & 2)
sbarVaRES0a <- sbarVaRES0b <- NULL ## scoring function (h=0), confidence levels a & b (1 & 2)

viol.a <- viol.b <-matrix(nrow=B,ncol=10)

## =======================================================================

for (j in 1:B){
  
  simdat=simDGP(N=n+w, seed=j)
  x=simdat$xt
  y=tail(x,n) # verifying observations
  ymat=matrix(y, nrow=n,ncol=10) # matrix of observations for computing % violations of VaR forecasts 
  mut=tail(simdat$mut,n)
  sigt=tail(simdat$sigt,n)
  # seeds[[j]] <- simdat$seed
  
  # =======================================================
  # Forecasting assuming normal innovations
  # =======================================================
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
  qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1)
  esmodel = NULL # expected shortfall for the assumed model/distribution
  for(i in inu)
    esmodel = c(esmodel, integrate(function(x) x*dnorm(x), qmodel[i], Inf)$value/(1-VaR.levels[i]))
  emodel = enorm(asy=tvec)
  
  for(i in 1:n)
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i)
    
    VaR.n[i,]=foc$VaRmodel
    VaRfhs.n[i,] = foc$VaRfhs
    VaRevt.n[i,] = foc$VaRevt
    
    EXP.n[i,]=foc$EXPmodel
    EXPfhs.n[i,] = foc$EXPfhs
    EXPevt.n[i,] = foc$EXPevt
    
    ES.n[i,]=foc$ESmodel
    ESfhs.n[i,] = foc$ESfhs
    ESevt.n[i,] = foc$ESevt
    mut.n[i] = foc$mut
    sigt.n[i] = foc$sigt
    }
  
  # =======================================================
  # Forecasting assuming t innovations
  # =======================================================
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="std")
  
  for(i in 1:n)
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    
    nu.hat = coef(fit)["shape"] #shape parameter, or degrees of freedom of t distribution
    qmodel = qdist("std",p=VaR.levels,mu=0,sigma=1,shape=nu.hat)
    esmodel = NULL # expected shortfall for the assumed model/distribution
    for(k in inu)
      esmodel = c(esmodel, integrate(function(x) x*dt(x, df=nu.hat), qmodel[k], Inf)$value/(1-VaR.levels[k]))
    esmodel=sqrt((nu.hat-2)/nu.hat)*esmodel # assuming nu>2
    emodel = sqrt((nu.hat-2)/nu.hat)*et(asy=tvec, df=nu.hat) # assuming nu>2
    if(nu.hat<=2) flag.t=flag.t+1
    
    foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i^2)
    
    VaR.t[i,]=foc$VaRmodel
    VaRfhs.t[i,] = foc$VaRfhs
    VaRevt.t[i,] = foc$VaRevt
    
    EXP.t[i,]=foc$EXPmodel
    EXPfhs.t[i,] = foc$EXPfhs
    EXPevt.t[i,] = foc$EXPevt
    
    ES.t[i,]=foc$ESmodel
    ESfhs.t[i,] = foc$ESfhs
    ESevt.t[i,] = foc$ESevt
    
    mut.t[i] = foc$mut
    sigt.t[i] = foc$sigt
  }
  
  # =======================================================
  # Forecasting assuming skew-t innovations
  # =======================================================
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")
  
  for(i in 1:n)
    
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    
    # mean and std dev of skewed t rv
    nu.hat=coef(fit)["shape"]; ga=coef(fit)["skew"]
    m = mean.st(shape=nu.hat,skew=ga)
    s = sqrt(var.st(shape=nu.hat,skew=ga))
    
    qmodel = qskt(p=VaR.levels,df=nu.hat,gamma=ga) # for skew t variable in standard form
    esmodel = NULL # expected shortfall for the assumed model/distribution
    for(k in inu)
      esmodel = c(esmodel, integrate(function(x) x*dskt(x, df=nu.hat, gamma=ga), qmodel[k], Inf)$value/(1-VaR.levels[k]))
    esmodel=(esmodel-m)/s
    qmodel = (qmodel-m)/s # adjustment to have mean zero, variance one
    emodel = (est(asy=tvec,shape=nu.hat,skew=ga)-m)/s
    
    foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i^3)
    
    VaR.st[i,]=foc$VaRmodel
    VaRfhs.st[i,] = foc$VaRfhs
    VaRevt.st[i,] = foc$VaRevt
    
    EXP.st[i,]=foc$EXPmodel
    EXPfhs.st[i,] = foc$EXPfhs
    EXPevt.st[i,] = foc$EXPevt
    
    ES.st[i,]=foc$ESmodel
    ESfhs.st[i,] = foc$ESfhs
    ESevt.st[i,] = foc$ESevt
    
    mut.st[i] = foc$mut
    sigt.st[i] = foc$sigt
  }
  
  ## sort data into a matrix with columns corresponding to methods 
  ## and then by risk measure confidence levels
  VaRmat<-NULL; EXPmat<-NULL; VaR2mat<-NULL; ESmat<-NULL
  for (k in 1:K)
  {
    VaRa=qsstd(p=avec[k],nu=nu,xi=xi)  # a-VaR for the true innovation distribution
    VaRopt=mut + sigt*VaRa
    VaRmat = cbind(VaRmat,VaR.n[,k],VaRfhs.n[,k],VaRevt.n[,k],VaR.t[,k],VaRfhs.t[,k],VaRevt.t[,k],VaR.st[,k],VaRfhs.st[,k],VaRevt.st[,k],VaRopt)
    
    #expectile
    EXPa = (est(asy=tvec[k],shape=nu,skew=xi)-mst)/sst # a-expectile for the innovation distribution
    EXPopt=mut + sigt*EXPa
    EXPmat = cbind(EXPmat,EXP.n[,k],EXPfhs.n[,k],EXPevt.n[,k],
                   EXP.t[,k],EXPfhs.t[,k],EXPevt.t[,k],
                   EXP.st[,k],EXPfhs.st[,k],EXPevt.st[,k],EXPopt)
    
    #ES
    VaRa = qskt(p=nvec[k],df=nu,gamma=xi) # for skew t variable in standard form of the density
    ESa=integrate(function(x) x*dskt(x, df=nu, gamma=xi), VaRa, Inf)$value/(1-nvec[k])
    
    VaRa = (VaRa-mst)/sst
    ESa=(ESa-mst)/sst # adjustment for skew t distribution with mean zero and sd=1
    VaRopt=mut + sigt*VaRa
    ESopt = mut + sigt*ESa
    VaR2mat = cbind(VaR2mat,VaR.n[,(K+k)],VaRfhs.n[,(K+k)],VaRevt.n[,(K+k)],
                    VaR.t[,(K+k)],VaRfhs.t[,(K+k)],VaRevt.t[,(K+k)],
                    VaR.st[,(K+k)],VaRfhs.st[,(K+k)],VaRevt.st[,(K+k)],VaRopt)
    ESmat = cbind(ESmat,ES.n[,k],ESfhs.n[,k],ESevt.n[,k],
                ES.t[,k],ESfhs.t[,k],ESevt.t[,k],
                ES.st[,k],ESfhs.st[,k],ESevt.st[,k],ESopt)
  }

      sVaRmat1a <- apply(VaRmat[,1:10],2,sfVaR,x=y,a=avec[1],h=1)  
      sVaRmat1b <- apply(VaRmat[,11:20],2,sfVaR,x=y,a=avec[2],h=1)  
      sVaRmat0a <- apply(VaRmat[,1:10],2,sfVaR,x=y,a=avec[1],h=0)  
      sVaRmat0b <- apply(VaRmat[,11:20],2,sfVaR,x=y,a=avec[2],h=0)  
      
      sbarVaR1a <- rbind(sbarVaR1a, apply(sVaRmat1a,2,mean))
      sbarVaR1b <- rbind(sbarVaR1b, apply(sVaRmat1b,2,mean))
      sbarVaR0a <- rbind(sbarVaR0a, apply(sVaRmat0a,2,mean))
      sbarVaR0b <- rbind(sbarVaR0b, apply(sVaRmat0b,2,mean))
      
      viol.a[j,] <- apply((VaRmat[,1:10]<ymat),2,sum)/n
      viol.b[j,] <- apply((VaRmat[,11:20]<ymat),2,sum)/n
      
      sEXPmat1a <- apply(EXPmat[,1:10],2,sfEXP,x=y,a=tvec[1],h=1)  
      sEXPmat1b <- apply(EXPmat[,11:20],2,sfEXP,x=y,a=tvec[2],h=1)  
      sEXPmat0a <- apply(EXPmat[,1:10],2,sfEXP,x=y,a=tvec[1],h=0)  
      sEXPmat0b <- apply(EXPmat[,11:20],2,sfEXP,x=y,a=tvec[2],h=0)  
      
      sbarEXP1a <- rbind(sbarEXP1a, apply(sEXPmat1a,2,mean))
      sbarEXP1b <- rbind(sbarEXP1b, apply(sEXPmat1b,2,mean))
      sbarEXP0a <- rbind(sbarEXP0a, apply(sEXPmat0a,2,mean))
      sbarEXP0b <- rbind(sbarEXP0b, apply(sEXPmat0b,2,mean))
      
      for (l in 1:10) ## loop over methods
      {
        sVaRESmat1a[,l] <- sfVaRES(r1=VaR2mat[,l],r2=ESmat[,l],x=y,a=nvec[1],h=1)
        sVaRESmat0a[,l] <- sfVaRES(r1=VaR2mat[,l],r2=ESmat[,l],x=y,a=nvec[1],h=0)
        sVaRESmat1b[,l] <- sfVaRES(r1=VaR2mat[,(10+l)],r2=ESmat[,(10+l)],x=y,a=nvec[2],h=1)
        sVaRESmat0b[,l] <- sfVaRES(r1=VaR2mat[,(10+l)],r2=ESmat[,(10+l)],x=y,a=nvec[2],h=0)
      }
      
      sbarVaRES1a <- rbind(sbarVaRES1a, apply(sVaRESmat1a,2,mean))
      sbarVaRES1b <- rbind(sbarVaRES1b, apply(sVaRESmat1b,2,mean))
      sbarVaRES0a <- rbind(sbarVaRES0a, apply(sVaRESmat0a,2,mean))
      sbarVaRES0b <- rbind(sbarVaRES0b, apply(sVaRESmat0b,2,mean))
      
      
  ##if(j %% 5 == 0) print(j)
      print(j)
  
}

##======================================================================
# saving the output

save(sbarVaR1a,sbarVaR1b,sbarVaR0a,sbarVaR0b, file="sim1000_sbarVaR.RDATA")
save(sbarEXP1a,sbarEXP1b,sbarEXP0a,sbarEXP0b, file="sim1000_sbarEXP.RDATA")
save(sbarVaRES1a,sbarVaRES1b,sbarVaRES0a,sbarVaRES0b, file="sim1000_sbarVaRES.RDATA")
save(viol.a,viol.b,file="sim1000_viol.RDATA")

##======================================================================
## Analysis of the output 

colnames(sbarVaR1.a) <- colnames(sbarVaR1.b) <-method
colnames(sbarVaR0.a) <- colnames(sbarVaR0.b) <- method
colnames(sbarVaRES1.a) <- colnames(sbarVaRES1.b) <-method
colnames(sbarVaRES0.a) <- colnames(sbarVaRES0.b) <- method
colnames(viol.a) <- colnames(viol.b) <- method

## boxplots of RANKS based on average scores for VaR and (VaR, ES) at standard Basel levels

rankVaR1a <- t(apply(sbarVaR1.a,1,"rank"))
rankVaR1b <- t(apply(sbarVaR1.b,1,"rank"))

rankVaR0a <- t(apply(sbarVaR0.a,1,"rank"))
rankVaR0b <- t(apply(sbarVaR0.b,1,"rank"))

rankVaRES1a <- t(apply(sbarVaRES1.a,1,"rank"))
rankVaRES1b <- t(apply(sbarVaRES1.b,1,"rank"))

rankVaRES0a <- t(apply(sbarVaRES0.a,1,"rank"))
rankVaRES0b <- t(apply(sbarVaRES0.b,1,"rank"))

colnames(rankVaR1a) <- colnames(rankVaR1b) <-method
colnames(rankVaR0a) <- colnames(rankVaR0b) <- method
colnames(rankVaRES1a) <- colnames(rankVaRES1b) <-method
colnames(rankVaRES0a) <- colnames(rankVaRES0b) <- method
colnames(viol.a) <- colnames(viol.b) <- method


## Merging data for two scores
score = c(rep('A',dim(rankVaR1a)[1]),rep('B',dim(rankVaR1a)[1]))

datVaRa = cbind(rankVaR1a,rankVaR0a)[,c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)]
datVaRb = cbind(rankVaR1b,rankVaR0b)[,c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)]
datviol = cbind(viol.a,viol.b)

datVaRESa = cbind(rankVaRES1a,rankVaRES0a)[,c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)]
datVaRESb = cbind(rankVaRES1b,rankVaRES0b)[,c(1,11,2,12,3,13,4,14,5,15,6,16,7,17,8,18,9,19,10,20)]

atvec=(1:29)[-c(3,6,9,12,15,18,21,24,27)]
bpcol=rep(c("white","gray"),10)

postscript(file="plot_sim250VaR.eps", width=14, height=8, horizontal=FALSE)

m <- rbind(c(1, 1), c(2,2), c(3, 4))
layout(m)

boxplot.matrix(datVaRa,main=expression('Ranks for VaR'[0.90]), cex.axis=0.7,col=bpcol, at=atvec)
boxplot.matrix(datVaRb,main=expression('Ranks for VaR'[0.99]), cex.axis=0.7,col=bpcol, at=atvec)

par(mar=c(3,2,3,1)+0.1)

boxplot.matrix(viol.a, main=expression('%-violations of VaR'[0.90]), cex.axis=0.6)
abline(h=.1, col="darkgray", lwd=2, lty=2)

boxplot.matrix(viol.b, main=expression('%-violations of VaR'[0.99]), cex.axis=0.6)
abline(h=.01, col="darkgray", lwd=2, lty=2)
dev.off()


postscript(file="plot_sim250VaRES_sbs.eps", width=12, height=8, horizontal=FALSE)
par(mfrow=c(2,1))
boxplot.matrix(datVaRESa,main=expression(paste("Ranks for (",VaR[0.754],", ",ES[0.754],")")), cex.axis=0.7,col=bpcol, at=atvec)
boxplot.matrix(datVaRESb,main=expression(paste("Ranks for (",VaR[0.975],", ",ES[0.975],")")), cex.axis=0.7,col=bpcol, at=atvec)
dev.off()


## VaR ranks

boxplot.matrix(rankVaR1a,main=expression(paste(alpha," = ",0.90)), cex.axis=0.6)


postscript(file="plot_sim250VaR.eps", width=14, height=8, horizontal=FALSE)
par(mfrow=c(3,2),mar=c(3,2,3,1)+0.1)

boxplot.matrix(rankVaR1a,main=expression(paste(alpha," = ",0.90)), cex.axis=0.6)
boxplot.matrix(rankVaR0a,main=expression(paste(alpha," = ",0.90)), cex.axis=0.6)

boxplot.matrix(rankVaR1b,main=expression(paste(alpha," = ",0.99)), cex.axis=0.6)
boxplot.matrix(rankVaR0b,main=expression(paste(alpha," = ",0.99)), cex.axis=0.6)

boxplot.matrix(viol.a, main=expression('%-violations of VaR'[0.90]), cex.axis=0.6)
abline(h=.1, col="darkgray", lwd=2, lty=2)

boxplot.matrix(viol.b, main=expression('%-violations of VaR'[0.99]), cex.axis=0.6)
abline(h=.01, col="darkgray", lwd=2, lty=2)

dev.off()

## (VaR,ES) ranks
postscript(file="plot_sim250VaRES.eps", width=14, height=5.5, horizontal=FALSE)
par(mfrow=c(2,2),mar=c(3,2,3,1)+0.1, cex.main=0.9)

boxplot.matrix(rankVaRES1a,main=expression(paste(nu," = ",0.754)), cex.axis=0.5)
boxplot.matrix(rankVaRES0a,main=expression(paste(nu," = ",0.754)), cex.axis=0.5)

boxplot.matrix(rankVaRES1b,main=expression(paste(nu," = ",0.975)), cex.axis=0.5)
boxplot.matrix(rankVaRES0b,main=expression(paste(nu," = ",0.975)), cex.axis=0.5)

dev.off()

##======================================================================
## Assessing variability in ranks across samples and methods

## ranks based on the large-sample-size simulation study
rankVaR1a.star <- c(9,8,7,10,6,5,2,4,3,1)
rankVaR0a.star <- c(7,9,8,10,5,6,2,3,4,1)
rankVaR1b.star <- c(10,8,5,9,7,4,3,6,2,1)
rankVaR0b.star <- c(10,7,5,9,8,4,3,6,2,1)

rankVaRES1a.star <- c(10,5,9,7,4,8,2,3,6,1)
rankVaRES0a.star <- c(9,4,8,10,5,7,3,2,6,1)
rankVaRES1b.star <- c(10,8,7,9,5,6,2,3,4,1)
rankVaRES0b.star <- c(10,8,7,9,6,5,2,4,3,1)

B=dim(rankVaR1a)[1] # number of samples generated

## Combined statistic
sum.abs <- function(x) sum(abs(x))

se.VaR1a = apply(rankVaR1a - matrix(rep(rankVaR1a.star,B),byrow=TRUE, nrow=B),1,"sum.abs")
se.VaR0a = apply(rankVaR0a - matrix(rep(rankVaR0a.star,B),byrow=TRUE, nrow=B),1,"sum.abs")
se.VaR1b = apply(rankVaR1b - matrix(rep(rankVaR1b.star,B),byrow=TRUE, nrow=B),1,"sum.abs")
se.VaR0b = apply(rankVaR0b - matrix(rep(rankVaR0b.star,B),byrow=TRUE, nrow=B),1,"sum.abs")

se.VaRES1a = apply(rankVaRES1a - matrix(rep(rankVaRES1a.star,B),byrow=TRUE, nrow=B),1,"sum.abs")
se.VaRES0a = apply(rankVaRES0a - matrix(rep(rankVaRES0a.star,B),byrow=TRUE, nrow=B),1,"sum.abs")
se.VaRES1b = apply(rankVaRES1b - matrix(rep(rankVaRES1b.star,B),byrow=TRUE, nrow=B),1,"sum.abs")
se.VaRES0b = apply(rankVaRES0b - matrix(rep(rankVaRES0b.star,1000),byrow=TRUE, nrow=1000),1,"sum.abs")


par(mfrow=c(2,2),mar=c(4,4,3,1)+0.1)
hist(se.VaR1a, main=expression(paste(alpha," = ",0.90,", 1-homog. score")))
hist(se.VaR0a, main=expression(paste(alpha," = ",0.90,", 0-homog. score")))
hist(se.VaR1b, main=expression(paste(alpha," = ",0.99,", 1-homog. score")))
hist(se.VaR0b, main=expression(paste(alpha," = ",0.99,", 0-homog. score")))

par(mfrow=c(2,2),mar=c(4,4,3,1)+0.1)
hist(se.VaRES1a, main=expression(paste(alpha," = ",0.90,", 1-homog. score")))
hist(se.VaRES0a, main=expression(paste(alpha," = ",0.90,", 0-homog. score")))
hist(se.VaRES1b, main=expression(paste(alpha," = ",0.99,", 1-homog. score")))
hist(se.VaRES0b, main=expression(paste(alpha," = ",0.99,", 0-homog. score")))


#### Just looking at rank differences
diff.VaR1a <- rankVaR1a - matrix(rep(rankVaR1a.star,B),byrow=TRUE, nrow=B)
diff.VaR0a <- rankVaR0a - matrix(rep(rankVaR0a.star,B),byrow=TRUE, nrow=B)
diff.VaR1b <- rankVaR1b - matrix(rep(rankVaR1b.star,B),byrow=TRUE, nrow=B)
diff.VaR0b <- rankVaR0b - matrix(rep(rankVaR0b.star,B),byrow=TRUE, nrow=B)

par(mfrow=c(4,3)) 
for(i in 1:10) hist(diff.VaR1a[,i], main=method[i])

par(mfrow=c(4,3)) 
for(i in 1:10) hist(diff.VaR0a[,i], main=method[i])

par(mfrow=c(4,3)) 
for(i in 1:10) hist(diff.VaR1b[,i], main=method[i])

par(mfrow=c(4,3)) 
for(i in 1:10) hist(diff.VaR0b[,i], main=method[i])


##======================================================================
## Analyse a "strange" case
j=12 #re-run above; consider alpha=0.99 nu=0.975


save(sbarVaR1a,sbarVaR1b,sbarVaR0a,sbarVaR0b, file="sim250j12_sbarVaR.RDATA")
save(sbarEXP1a,sbarEXP1b,sbarEXP0a,sbarEXP0b, file="sim250j12_sbarEXP.RDATA")
save(sbarVaRES1a,sbarVaRES1b,sbarVaRES0a,sbarVaRES0b, file="sim250j12_sbarVaRES.RDATA")
save(viol.a,viol.b,file="sim250j12_viol.RDATA")

load("sim250j12_sbarVaR.RDATA"); load("sim250j12_sbarEXP.RDATA"); load("sim250j12_sbarVaRES.RDATA")
load("sim250j12_viol.RDATA")

# traffic light matrices

rm=1
tlm1VaR <- TLMfn(r=t(VaRmat[,11:20]), y=y, rm=rm, lev=avec[2],h=1)
tlm0VaR <- TLMfn(r=t(VaRmat[,11:20]), y=y, rm=rm, lev=avec[2],h=0)

rm=3
tlm1ES <- TLMfn2(r=t(VaR2mat[,11:20]), r2 = t(ESmat[,11:20]), y=y, lev=nvec[2],h=1)
tlm0ES <- TLMfn2(r=t(VaR2mat[,11:20]), r2 = t(ESmat[,11:20]), y=y, lev=nvec[2],h=0)


## Plots to illustrate the case
postscript(file="plot_sim250j12.eps", width=12, height=8, horizontal=FALSE)

m <- rbind(c(1, 1), c(2,3))
layout(m)

plot(x,type="l", ylab="")
lines(head(x,500), col="grey")
abline(v=500,lty=3, lwd=2)
lines(501:750,VaRout[1,], col=gray(0.7), lty=2, lwd=1.5)
lines(501:750,VaRout[10,], col=gray(0.3), lty=3, lwd=1.5)
legend("topright",legend=c(method[1],method[10]), lty=c(2,3),lwd=1.5,col=c(gray(0.7),gray(0.3)), cex=0.9)

plotTLM(tlm0VaR$TLM05, rm=1, lev=avec[2])
plotTLM(tlm0ES$TLM05, rm=3, lev=nvec[2])

dev.off()

# traditional backtesting
## VaR
nm=length(method)
VaRout=t(VaRmat[,11:20])
hac=FALSE
cctVaR = matrix(nrow=nm, ncol=4) # p-values for average/conditional calibration tests for VaR

for(i in 1:nm)
{
  tmp1=cct.2s.VaR(x=y, r=VaRout[i,],lev=avec[2],hac=hac)
  cctVaR[i,1:2] = c(tmp1$pv.avg, tmp1$pv.cond)
  
  tmp2=cct.1s.VaR(x=y, r=VaRout[i,],lev=avec[2],hac=hac)
  cctVaR[i,3:4] = c(tmp2$pv.avg, tmp2$pv.cond)	
}
(i=i+1)
round(cctVaR,3)

## output for LaTex table
lbr = rep("(",10); rbr = rep(")",10)
and = rep("&",10); and3 = rep(and,3)
nm=10 # number of methods considered

outVaR = cbind(method,and,round(apply(VaRout,1,"mean"),3),and,round(viol.b*100,2),and,
               round(sbarVaR1b,4)[1,],lbr,rank(sbarVaR1b),rbr,and,
               round(sbarVaR0b,4)[1,],lbr,rank(sbarVaR0b),rbr)
for(i in 1:4) outVaR=cbind(outVaR,and,round(cctVaR[,i],3))

