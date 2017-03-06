# ===================================================
# Comparative backtesting analysis of the simulation data
# Data generating process: AR(1)-GARCH(1,1) with skewed t innovations
# ===================================================
library(plotrix) #color2D.matplot

source("Rcode_sim_RMdataprep.R")
source("Rfns.R")

### ===================================================
### Traffic light matrices (method in row is the standard one)
### ===================================================

# VaR
rm=1
tlm1a <- TLMfn(r=VaRout1, y=y, rm=rm, lev=avec[1],h=1)
tlm1b <- TLMfn(r=VaRout1, y=y, rm=rm, lev=avec[1],h=0)
tlm3a <- TLMfn(r=VaRout3, y=y, rm=rm, lev=avec[3],h=1)
tlm3b <- TLMfn(r=VaRout3, y=y, rm=rm, lev=avec[3],h=0)

postscript(file="plot_sim_tlm_VaR2.eps", width=8, height=8, horizontal=FALSE)
par(mfrow=c(2,2))
plotTLM(tlm1a$TLM05, rm=rm, lev=avec[1])
plotTLM(tlm3a$TLM05, rm=rm, lev=avec[3])
plotTLM(tlm1b$TLM05, rm=rm, lev=avec[1])
plotTLM(tlm3b$TLM05, rm=rm, lev=avec[3])
dev.off()


# Expectiles
rm=2
tlm1a <- TLMfn(r=EXPout1, y=y, rm=rm, lev=tvec[1],h=1)
tlm1b <- TLMfn(r=EXPout1, y=y, rm=rm, lev=tvec[1],h=0)
tlm3a <- TLMfn(r=EXPout3, y=y, rm=rm, lev=tvec[3],h=1)
tlm3b <- TLMfn(r=EXPout3, y=y, rm=rm, lev=tvec[3],h=0)


postscript(file="plot_sim_tlm_EXP2.eps", width=8, height=8, horizontal=FALSE)
par(mfrow=c(2,2))
plotTLM(tlm1a$TLM05, rm=rm, lev=tvec[1])
plotTLM(tlm3a$TLM05, rm=rm, lev=tvec[3])
plotTLM(tlm1b$TLM05, rm=rm, lev=tvec[1])
plotTLM(tlm3b$TLM05, rm=rm, lev=tvec[3])
dev.off()


# (VaR, ES)
rm = 3
tlm1a <- TLMfn2(r=VaRout1b, r2 = ESout1, y=y, lev=nvec[1],h=1)
tlm1b <- TLMfn2(r=VaRout1b, r2 = ESout1, y=y, lev=nvec[1],h=0)
tlm3a <- TLMfn2(r=VaRout3b, r2 = ESout3, y=y, lev=nvec[3],h=1)
tlm3b <- TLMfn2(r=VaRout3b, r2 = ESout3, y=y, lev=nvec[3],h=0)

postscript(file="plot_sim_tlm_VaRES2.eps", width=8, height=8, horizontal=FALSE)
par(mfrow=c(2,2))
plotTLM(tlm1a$TLM05, rm=rm, lev=nvec[1])
plotTLM(tlm3a$TLM05, rm=rm, lev=nvec[3])
plotTLM(tlm1b$TLM05, rm=rm, lev=nvec[1])
plotTLM(tlm3b$TLM05, rm=rm, lev=nvec[3])
dev.off()


### ===================================================
### Summary table of average risk measure estimates with scores and corresponding rankings
### ===================================================

out1=rmf_comp(VaRout1, EXPout1,VaRout1b, ESout1, y, aa=avec[1], tt=tvec[1], nn=nvec[1])
out2=rmf_comp(VaRout2, EXPout2,VaRout2b, ESout2, y, aa=avec[2], tt=tvec[2], nn=nvec[2])
out3=rmf_comp(VaRout3, EXPout3,VaRout3b, ESout3, y, aa=avec[3], tt=tvec[3], nn=nvec[3])
out=rbind(out1,out2,out3)
