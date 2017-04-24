# ===================================================
# Traditional backtesting analysis of the simulation data
# Data generating process: AR(1)-GARCH(1,1) with skewed t innovations
# ===================================================
source("Rcode_sim_RMdataprep.R")

### summary for VaR and expectiles
hac=FALSE
cctVaR = matrix(nrow=3*nm, ncol=4) # p-values for average/conditional calibration tests for VaR

for(i in 1:nm)
{
	tmp1=cct.2s.VaR(x=y, r=VaRout1[i,],lev=avec[1],hac=hac)
	tmp2=cct.1s.VaR(x=y, r=VaRout1[i,],lev=avec[1],hac=hac)
	cctVaR[i,] = c(tmp1$pv.avg, tmp1$pv.cond,tmp2$pv.avg, tmp2$pv.cond)	
	
	tmp1=cct.2s.VaR(x=y, r=VaRout2[i,],lev=avec[2],hac=hac)
	tmp2=cct.1s.VaR(x=y, r=VaRout2[i,],lev=avec[2],hac=hac)
	cctVaR[(nm+i),] = c(tmp1$pv.avg, tmp1$pv.cond,tmp2$pv.avg, tmp2$pv.cond)

	tmp1=cct.2s.VaR(x=y, r=VaRout3[i,],lev=avec[3],hac=hac)
	tmp2=cct.1s.VaR(x=y, r=VaRout3[i,],lev=avec[3],hac=hac)
	cctVaR[(2*nm+i),] = c(tmp1$pv.avg, tmp1$pv.cond,tmp2$pv.avg, tmp2$pv.cond)
}
round(cctVaR,3)
# ===================================================

cctEXP = matrix(nrow=3*nm, ncol=4) # p-values for average/conditional calibration tests for expectile

for(i in 1:nm)
{
	tmp1=cct.2s.EXP(x=y, r=EXPout1[i,],lev=tvec[1],sigt=sigt.mat[i,],hac=hac)
	tmp2=cct.1s.EXP(x=y, r=EXPout1[i,],lev=tvec[1],sigt=sigt.mat[i,],hac=hac)
	cctEXP[i,] = c(tmp1$pv.avg, tmp1$pv.cond,tmp2$pv.avg, tmp2$pv.cond)	
	
	tmp1=cct.2s.EXP(x=y, r=EXPout2[i,],lev=tvec[2],sigt=sigt.mat[i,],hac=hac)
	tmp2=cct.1s.EXP(x=y, r=EXPout2[i,],lev=tvec[2],sigt=sigt.mat[i,],hac=hac)
	cctEXP[(nm+i),] = c(tmp1$pv.avg, tmp1$pv.cond,tmp2$pv.avg, tmp2$pv.cond)

	tmp1=cct.2s.EXP(x=y, r=EXPout3[i,],lev=tvec[3],sigt=sigt.mat[i,],hac=hac)
	tmp2=cct.1s.EXP(x=y, r=EXPout3[i,],lev=tvec[3],sigt=sigt.mat[i,],hac=hac)
	cctEXP[(2*nm+i),] = c(tmp1$pv.avg, tmp1$pv.cond,tmp2$pv.avg, tmp2$pv.cond)
}
round(cctEXP,3)

# ===================================================
# Summary for (VaR,ES)

cct = matrix(nrow=3*nm, ncol=4) # p-values for conditional calibration tests for (VaR,ES) pairs

hac=FALSE

for(i in 1:nm)
{
	rr=cbind(VaRout1b[i,], ESout1[i,])
  tmp1=cct.twosided2(x=y, r=rr, lev=nvec[1],sigt=sigt.mat[i,],hac=hac)
	tmp2=cct.onesided2(x=y, r=rr, lev=nvec[1],sigt=sigt.mat[i,],hac=hac)
	cct[i,] = c(tmp1$pv.avg,tmp1$pv.cond,tmp2$pv.avg,tmp2$pv.cond)
	
	rr=cbind(VaRout2b[i,], ESout2[i,])
	tmp1=cct.twosided2(x=y, r=rr, lev=nvec[2],sigt=sigt.mat[i,],hac=hac)
	tmp2=cct.onesided2(x=y, r=rr, lev=nvec[2],sigt=sigt.mat[i,],hac=hac)
	cct[(nm+i),] = c(tmp1$pv.avg,tmp1$pv.cond,tmp2$pv.avg,tmp2$pv.cond)

	rr=cbind(VaRout3b[i,], ESout3[i,])
	tmp1=cct.twosided2(x=y, r=rr, lev=nvec[3],sigt=sigt.mat[i,],hac=hac)
	tmp2=cct.onesided2(x=y, r=rr, lev=nvec[3],sigt=sigt.mat[i,],hac=hac)
	cct[(2*nm+i),] = c(tmp1$pv.avg,tmp1$pv.cond,tmp2$pv.avg,tmp2$pv.cond)
	
}

round(cct,3)

