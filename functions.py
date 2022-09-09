import math
import numpy as np
from scipy.optimize import minimize

k=1.38e-23
T=300

def Cc(dp):
	lamda=67.0e-9
	ap=dp/2.0
	return 1+lamda/ap*(1.257+0.4*np.exp(-1.1*ap/lamda))

def Zp_from_dp(dp):
	lamda=67.0e-9
	ap=dp/2.0
	Cc=1+lamda/ap*(1.257+0.4*np.exp(-1.1*ap/lamda))
	myu=1.822e-5
	return Cc*1.6e-19/(3*math.pi*myu*dp)


def Zp_from_V(V,L,Rin,Rout,Qs):
	return Qs/V*0.5/L/math.pi*np.log(Rout/Rin)

def V_from_Zp(Zp,L,Rin,Rout,Qs):
	return Qs/Zp*0.5/L/math.pi*np.log(Rout/Rin)

def dp_from_Zp(Zp):
	def dZp2(dp):
		return ((Zp_from_dp(dp)-Zp)*1e10)**2
	result=minimize(dZp2,1e-9,method="Nelder-Mead")
	return(result.x[0])

def mp(dp,rho):
	return(np.pi/6.0*dp*dp*dp*rho)

def lamda(D,m,T):
	return(1.329*D*(m/k/T)**0.5)

def delta(ap,lamda):
	r_a=lamda/ap
	r_aSQ=r_a*r_a
	return(ap/r_aSQ*(0.2*(1+r_a)**5-1/3.0*(1+r_aSQ)*(1+r_a)**3+2/15.0*(1+r_aSQ)**2.5))

def Diff(myu,dp,T):
	return (Cc(dp)*1.38e-23*T/(3.0*np.pi*myu*dp))



