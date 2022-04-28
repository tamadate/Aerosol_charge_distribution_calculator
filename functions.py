import math
import numpy as np
from scipy.optimize import minimize

k=1.38e-23
T=300

def Cc(dp):
	ramda=67.0e-9
	ap=dp/2.0
	return 1+ramda/ap*(1.257+0.4*np.exp(-1.1*ap/ramda))

def Zp_from_dp(dp):
	ramda=67.0e-9
	ap=dp/2.0
	Cc=1+ramda/ap*(1.257+0.4*np.exp(-1.1*ap/ramda))
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



