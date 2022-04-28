import math
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import optimize
import matplotlib.pylab as plt
import functions

def V_dp(Qs,Qa,L,rin,rout,Vmin,Vmax):
	T=300.0
	myu=1.8e-5
	ram=6.6e-8
	e=1.6e-19
	nd=1000.0
	V=np.arange(Vmin,Vmax,(Vmax-Vmin)/nd)
	Dp=np.arange(int(nd))*1e-9

	Qs=Qs/1000.0/60.0
	Qa=Qa/1000.0/60.0
	rin/=1000.0
	rout/=1000.0
	L/=1000.0
	beta=Qa/Qs
	delta=(Qs-Qa)/(Qs+Qa)

	Zp=functions.Zp_from_V(V,L,rin,rout,Qs)
	for i in np.arange(int(nd)):
		def dZp2(dp):
			return ((functions.Zp_from_dp(dp*1e-9)-Zp[i])*1e7)**2
		result=minimize(dZp2,6,method="Nelder-Mead")
		Dp[i]=(result.x[0])

	plt.plot(V,Dp)
	plt.yscale('log')
	plt.xscale('log')
	plt.xlabel("Voltage [V]")
	plt.ylabel("Particle diameter [nm]")
	plt.tick_params(direction='in', length=6)
	plt.title("Q$\mathsf{_{a}}$="+str(Qa*1000*60)+"L/min, Q$\mathsf{_{s}}$="+str(Qs*1000*60)+"L/min,\nL="+str(L*1000)+"mm, R$\mathsf{_{in}}$="+str(rin*1000)+"mm, R$\mathsf{_{out}}$="+str(rout*1000))
	plt.show()
	return 0



	
	
	
