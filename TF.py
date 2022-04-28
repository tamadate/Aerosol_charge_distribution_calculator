import math
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import optimize
from scipy import special
from scipy import integrate
import matplotlib.pylab as plt
import functions

T=300.0
kb=1.38e-23
myu=1.8e-5
ram=6.6e-8
e=1.6e-19

def Ig(w,gamma,A):
	logg=np.log(gamma)
	logw=np.log(w)
	ter1=-w**2*0.5*((1-gamma)*logw-(1-w)*logg)**2
	ter2=(w**2*0.5*(1-gamma)+w**3/3.0*logg)*((1-gamma)*logw-(1-w)*logg)
	ter3=0.25*(1-gamma)**2*(1-w**2)
	ter4=5/18.0*(1-w**3)*(1-gamma)*logg
	ter5=1/12.0*(1-w**4)*(logg)**2
	return A**2*(ter1+ter2+ter3+ter4+ter5)/(1-gamma)

def Fw(w,A,gamma,F):
	return (-A*(w*np.log(w)+(1-w)+0.5*(1-w)**2/(1-gamma)*np.log(gamma))-F)**2*1e10

def Omega(Zp, beta, sigma, delta):
	def epsilon(x):
		return x*special.erf(x)+np.exp(-x*x)/math.pi**0.5
	efunc1=epsilon((Zp-(1+beta))*0.5**0.5/sigma)
	efunc2=epsilon((Zp-(1-beta))*0.5**0.5/sigma)
	efunc3=epsilon((Zp-(1+beta*delta))*0.5**0.5/sigma)
	efunc4=epsilon((Zp-(1-beta*delta))*0.5**0.5/sigma)
	return sigma/2**0.5/beta/(1-delta)*(efunc1+efunc2-efunc3-efunc4)


def FWHM(beta,sigma,delta,peak):
	def F(x):
		return (peak*0.5-Omega(x, beta, sigma, delta))**2
	return 	(minimize(F,1,method="Nelder-Mead").x[0]-1)*2.0


def OMEGA(Qs,Qa,L,rin,rout,dp,file_name):
	f=open(file_name,'w')
	f.close
	Qs=Qs/1000.0/60.0
	Qa=Qa/1000.0/60.0
	QT=Qs+Qa
	rin/=1000.0
	rout/=1000.0
	L/=1000.0
	beta=Qa/Qs
	delta=0
	gamma=(rin/rout)**2
	A=-1/(0.5*(1+gamma)*np.log(gamma)+1-gamma)
	kappa=L/rout
	Fwa=Qa/2.0/QT
	Fws=1-Qa/2.0/QT
	wa=minimize(Fw,0.7,args=(A,gamma,Fwa),method="Nelder-Mead").x[0]
	ws=minimize(Fw,0.3,args=(A,gamma,Fws),method="Nelder-Mead").x[0]
	GDMA=(4*(1+beta)**2/(1-gamma)*(Ig(ws,gamma,A)-Ig(wa,gamma,A))+(wa-ws)/kappa**2)*-np.log(gamma)*0.5

	dp=dp*1e-9
	Zp_mean=functions.Zp_from_dp(dp)
	V=functions.V_from_Zp(Zp_mean,L,rin,rout,Qs+Qa)


	sigma=((kb*T/e/V)*GDMA)**0.5
	peak=Omega(1,beta,sigma,delta)
	fwhm=FWHM(beta,sigma,delta,peak)

	Zp_tilde=np.arange(0.5,1.5,0.001)
	omega1=Omega(Zp_tilde,beta,sigma,delta)
	omega2=Omega(Zp_tilde,beta,0.000000001,delta)
	Zp=Zp_tilde*Zp_mean
	Dp=np.arange(0.5,1.5,0.001)
	for i in np.arange(len(Dp)):
		f=open(file_name,'a')
		Dp[i]=functions.dp_from_Zp(Zp[i])
		f.write(str(Zp_tilde[i])+"\t"+str(Dp[i])+"\t"+str(Zp[i])+"\t"+str(omega1[i])+"\t"+str(omega2[i])+"\n")
		f.close
	
	plt.plot(Dp*1e9,omega1)
	plt.plot(Dp*1e9,omega2, "--")
	plt.xlabel("Particle diameter [nm]")
	plt.ylabel("Transfer function [-]")
	plt.tick_params(direction='in', length=6)
	plt.text(Dp[400]*1e9,0.85,"NFWHM={:0.3f}".format(fwhm))
	plt.text(Dp[400]*1e9,0.8,"Integral={:0.3f}".format(integrate.simps(omega1,Zp_tilde)))
	plt.show()
	return 0



def peak_FWHM(Qs,Qa,L,rin,rout,Vmin,Vmax,dp,file_name):
	f=open(file_name,'w')
	f.close
	Qs=Qs/1000.0/60.0
	Qa=Qa/1000.0/60.0
	QT=Qs+Qa
	rin/=1000.0
	rout/=1000.0
	L/=1000.0
	beta=Qa/Qs
	delta=0
	gamma=(rin/rout)**2
	A=-1/(0.5*(1+gamma)*np.log(gamma)+1-gamma)
	kappa=L/rout
	Fwa=Qa/2.0/QT
	Fws=1-Qa/2.0/QT
	wa=minimize(Fw,0.7,args=(A,gamma,Fwa),method="Nelder-Mead").x[0]
	ws=minimize(Fw,0.3,args=(A,gamma,Fws),method="Nelder-Mead").x[0]
	GDMA=(4*(1+beta)**2/(1-gamma)*(Ig(ws,gamma,A)-Ig(wa,gamma,A))+(wa-ws)/kappa**2)*-np.log(gamma)*0.5

	V=np.arange(Vmin,Vmax,(Vmax-Vmin)/1000.0)

	sigma=((kb*T/e/V)*GDMA)**0.5
	peak=Omega(1,beta,sigma,delta)
	Zp=functions.Zp_from_V(V,L,rin,rout,Qs+Qa)
	Dp=np.arange(0.5,1.5,0.001)
	fwhm=np.arange(0.5,1.5,0.001)
	for i in np.arange(len(Dp)):
		f=open(file_name,'a')
		Dp[i]=functions.dp_from_Zp(Zp[i])
		fwhm[i]=FWHM(beta,sigma[i],delta,peak[i])
		f.write(str(V[i])+"\t"+str(Dp[i])+"\t"+str(Zp[i])+"\t"+str(sigma[i])+"\t"+str(peak[i])+"\t"+str(fwhm[i])+"\n")
		f.close

	fig=plt.figure()
	ax1=fig.add_subplot(2,2,1, xlabel="Voltage [V]", ylabel="Particle diameter [nm]", )
	ax2=fig.add_subplot(2,2,2, xlabel="Particle diameter [nm]", ylabel="Peak height [-]")
	ax3=fig.add_subplot(2,2,3, xlabel="Particle diameter [nm]", ylabel="NFWHM [-]")
	ax4=fig.add_subplot(2,2,4, xlabel="Zp/Zp* [-]", ylabel="Transfer function [-]")

	ax1.plot(V,Dp*1e9)
	ax2.plot(Dp*1e9,peak)
	ax3.plot(Dp*1e9,fwhm)



	dp=dp*1e-9
	Zp_mean=functions.Zp_from_dp(dp)
	V=functions.V_from_Zp(Zp_mean,L,rin,rout,Qs+Qa)
	sigma=((kb*T/e/V)*GDMA)**0.5
	Zp=np.arange(0.5,1.5,0.001)
	Dp=np.arange(0.5,1.5,0.001)
	for i in np.arange(len(Dp)):
		Dp[i]=functions.dp_from_Zp(Zp[i]*Zp_mean)


	omega=Omega(Zp,beta,sigma,delta)
	ax4.plot(Zp,omega)

	sigma=0.00000001
	omega=Omega(Zp,beta,sigma,delta)
	ax4.plot(Zp,omega, "--")

	ax1.set_xscale('log')
	ax1.set_yscale('log')
	ax2.set_xscale('log')
	ax3.set_xscale('log')


	plt.tick_params(direction='in', length=6)
	plt.tight_layout()

	plt.show()
	return 0



	
	
	
