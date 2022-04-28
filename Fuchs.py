import Tkinter as tk
import numpy as np
import pandas as pd
import math
import matplotlib.pylab as plt
import functions

k=1.38e-23
Nw=6.02e23
e=1.6e-19
e0 = 8.85e-12


T=300
p=1e5
dt=1e-8
rhop=1e3
ep=1000000
tmax=1

Mgas=28*0.8+32*0.2
mgas=Mgas/Nw

N=7
C=np.zeros(N*2+1)
beta=np.zeros((N*2+1,2))
Dp=np.zeros((N*2+1))
z=np.arange(-1*N,N+1)

Cion=np.array((1e15,1e15))
Mion=np.array((0.109,0.05))
mion=Mion/Nw
mig=1/(1/mion+1/mgas)
Kion=np.array((1.4e-4,1.9e-4))
Dion=Kion/e*k*T
cion=(8*k*T/math.pi/mion)**0.5
zion=np.array((1,-1))

dp=10e-9
ap=dp*0.5
Cc=functions.Cc(dp)
mp=rhop*math.pi*dp**3/6.0
ramda = 1.329*Dion*(mig/k/T)**0.5
delta=ap**3.0/ramda**2.0*(0.2*(1+ramda/ap)**5-1/3.0*(1+ramda*ramda/ap/ap)*(1+ramda/ap)**3+2/15.0*(1+ramda*ramda/ap/ap)**2.5)

def calcDiff():
	myu=1.8e-5
	Cc=functions.Cc(dp)
	Dp=Cc*k*T/3.0/math.pi/myu/dp
	return Dp
	

## bm

def phi(r,a,zp,zi):
	KE=1/4.0/math.pi/e0*e*e
	return KE/r*zp*zi-KE*a*a*a*0.5/r/r/(r*r-a*a)*(ep-1.0)/(ep+1.0)

def F(x,a,zp,zi):
	KE=1/4.0/math.pi/e0*e*e
	Phi=KE*x/a*zp*zi-KE*x*x*x*x/2.0/a/(1-x*x)*(ep-1.0)/(ep+1.0)
	return np.exp(Phi/k/T)

def Fuchs(IF_flag):
	mgas=Mgas/Nw*1e-3
	mp=rhop*math.pi*dp**3/6.0
	mion=Mion/Nw*1e-3
	mig=1/(1/mion+1/mgas)
	mion=mion*mp/(mion+mp)
	Dion=Kion/e*k*T+Dp
	cion=(8*k*T/math.pi/mion)**0.5
	ap=dp*0.5
	Cc=functions.Cc(dp)
	ramda = 1.329*Dion*(mig/k/T)**0.5
	delta=ap**3.0/ramda**2.0*(0.2*(1+ramda/ap)**5-1/3.0*(1+ramda*ramda/ap/ap)*(1+ramda/ap)**3+2/15.0*(1+ramda*ramda/ap/ap)**2.5)
	print(delta)
	bm=np.zeros((2*N+1,2))+1
	integral=np.zeros((2*N+1,2))
	if(IF_flag==1):
		for i in np.arange(2):
			dr=(delta[i]-ap)/100000.0
			for rm in np.arange(ap+dr,delta[i]-dr,dr):
				dphi=1+(phi(delta[i],ap,z,zion[i])-phi(rm,ap,z,zion[i]))/(k*T*1.5)
				for j in np.arange(2*N+1):
					if(dphi[j]<0):
						dphi[j]=1000000
				b=rm*dphi**0.5
				for j in np.arange(N*2+1):
					if (b[j]<bm[j][i]):
						bm[j][i]=b[j]
						if (bm[j][i]<0.0):
							bm[j][i]=0.0
			M=100000
			h=ap/delta[i]/M
			integral.T[i]=F(0.0,ap,z,zion[i])+F(ap/delta[i],ap,z,zion[i])
			for j in np.arange(M*0.5):
				I = i
				integral.T[i]+=4*F((2*j+1)*h,ap,z,zion[i])+2*F((2*j)*h,ap,z,zion[i])
			integral.T[i]-=2*F(0.0,ap,z,zion[i])
			integral.T[i]=integral.T[i]/3.0*h

		pd=bm*bm/delta/delta
		for i in np.arange(N*2+1):
			beta[i]=math.pi*cion*pd[i]*delta*delta*np.exp(-phi(delta,ap,z[i],zion)/k/T)/(1+np.exp(-phi(delta,ap,z[i],zion)/k/T)*cion*pd[i]*delta*delta/4.0/Dion/ap*integral[i])

	if(IF_flag==0):
		for i in np.arange(2):
			for j in np.arange(N*2+1):
				if(z[j]==0):
					beta_d=math.pi*ap*ap*cion[i]
					fpdD=4*math.pi*delta[i]*Dion[i]
					beta[j][i]=fpdD/(1+fpdD/beta_d)
				else:
					gamma=1-zion[i]*z[j]*e*e/4.0/np.pi/e0/k/T*(1/ap-1/delta[i])
					if gamma>0:
						beta_d=math.pi*ap*ap*cion[i]*gamma
						PHI=-zion[i]*z[j]*e*e/4.0/np.pi/e0/k/T/delta[i]
						first_term=z[j]*zion[i]*Dion[i]*e*e/e0/k/T/(np.exp(-PHI)-1)
						second_term=1-Dion[i]*z[j]*zion[i]*e*e/beta_d/e0/k/T/(np.exp(PHI)-1)
						beta[j][i]=first_term/second_term
					else:
						beta[j][i]=1e-200


def CG():
	mgas=Mgas/Nw*1e-3
	mp=rhop*math.pi*dp**3/6.0
	mion=Mion/Nw*1e-3
	mig=1/(1/mion+1/mgas)
	mion=mion*mp/(mion+mp)
	Dion=Kion/e*k*T+Dp
	print(Dion)
	cion=(8*k*T/math.pi/mion)**0.5
	ap=dp*0.5
	Cc=functions.Cc(dp)
	ramda = 1.329*Dion*(mig/k/T)**0.5
	for i in np.arange(2):
		for j in np.arange(N*2+1):
			PsiE=-zion[i]*z[j]*e*e/4.0/math.pi/e0/k/T/ap
			Psid=-zion[i]*z[j]*e*e/4.0/math.pi/e0/k/T/delta[i]
			if(PsiE==0):
				itaC=1
			else:
				itaC=PsiE/(1-np.exp(-PsiE))
			if(PsiE>0):
				itaFM=1+PsiE
			else:
				itaFM=np.exp(PsiE)
			friction=k*T/Dion[i]
			KnD=(mion[i]*k*T)**0.5*itaC/friction/ap/itaFM
			H=(4*math.pi*KnD**2+25.836*KnD**3+(8*math.pi)**0.5*11.211*KnD**4)/(1+3.502*KnD+7.211*KnD**2+11.211*KnD**3)
			if(PsiE>0):
				myuA=2.5
				myuB=4.528*np.exp(-1.088*PsiE)+0.7091*np.log(1+1.537*PsiE)
				myuC=11.36*PsiE**0.272-10.33
				myuk=-0.003533*PsiE+0.05971
				inmyu=1+myuk*(np.log(KnD)-myuB)/myuA
				Myu=myuC/myuA*(inmyu)**(-1/myuk-1)*np.exp(-((inmyu)**(-1/myuk)))
				H=np.exp(Myu)*H
			print(zion[i])
			print(z[j])
			print friction
			beta[j][i]=H*friction*ap**3*itaFM**2/mion[i]/itaC
		
		

def electron():
	Mion=5.486e-4
	mion=Mion/Nw*1e-3
	mion=mion*mp/(mion+mp)
	Te=2*e/k/1.5
	cion=(8*k*Te/math.pi/mion)**0.5
	ap=dp*0.5-1.5e-10
	i=1
	for j in np.arange(N*2+1):
		PsiE=-zion[i]*z[j]*e*e/4.0/math.pi/e0/k/Te/ap
		if(PsiE>0):
			itaFM=1+PsiE
		else:
			itaFM=np.exp(PsiE)
		beta[j][i]=cion*math.pi*ap*ap*itaFM


	

def func_pop(C0,C1,C_1,i):
	return -C0*Cion[0]*beta[i][0]-C0*Cion[1]*beta[i][1]+C_1*Cion[0]*beta[i-1][0]+C1*Cion[1]*beta[i+1][1]


def Time_develop():
	K1=np.zeros(N*2+1)
	K2=np.zeros(N*2+1)
	K3=np.zeros(N*2+1)	
	K4=np.zeros(N*2+1)
	Cs=np.zeros((1,15))
	CT=0
	for i in np.arange(N*2+1):
		CT += C[i]
	Cs[0]=C/CT
	ts=np.zeros(1)
	loop=0
	for t in np.arange(0,tmax,dt):
		for i in np.arange(1,N*2):
			K1[i] = dt*func_pop(C[i], C[i+1], C[i-1], i)
		for i in np.arange(1,N*2):
			K2[i] = dt*func_pop(C[i]+0.5*K1[i], C[i+1]+0.5*K1[i+1], C[i-1]+0.5*K1[i-1], i)
		for i in np.arange(1,N*2):
			K3[i] = dt*func_pop(C[i]+0.5*K2[i], C[i+1]+0.5*K2[i+1], C[i-1]+0.5*K2[i-1], i)
		for i in np.arange(1,N*2):
			K4[i] = dt*func_pop(C[i]+K3[i], C[i+1]+K3[i+1], C[i-1]+K3[i-1], i)
		for i in np.arange(1,N*2):
			C[i] += (K1[i]+2*K2[i]+2*K3[i]+K4[i])/6.0
		CT=0
		for i in np.arange(N*2+1):
			CT += C[i]
		if loop%1000==0:
			Cs=np.append(Cs,[C]/CT,axis=0)
			ts=np.append(ts,t)
			plt.cla()
			plt.plot(ts,Cs.T[7],label="z=0")
			plt.plot(ts,Cs.T[6],label="z=-1")
			plt.plot(ts,Cs.T[8],label="z=+1")
			plt.plot(ts,Cs.T[5],label="z=-2")
			plt.plot(ts,Cs.T[9],label="z=+2")
			plt.xlabel("Time [s]")
			plt.ylabel("Normalized concentration [-]")
			plt.legend(loc='upper right')
			plt.pause(0.01)
		loop+=1
	np.savetxt("chargeDistribution.dat", C.T)

		
