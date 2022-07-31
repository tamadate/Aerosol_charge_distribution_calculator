import numpy as np
import functions

class PP:
	k=1.38e-23
	Nw=6.02e23
	e=1.6e-19
	e0 = 8.85e-12
	e_fourPiE0=e/(4*np.pi*e0)
	ee_fourPiE0=e*e_fourPiE0

	T=300
	p=1e5
	rhop=1e3
	ep=1000000

	Mgas=28*0.8+32*0.2
	mgas=Mgas/Nw

	N=7
	Dp=7.08e-6
	z=np.arange(-1*N,N+1)
	beta=np.zeros((N*2+1,2))
	Cion=np.array((1e15,1e15))
	Mion=np.array((0.109,0.05))
	mion=Mion/Nw
	mig=1/(1/mion+1/mgas)
	Kion=np.array((1.4e-4,1.9e-4))
	Dion=Kion/e*k*T
	cion=(8*k*T/np.pi/mion)**0.5
	zion=np.array((1,-1))

	dp=10e-9
	ap=dp*0.5
	Cc=functions.Cc(dp)
	mp=functions.mp(dp,rhop)
	ramda=functions.ramda(Dion,mig,T)
	delta=functions.delta(ap,ramda)

	def calculatePPs(self):
		self.mgas=self.Mgas/self.Nw*1e-3
		self.mion=self.Mion/self.Nw*1e-3
		self.mp=functions.mp(self.dp,self.rhop)
		self.mion=self.mion*self.mp/(self.mion+self.mp)
		self.mig=1/(1/self.mion+1/self.mgas)
		self.Dion=self.Kion/self.e*self.k*self.T+self.Dp
		self.cion=(8*self.k*self.T/np.pi/self.mion)**0.5
		self.ap=self.dp*0.5
		self.Cc=functions.Cc(self.dp)
		self.ramda=functions.ramda(self.Dion,self.mig,self.T)
		self.delta=functions.delta(self.ap,self.ramda)
		self.bm=np.zeros((2*self.N+1,2))+1
		self.integral=np.zeros((2*self.N+1,2))


