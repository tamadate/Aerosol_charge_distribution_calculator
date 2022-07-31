import numpy as np
import matplotlib.pyplot as plt
import functions
import physicalProp
import Fuchs
import LD



class collisionRate:
	PPs=physicalProp.PP()
	N=7
	C=np.zeros(N*2+1)
	dt=1e-8
	tmax=1

	Fuchs=Fuchs.Fuchs()
	LD=LD.LD()

	def computeFuchsBeta(self):
		self.Fuchs.computeBeta(self.PPs)

	def computeLDBeta(self):
		self.LD.computeBetaLD(self.PPs)


	## Population balence equation
	def func_pop(self,C0,C1,C_1,i):
		term1=-C0*self.PPs.Cion[0]*self.PPs.beta[i][0]
		term2=-C0*self.PPs.Cion[1]*self.PPs.beta[i][1]
		term3=C_1*self.PPs.Cion[0]*self.PPs.beta[i-1][0]
		term4=C1*self.PPs.Cion[1]*self.PPs.beta[i+1][1]
		return term1+term2+term3+term4

	def Time_develop(self):
		N=7
		K1=np.zeros(N*2+1)
		K2=np.zeros(N*2+1)
		K3=np.zeros(N*2+1)	
		K4=np.zeros(N*2+1)
		Cs=np.zeros((1,15))
		CT=0
		for i in np.arange(N*2+1):
			CT += self.C[i]
		Cs[0]=self.C/CT
		ts=np.zeros(1)
		loop=0
		for t in np.arange(0,self.tmax,self.dt):
			for i in np.arange(1,N*2):
				K1[i] = self.dt*self.func_pop(self.C[i], self.C[i+1], self.C[i-1], i)
			for i in np.arange(1,N*2):
				K2[i] = self.dt*self.func_pop(self.C[i]+0.5*K1[i], self.C[i+1]+0.5*K1[i+1], self.C[i-1]+0.5*K1[i-1], i)
			for i in np.arange(1,N*2):
				K3[i] = self.dt*self.func_pop(self.C[i]+0.5*K2[i], self.C[i+1]+0.5*K2[i+1], self.C[i-1]+0.5*K2[i-1], i)
			for i in np.arange(1,N*2):
				K4[i] = self.dt*self.func_pop(self.C[i]+K3[i], self.C[i+1]+K3[i+1], self.C[i-1]+K3[i-1], i)
			for i in np.arange(1,N*2):
				self.C[i] += (K1[i]+2*K2[i]+2*K3[i]+K4[i])/6.0
			CT=0
			for i in np.arange(N*2+1):
				CT += self.C[i]
			if loop%1000==0:
				Cs=np.append(Cs,[self.C]/CT,axis=0)
				ts=np.append(ts,t)
				plt.cla()
				for i in np.arange(N*2+1):
					if(Cs.T[i][int(loop/1000)]>1e-3):
						plt.plot(ts,Cs.T[i],label="z="+str(i-7))
				plt.xlabel("Time [s]")
				plt.ylabel("Normalized concentration [-]")
				plt.legend(loc='upper right')
				plt.pause(0.01)
			loop+=1
		np.savetxt("chargeDistribution.dat", self.C.T)

