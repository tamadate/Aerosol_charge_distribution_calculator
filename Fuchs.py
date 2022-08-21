import numpy as np
import math
import functions
import physicalProp
from scipy.integrate import quad

class Fuchs():
	def phi(self,r,zp,zi,PPs):
		rinv=1/r
		return PPs.ee_fourPiE0*(rinv*zp*zi-PPs.ap*PPs.ap*PPs.ap*0.5*rinv*rinv/(r*r-PPs.ap*PPs.ap)*(PPs.ep-1.0)/(PPs.ep+1.0))

	def F(self,x,zp,zi,PPs):
		Phi=self.phi(1/x,zp,zi,PPs)
		return np.exp(Phi/PPs.k/PPs.T)


	def computeBeta(self,PPs):
		PPs.calculatePPs()
		bm=np.zeros((15,2))+1
		integral=np.zeros((15,2))
		beta=np.zeros((15,2))
		for i in np.arange(2):
			dr=PPs.delta[i]-PPs.ap
			for j in np.arange(15):			
				direction=-1
				bm[j][i]=bb=rm=PPs.delta[i]
				for loop in np.arange(1,100):
					rm+=dr*direction*(0.9**loop)
					if(rm>PPs.delta[i]):
						rm=PPs.delta[i]
					dphi=1+(self.phi(PPs.delta[i],PPs.z[j],PPs.zion[i],PPs)-self.phi(rm,PPs.z[j],PPs.zion[i],PPs))/(PPs.k*PPs.T*1.5)
					if(dphi<0):
						direction*=-1
						continue
					b=rm*dphi**0.5
					if (b<bm[j][i]):
						bm[j][i]=b
					if (bb<b):
						direction*=-1
					bb=b
				integral[j][i],err=quad(self.F,0,PPs.ap/PPs.delta[i],args=(PPs.z[j],PPs.zion[i],PPs))
		pd=bm*bm/PPs.delta/PPs.delta
		for i in np.arange(15):
			phid_kT=self.phi(PPs.delta,PPs.z[i],PPs.zion,PPs)/PPs.k/PPs.T
			cpdd=PPs.cion*pd[i]*PPs.delta*PPs.delta
			PPs.beta[i]=math.pi*cpdd*np.exp(-phid_kT)/(1+np.exp(-phid_kT)*cpdd/4.0/PPs.Dion/PPs.ap*integral[i])

