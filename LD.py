import numpy as np
import math
import functions
import physicalProp
import Fuchs
from scipy.integrate import quad


class LD(Fuchs.Fuchs):
	def F2(self,x,zp,zion,mion,PPs):
		direction=1
		bm=bb=rm=1000e-9
		dr=500e-9
		x2=x*x
		vk=(2*PPs.k*PPs.T/mion)**0.5
		v=vk*x
		v2=v*v
		for loop in np.arange(1,100):
			rm+=dr*direction*(0.9**loop)
			if(rm<PPs.ap):
				rm=PPs.ap+1e-20
				direction*=-1
			dphi=1-self.phi(rm,zp,zion,PPs)/(0.5*mion*v2)
			if(dphi<0):
				bm=0
				break
			b=rm*dphi**0.5
			if (b<bm):
				bm=b
			if (bb<b):
				direction*=-1
			bb=b
		BM=bm/PPs.ap
		return 2*x2*x*np.exp(-x2)*BM*BM

	def computeBetaLD(self,PPs):
		PPs.calculatePPs()
		bm=np.zeros((15,2))+1
		integral=np.zeros((15,2))
		for i in np.arange(2):
			for j in np.arange(15):
				PsiE=-PPs.zion[i]*PPs.z[j]*PPs.ee_fourPiE0/PPs.k/PPs.T/PPs.ap
				if(PPs.ep!=0):
					integral,err=quad(self.F,0,1,args=(PPs.z[j],PPs.zion[i],PPs))
					itaC=integral**-1
					integral,err=quad(self.F2,0,np.inf,args=(PPs.z[j],PPs.zion[i],PPs.mion[i],PPs))
					itaFM=integral
				if(PPs.ep==0):
					if(PsiE==0):
						itaC=1
					else:
						itaC=PsiE/(1-np.exp(-PsiE))
					if(PsiE>0):
						itaFM=1+PsiE
					else:
						itaFM=np.exp(PsiE)
				friction=PPs.k*PPs.T/PPs.Dion[i]
				KnD=(PPs.mion[i]*PPs.k*PPs.T)**0.5*itaC/friction/PPs.ap/itaFM
				H=(4*math.pi*KnD**2+25.836*KnD**3+(8*math.pi)**0.5*11.211*KnD**4)/(1+3.502*KnD+7.211*KnD**2+11.211*KnD**3)
				if(PsiE>0):
					myuA=2.5
					myuB=4.528*np.exp(-1.088*PsiE)+0.7091*np.log(1+1.537*PsiE)
					myuC=11.36*PsiE**0.272-10.33
					myuk=-0.003533*PsiE+0.05971
					inmyu=1+myuk*(np.log(KnD)-myuB)/myuA
					Myu=myuC/myuA*(inmyu)**(-1/myuk-1)*np.exp(-((inmyu)**(-1/myuk)))
					H=np.exp(Myu)*H
				PPs.beta[j][i]=H*friction*PPs.ap**3*itaFM**2/PPs.mion[i]/itaC
			

