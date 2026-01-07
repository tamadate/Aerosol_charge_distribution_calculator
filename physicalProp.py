import numpy as np
import functions

class PP:
	def __init__(self):
		# constants
		self.k = 1.38e-23
		self.Nw = 6.02e23
		self.e = 1.6e-19
		self.e0 = 8.85e-12
		self.e_fourPiE0 = self.e / (4 * np.pi * self.e0)
		self.ee_fourPiE0 = self.e * self.e_fourPiE0

		# initialization
		# Gas
		# Temperature, Pressure, Gas molar mass
		self.T = 300
		self.p = 1e5
		self.Mgas = 28 * 0.8 + 32 * 0.2

		# Ion (+, -)
		# Charges, Concentration, Molar mass, Mobility
		self.zion = np.array((1, -1))
		self.nIon = len(self.zion)   
		self.Cion = np.array((1e15, 1e15))
		self.Mion = np.array((109, 50))
		self.Kion = np.array((1.4e-4, 1.9e-4))

		# Particle
		# Density, Permittivity, Limit of charges, Diameter, Diffusion coefficient
		self.rhop = 1e3
		self.ep = 1000000
		self.N = 7
		self.nZ = 2 * self.N + 1     
		self.dp = 10e-9
		self.Dp = 7.08e-6
		# Charges
		self.z = np.arange(-self.N, self.N + 1)

		# Brank arrays
		# Collision rate coefficient
		self.beta = np.zeros((2 * self.N + 1, 2))

		self.calculatePPs()

	def calculatePPs(self):
		self.mgas = self.Mgas / self.Nw * 1e-3
		self.mion = self.Mion / self.Nw * 1e-3
		self.mp = functions.mp(self.dp, self.rhop)
		self.mig = 1 / (1 / self.mion + 1 / self.mgas)
		self.Dion = self.Kion / self.e * self.k * self.T
		self.cion = (8 * self.k * self.T / np.pi / self.mion)**0.5
		self.ap = self.dp*0.5
		self.Cc = functions.Cc(self.dp)
		self.lamda = 1.329 * self.Dion * (self.mig / self.k / self.T)**0.5
		self.delta = functions.delta(self.ap,self.lamda)
		self.bm = np.zeros((2*self.N+1,2))+1
		self.integral = np.zeros((2*self.N+1,2))

	def resetPPs(self, p, T, Mgas, myu, dp, rhop, Dp, Kion, Mion, Cion, ep):
		self.p=p
		self.T=T
		self.Mgas=Mgas
		self.myu=myu
		self.dp=dp
		self.rhop=rhop
		self.Dp=Dp
		self.Kion=Kion
		self.Mion=Mion
		self.Cion=Cion
		self.ep=ep
		self.calculatePPs()
			


