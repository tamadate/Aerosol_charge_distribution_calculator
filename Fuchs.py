import numpy as np
import math
import functions
import physicalProp
from scipy.integrate import quad
import numpy as np
from scipy.optimize import minimize_scalar

class Fuchs():
	def phi(self, r, zp, zi, PPs):
		rinv = 1 / r
		firstTerm = rinv*zp*zi
		secondTerm = -PPs.ap * PPs.ap * PPs.ap * 0.5 * rinv * rinv / (r * r - PPs.ap * PPs.ap) * (PPs.ep - 1.0) / (PPs.ep + 1.0)
		return PPs.ee_fourPiE0 * (firstTerm + secondTerm)

	def F(self, x, zp, zi, PPs):
		Phi=self.phi(1/x,zp,zi,PPs)
		return np.exp(Phi/PPs.k/PPs.T)

	def computeBeta(self, PPs):
		def computeDphi(r, delta_i, z_j, zion_i, PPs):
			return 1.0 + (self.phi(delta_i, z_j, zion_i, PPs) - self.phi(r, z_j, zion_i, PPs)) / (PPs.k * PPs.T * 1.5)
		
		def b_of_r(r, delta_i, z_j, zion_i, PPs):
			dphi = computeDphi(r, delta_i, z_j, zion_i, PPs)
			if dphi <= 0.0 or r <= 0.0:
				return np.inf
			return r * np.sqrt(dphi)

		def minimize_b(delta_i, z_j, zion_i, PPs, rmin, rmax):
			# bounded は [rmin, rmax] 内で探索
			res = minimize_scalar(
				lambda r: b_of_r(r, delta_i, z_j, zion_i, PPs),
				bounds=(rmin, rmax),
				method="bounded",
				options={"xatol": 1e-12}  # 必要なら調整
			)
			return res.x, res.fun
		
		PPs.calculatePPs()

		bm = np.ones((PPs.nZ, PPs.nIon))
		integral = np.zeros((PPs.nZ, PPs.nIon))

		for i in range(PPs.nIon):
			delta_i = PPs.delta[i]
			zion_i = PPs.zion[i]
			cion_i = PPs.cion[i]
			Dion_i = PPs.Dion[i]

			for j in range(PPs.nZ):
				z_j = PPs.z[j]
				
				rm, bm = minimize_b(delta_i, z_j, zion_i, PPs, PPs.ap, delta_i)

				integral, err = quad(self.F, 0, PPs.ap / delta_i, args=(z_j, zion_i, PPs))

				pd = (bm / delta_i) ** 2   # bm*bm/PPs.delta/PPs.delta と同じ

				phid_kT = self.phi(delta_i, z_j, zion_i, PPs) / PPs.k / PPs.T

				cpdd =  cion_i * pd * delta_i * delta_i

				PPs.beta[j][i] = math.pi * cpdd * np.exp(-phid_kT) / (1 + np.exp(-phid_kT) * cpdd / 4.0 / Dion_i / PPs.ap * integral)

