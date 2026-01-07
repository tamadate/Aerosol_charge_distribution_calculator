import numpy as np
import matplotlib.pyplot as plt
import functions
import physicalProp
import Fuchs
import LD



class collisionRate:
	def __init__(self):
		self.PPs = physicalProp.PP()
		self.C = np.zeros(self.PPs.nZ)
		self.dt = 1e-8
		self.tmax = 1

		self.Fuchs = Fuchs.Fuchs()
		self.LD = LD.LD()

	def computeFuchsBeta(self):
		self.Fuchs.computeBeta(self.PPs)

	def computeLDBeta(self):
		self.LD.computeBetaLD(self.PPs)


	def showParameters(self):
		print("=== PPs parameters ===")
		for k, v in vars(self.PPs).items():
			print(f"{k:15s} : {v}")



	## Population balence equation
	def func_pop(self, C0, C1, C_1, i):
		term1 = - C0 * self.PPs.Cion[0] * self.PPs.beta[i][0]
		term2 = - C0 * self.PPs.Cion[1] * self.PPs.beta[i][1]
		term3 = C_1 * self.PPs.Cion[0] * self.PPs.beta[i - 1][0]
		term4 = C1 * self.PPs.Cion[1] * self.PPs.beta[i + 1][1]
		return term1 + term2 + term3 + term4

	def Time_develop(self, save_every=1000, plot=True, out_prefix="chargeDistribution"):
		nZ = self.PPs.nZ       
		idx = range(1, nZ - 1)   

		# --- 初期正規化 ---
		CT = float(np.sum(self.C))
		Cs_list = [self.C / CT]
		ts_list = [0.0]

		# --- 事前確保（毎ステップ確保を避ける） ---
		K1 = np.zeros(nZ)
		K2 = np.zeros(nZ)
		K3 = np.zeros(nZ)
		K4 = np.zeros(nZ)

		if plot:
			plt.ion()
			fig, ax = plt.subplots()

		nsteps = int(np.floor(self.tmax / self.dt))
		for step in range(nsteps):
			t = step * self.dt

			# ----- RK4（4回のforを1回にまとめる） -----
			for i in idx:
				K1[i] = self.dt * self.func_pop(C0 = self.C[i], C1 = self.C[i+1], C_1 = self.C[i-1], i = i)

			for i in idx:
				ci  = self.C[i]   + 0.5 * K1[i]
				cip = self.C[i + 1] + 0.5 * K1[i + 1]
				cim = self.C[i - 1] + 0.5 * K1[i - 1]
				K2[i] = self.dt * self.func_pop(ci, cip, cim, i)

			for i in idx:
				ci  = self.C[i]   + 0.5 * K2[i]
				cip = self.C[i + 1 ] + 0.5 * K2[i + 1]
				cim = self.C[i - 1] + 0.5 * K2[i - 1]
				K3[i] = self.dt * self.func_pop(ci, cip, cim, i)

			for i in idx:
				ci  = self.C[i]   + K3[i]
				cip = self.C[i+1] + K3[i+1]
				cim = self.C[i-1] + K3[i-1]
				K4[i] = self.dt * self.func_pop(ci, cip, cim, i)

			for i in idx:
				self.C[i] += (K1[i] + 2*K2[i] + 2*K3[i] + K4[i]) / 6.0

			# ----- 出力（間引き） -----
			if (step % save_every) == 0:
				CT = float(np.sum(self.C))
				Cs_list.append(self.C / CT)
				ts_list.append(t)

				if plot:
					ax.cla()
					Cs_arr = np.asarray(Cs_list)     # shape: (nsave, nZ)
					ts_arr = np.asarray(ts_list)

					for i in range(nZ):
						# 最新時刻で 1e-3 以上だけ描く（あなたの条件を維持）
						if Cs_arr[-1, i] > 1e-3:
							ax.plot(ts_arr, Cs_arr[:, i], label=f"z={i-(nZ//2)}")

					ax.set_xlabel("Time [s]")
					ax.set_ylabel("Normalized concentration [-]")
					ax.legend(loc="upper right")
					plt.pause(0.01)

		# --- numpy 配列として返す & 保存 ---
		ts = np.asarray(ts_list)          # (nsave,)
		Cs = np.asarray(Cs_list)          # (nsave, nZ)

		# time + Cs を横に結合
		out = np.column_stack((ts, Cs))

		# ヘッダー作成
		z0 = self.PPs.nZ // 2
		header = "time," + ",".join([f"z={i-z0}" for i in range(self.PPs.nZ)])

		np.savetxt(
			"chargeDistribution.csv",
			out,
			delimiter=",",
			header=header,
			comments=""
		)

		return ts, Cs

'''
	def equilibrium(self):   
		idx = range(0, self.PPs.N - 1)   

		# --- 初期正規化 ---
		CT = float(np.sum(self.C))
		Cs_list = [self.C / CT]
		ts_list = [0.0]

		# ----- RK4（4回のforを1回にまとめる） -----
		self.A_p = np.zeros(self.PPs.N)
		self.A_m = np.zeros(self.PPs.N)
		self.A_0 = 0
		offset_p = 8
		offset_m = 6
		
		for i in np.arange(10):
			for i in idx:
				ip = offset_p + i
				im = offset_m - i
				ratio_p = self.PPs.Cion[0] / self.PPs.Cion[1]
				self.A_p[i] = ratio_p * self.PPs.beta[ip - 1][0] / (self.PPs.beta[ip][1] + ratio_p * 
														self.PPs.beta[ip][0] - self.PPs.beta[ip + 1][1] * self.C[ip + 1] / (self.C[ip] + 1e-10))
				ratio_m = self.PPs.Cion[1] / self.PPs.Cion[0]
				self.A_m[i] = ratio_m * self.PPs.beta[im - 1][1] / (self.PPs.beta[im][0] + ratio_m * 
														self.PPs.beta[im][1] - self.PPs.beta[im + 1][0] * self.C[im + 0] / (self.C[im] + 1e-10))
			A_p_fact = np.array([1])
			A_m_fact = np.array([1])
			Csum = 0
			for i in idx:
				ip = offset_p + i
				im = offset_m - i
				A_p_fact = np.append(A_p_fact, A_p_fact[-1] * self.A_p[i])
				A_m_fact = np.append(A_m_fact, A_m_fact[-1] * self.A_m[i])
				self.C[ip] = A_p_fact[-1] / (1 + np.sum(A_p_fact[1:] + A_m_fact[1:]))
				self.C[im] = A_m_fact[-1] / (1 + np.sum(A_p_fact[1:] + A_m_fact[1:]))
				Csum += self.C[ip] + self.C[im]
			self.C[7] = 1 - Csum
			print(self.C)
'''


