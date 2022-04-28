import Tkinter as tk
import numpy as np
import pandas as pd
import math
from decimal import Decimal, ROUND_HALF_UP
import Fuchs


window=tk.Tk()
window.title("Charge distribution calculater")
window.geometry("940x470")

Base=40
A=20
x1=10
x2=150
x3=250
Boxsize=10

label1=tk.Label(text="   Physical properties")
label1.place(x=x1,y=Base-20)

p_l=tk.Label(text="Pressure")
p_e=tk.Entry(width=Boxsize)
p_u=tk.Label(text="Pa")
p_l.place(x=x1,y=Base)
p_e.place(x=x2,y=Base)
p_u.place(x=x3,y=Base)

T_l=tk.Label(text="Temperature")
T_e=tk.Entry(width=Boxsize)
T_u=tk.Label(text="K")
T_l.place(x=x1,y=Base+A*1)
T_e.place(x=x2,y=Base+A*1)
T_u.place(x=x3,y=Base+A*1)

Mgas_l=tk.Label(text="Mass of gas")
Mgas_e=tk.Entry(width=Boxsize)
Mgas_u=tk.Label(text="g/mol")
Mgas_l.place(x=x1,y=Base+A*2)
Mgas_e.place(x=x2,y=Base+A*2)
Mgas_u.place(x=x3,y=Base+A*2)

myu_l=tk.Label(text="Gas viscosity")
myu_e=tk.Entry(width=Boxsize)
myu_u=tk.Label(text="Pas")
myu_l.place(x=x1,y=Base+A*3)
myu_e.place(x=x2,y=Base+A*3)
myu_u.place(x=x3,y=Base+A*3)

dp_l=tk.Label(text="Particle diameter")
dp_e=tk.Entry(width=Boxsize)
dp_u=tk.Label(text="m")
dp_l.place(x=x1,y=Base+A*4)
dp_e.place(x=x2,y=Base+A*4)
dp_u.place(x=x3,y=Base+A*4)

rhop_l=tk.Label(text="Particle density")
rhop_e=tk.Entry(width=Boxsize)
rhop_u=tk.Label(text="kg/m3")
rhop_l.place(x=x1,y=Base+A*5)
rhop_e.place(x=x2,y=Base+A*5)
rhop_u.place(x=x3,y=Base+A*5)

Dp_l=tk.Label(text="Particle diffusion")
Dp_e=tk.Entry(width=Boxsize)
Dp_u=tk.Label(text="m2/s")
Dp_l.place(x=x1,y=Base+A*6)
Dp_e.place(x=x2,y=Base+A*6)
Dp_u.place(x=x3,y=Base+A*6)

ep_l=tk.Label(text="Particle permittivity")
ep_e=tk.Entry(width=Boxsize)
ep_u=tk.Label(text="F/m")
ep_l.place(x=x1,y=Base+A*7)
ep_e.place(x=x2,y=Base+A*7)
ep_u.place(x=x3,y=Base+A*7)

dt_l=tk.Label(text="Time step")
dt_e=tk.Entry(width=Boxsize)
dt_u=tk.Label(text="s")
dt_l.place(x=x1,y=Base+A*8)
dt_e.place(x=x2,y=Base+A*8)
dt_u.place(x=x3,y=Base+A*8)

tmax_l=tk.Label(text="tmax")
tmax_e=tk.Entry(width=Boxsize)
tmax_u=tk.Label(text="s")
tmax_l.place(x=x1,y=Base+A*9)
tmax_e.place(x=x2,y=Base+A*9)
tmax_u.place(x=x3,y=Base+A*9)

Zpionm_l=tk.Label(text="Ion mobility (-)")
Zpionm_e=tk.Entry(width=Boxsize)
Zpionm_u=tk.Label(text="m2/Vs")
Zpionm_l.place(x=x1,y=Base+A*10)
Zpionm_e.place(x=x2,y=Base+A*10)
Zpionm_u.place(x=x3,y=Base+A*10)

Zpionp_l=tk.Label(text="Ion mobility (+)")
Zpionp_e=tk.Entry(width=Boxsize)
Zpionp_u=tk.Label(text="m2/Vs")
Zpionp_l.place(x=x1,y=Base+A*11)
Zpionp_e.place(x=x2,y=Base+A*11)
Zpionp_u.place(x=x3,y=Base+A*11)

mionm_l=tk.Label(text="Ion mass (-)")
mionm_e=tk.Entry(width=Boxsize)
mionm_u=tk.Label(text="g/mol")
mionm_l.place(x=x1,y=Base+A*12)
mionm_e.place(x=x2,y=Base+A*12)
mionm_u.place(x=x3,y=Base+A*12)

mionp_l=tk.Label(text="Ion mass (+)")
mionp_e=tk.Entry(width=Boxsize)
mionp_u=tk.Label(text="g/mol")
mionp_l.place(x=x1,y=Base+A*13)
mionp_e.place(x=x2,y=Base+A*13)
mionp_u.place(x=x3,y=Base+A*13)

Cionm_l=tk.Label(text="Ion concentration (-)")
Cionm_e=tk.Entry(width=Boxsize)
Cionm_u=tk.Label(text="cc-1")
Cionm_l.place(x=x1,y=Base+A*14)
Cionm_e.place(x=x2,y=Base+A*14)
Cionm_u.place(x=x3,y=Base+A*14)

Cionp_l=tk.Label(text="Ion concentration (+)")
Cionp_e=tk.Entry(width=Boxsize)
Cionp_u=tk.Label(text="cc-1")
Cionp_l.place(x=x1,y=Base+A*15)
Cionp_e.place(x=x2,y=Base+A*15)
Cionp_u.place(x=x3,y=Base+A*15)


Base=40
A=20
x1=350
x2=400
x3=500
x4=720
Boxsize=25

label1=tk.Label(text="C0 [a.u.]")
label1.place(x=x2,y=Base-20)
label1=tk.Label(text="Beta(-) [m3/s]")
label1.place(x=x3,y=Base-20)
label1=tk.Label(text="Beta(+) [m3/s]")
label1.place(x=x4,y=Base-20)

C0_l=tk.Label(text="z=-7")
C0_e=tk.Entry(width=10)
beta0m_e=tk.Entry(width=Boxsize)
beta0p_e=tk.Entry(width=Boxsize)
C0_l.place(x=x1,y=Base)
C0_e.place(x=x2,y=Base)
beta0m_e.place(x=x3,y=Base)
beta0p_e.place(x=x4,y=Base)

C1_l=tk.Label(text="z=-6")
C1_e=tk.Entry(width=10)
beta1m_e=tk.Entry(width=Boxsize)
beta1p_e=tk.Entry(width=Boxsize)
C1_l.place(x=x1,y=Base+A*1)
C1_e.place(x=x2,y=Base+A*1)
beta1m_e.place(x=x3,y=Base+A*1)
beta1p_e.place(x=x4,y=Base+A*1)

C2_l=tk.Label(text="z=-5")
C2_e=tk.Entry(width=10)
beta2m_e=tk.Entry(width=Boxsize)
beta2p_e=tk.Entry(width=Boxsize)
C2_l.place(x=x1,y=Base+A*2)
C2_e.place(x=x2,y=Base+A*2)
beta2m_e.place(x=x3,y=Base+A*2)
beta2p_e.place(x=x4,y=Base+A*2)

C3_l=tk.Label(text="z=-4")
C3_e=tk.Entry(width=10)
beta3m_e=tk.Entry(width=Boxsize)
beta3p_e=tk.Entry(width=Boxsize)
C3_l.place(x=x1,y=Base+A*3)
C3_e.place(x=x2,y=Base+A*3)
beta3m_e.place(x=x3,y=Base+A*3)
beta3p_e.place(x=x4,y=Base+A*3)

C4_l=tk.Label(text="z=-3")
C4_e=tk.Entry(width=10)
beta4m_e=tk.Entry(width=Boxsize)
beta4p_e=tk.Entry(width=Boxsize)
C4_l.place(x=x1,y=Base+A*4)
C4_e.place(x=x2,y=Base+A*4)
beta4m_e.place(x=x3,y=Base+A*4)
beta4p_e.place(x=x4,y=Base+A*4)

C5_l=tk.Label(text="z=-2")
C5_e=tk.Entry(width=10)
beta5m_e=tk.Entry(width=Boxsize)
beta5p_e=tk.Entry(width=Boxsize)
C5_l.place(x=x1,y=Base+A*5)
C5_e.place(x=x2,y=Base+A*5)
beta5m_e.place(x=x3,y=Base+A*5)
beta5p_e.place(x=x4,y=Base+A*5)

C6_l=tk.Label(text="z=-1")
C6_e=tk.Entry(width=10)
beta6m_e=tk.Entry(width=Boxsize)
beta6p_e=tk.Entry(width=Boxsize)
C6_l.place(x=x1,y=Base+A*6)
C6_e.place(x=x2,y=Base+A*6)
beta6m_e.place(x=x3,y=Base+A*6)
beta6p_e.place(x=x4,y=Base+A*6)

C7_l=tk.Label(text="z=0")
C7_e=tk.Entry(width=10)
beta7m_e=tk.Entry(width=Boxsize)
beta7p_e=tk.Entry(width=Boxsize)
C7_l.place(x=x1,y=Base+A*7)
C7_e.place(x=x2,y=Base+A*7)
beta7m_e.place(x=x3,y=Base+A*7)
beta7p_e.place(x=x4,y=Base+A*7)

C8_l=tk.Label(text="z=+1")
C8_e=tk.Entry(width=10)
beta8m_e=tk.Entry(width=Boxsize)
beta8p_e=tk.Entry(width=Boxsize)
C8_l.place(x=x1,y=Base+A*8)
C8_e.place(x=x2,y=Base+A*8)
beta8m_e.place(x=x3,y=Base+A*8)
beta8p_e.place(x=x4,y=Base+A*8)

C9_l=tk.Label(text="z=+2")
C9_e=tk.Entry(width=10)
beta9m_e=tk.Entry(width=Boxsize)
beta9p_e=tk.Entry(width=Boxsize)
C9_l.place(x=x1,y=Base+A*9)
C9_e.place(x=x2,y=Base+A*9)
beta9m_e.place(x=x3,y=Base+A*9)
beta9p_e.place(x=x4,y=Base+A*9)

C10_l=tk.Label(text="z=+3")
C10_e=tk.Entry(width=10)
beta10m_e=tk.Entry(width=Boxsize)
beta10p_e=tk.Entry(width=Boxsize)
C10_l.place(x=x1,y=Base+A*10)
C10_e.place(x=x2,y=Base+A*10)
beta10m_e.place(x=x3,y=Base+A*10)
beta10p_e.place(x=x4,y=Base+A*10)

C11_l=tk.Label(text="z=+4")
C11_e=tk.Entry(width=10)
beta11m_e=tk.Entry(width=Boxsize)
beta11p_e=tk.Entry(width=Boxsize)
C11_l.place(x=x1,y=Base+A*11)
C11_e.place(x=x2,y=Base+A*11)
beta11m_e.place(x=x3,y=Base+A*11)
beta11p_e.place(x=x4,y=Base+A*11)

C12_l=tk.Label(text="z=+5")
C12_e=tk.Entry(width=10)
beta12m_e=tk.Entry(width=Boxsize)
beta12p_e=tk.Entry(width=Boxsize)
C12_l.place(x=x1,y=Base+A*12)
C12_e.place(x=x2,y=Base+A*12)
beta12m_e.place(x=x3,y=Base+A*12)
beta12p_e.place(x=x4,y=Base+A*12)

C13_l=tk.Label(text="z=+6")
C13_e=tk.Entry(width=10)
beta13m_e=tk.Entry(width=Boxsize)
beta13p_e=tk.Entry(width=Boxsize)
C13_l.place(x=x1,y=Base+A*13)
C13_e.place(x=x2,y=Base+A*13)
beta13m_e.place(x=x3,y=Base+A*13)
beta13p_e.place(x=x4,y=Base+A*13)

C14_l=tk.Label(text="z=+7")
C14_e=tk.Entry(width=10)
beta14m_e=tk.Entry(width=Boxsize)
beta14p_e=tk.Entry(width=Boxsize)
C14_l.place(x=x1,y=Base+A*14)
C14_e.place(x=x2,y=Base+A*14)
beta14m_e.place(x=x3,y=Base+A*14)
beta14p_e.place(x=x4,y=Base+A*14)


def setPP():
	Fuchs.p=float(p_e.get())
	Fuchs.T=float(T_e.get())
	Fuchs.Mgas=float(Mgas_e.get())
	Fuchs.myu=float(myu_e.get())
	Fuchs.dp=float(dp_e.get())
	Fuchs.rhop=float(rhop_e.get())
	Fuchs.Dp=float(Dp_e.get())
	Fuchs.ep=float(ep_e.get())
	Fuchs.dt=float(dt_e.get())
	Fuchs.tmax=float(tmax_e.get())
	Fuchs.Kion[0]=float(Zpionp_e.get())
	Fuchs.Kion[1]=float(Zpionm_e.get())
	Fuchs.Mion[0]=float(mionp_e.get())
	Fuchs.Mion[1]=float(mionm_e.get())
	Fuchs.Cion[0]=float(Cionp_e.get())
	Fuchs.Cion[1]=float(Cionm_e.get())

def	NormalCond():
	p_e.delete(0,tk.END)
	T_e.delete(0,tk.END)
	Mgas_e.delete(0,tk.END)
	myu_e.delete(0,tk.END)
	dp_e.delete(0,tk.END)
	rhop_e.delete(0,tk.END)
	Dp_e.delete(0,tk.END)
	ep_e.delete(0,tk.END)
	dt_e.delete(0,tk.END)
	tmax_e.delete(0,tk.END)
	Zpionm_e.delete(0,tk.END)
	Zpionp_e.delete(0,tk.END)
	mionm_e.delete(0,tk.END)
	mionp_e.delete(0,tk.END)
	Cionm_e.delete(0,tk.END)
	Cionp_e.delete(0,tk.END)

	p_e.insert(tk.END,101300)
	T_e.insert(tk.END,298.15)
	Mgas_e.insert(tk.END,28.8)
	myu_e.insert(tk.END,1.8e-5)
	dp_e.insert(tk.END,100e-9)
	rhop_e.insert(tk.END,1000)
	Dp_e.insert(tk.END,7.0825e-10)
	ep_e.insert(tk.END,100000)
	dt_e.insert(tk.END,1e-3)
	tmax_e.insert(tk.END,300)
	Zpionm_e.insert(tk.END,1.9e-4)
	Zpionp_e.insert(tk.END,1.4e-4)
	mionm_e.insert(tk.END,50)
	mionp_e.insert(tk.END,100)
	Cionm_e.insert(tk.END,1e10)
	Cionp_e.insert(tk.END,1e10)


def	setC0():
	C0_e.delete(0,tk.END)
	C1_e.delete(0,tk.END)
	C2_e.delete(0,tk.END)
	C3_e.delete(0,tk.END)
	C4_e.delete(0,tk.END)
	C5_e.delete(0,tk.END)
	C6_e.delete(0,tk.END)
	C7_e.delete(0,tk.END)
	C8_e.delete(0,tk.END)
	C9_e.delete(0,tk.END)
	C10_e.delete(0,tk.END)
	C11_e.delete(0,tk.END)
	C12_e.delete(0,tk.END)
	C13_e.delete(0,tk.END)
	C14_e.delete(0,tk.END)

	C0_e.insert(tk.END,0)
	C1_e.insert(tk.END,0)
	C2_e.insert(tk.END,0)
	C3_e.insert(tk.END,0)
	C4_e.insert(tk.END,0)
	C5_e.insert(tk.END,0)
	C6_e.insert(tk.END,0)
	C7_e.insert(tk.END,1)
	C8_e.insert(tk.END,0)
	C9_e.insert(tk.END,0)
	C10_e.insert(tk.END,0)
	C11_e.insert(tk.END,0)
	C12_e.insert(tk.END,0)
	C13_e.insert(tk.END,0)
	C14_e.insert(tk.END,0)

def	setbetaNeg():
	beta0m_e.delete(0,tk.END)
	beta1m_e.delete(0,tk.END)
	beta2m_e.delete(0,tk.END)
	beta3m_e.delete(0,tk.END)
	beta4m_e.delete(0,tk.END)
	beta5m_e.delete(0,tk.END)
	beta6m_e.delete(0,tk.END)
	beta7m_e.delete(0,tk.END)
	beta8m_e.delete(0,tk.END)
	beta9m_e.delete(0,tk.END)
	beta10m_e.delete(0,tk.END)
	beta11m_e.delete(0,tk.END)
	beta12m_e.delete(0,tk.END)
	beta13m_e.delete(0,tk.END)
	beta14m_e.delete(0,tk.END)

	beta0m_e.insert(tk.END,Fuchs.beta[0][1])
	beta1m_e.insert(tk.END,Fuchs.beta[1][1])
	beta2m_e.insert(tk.END,Fuchs.beta[2][1])
	beta3m_e.insert(tk.END,Fuchs.beta[3][1])
	beta4m_e.insert(tk.END,Fuchs.beta[4][1])
	beta5m_e.insert(tk.END,Fuchs.beta[5][1])
	beta6m_e.insert(tk.END,Fuchs.beta[6][1])
	beta7m_e.insert(tk.END,Fuchs.beta[7][1])
	beta8m_e.insert(tk.END,Fuchs.beta[8][1])
	beta9m_e.insert(tk.END,Fuchs.beta[9][1])
	beta10m_e.insert(tk.END,Fuchs.beta[10][1])
	beta11m_e.insert(tk.END,Fuchs.beta[11][1])
	beta12m_e.insert(tk.END,Fuchs.beta[12][1])
	beta13m_e.insert(tk.END,Fuchs.beta[13][1])
	beta14m_e.insert(tk.END,Fuchs.beta[14][1])


def	setbetaPos():
	beta0p_e.delete(0,tk.END)
	beta1p_e.delete(0,tk.END)
	beta2p_e.delete(0,tk.END)
	beta3p_e.delete(0,tk.END)
	beta4p_e.delete(0,tk.END)
	beta5p_e.delete(0,tk.END)
	beta6p_e.delete(0,tk.END)
	beta7p_e.delete(0,tk.END)
	beta8p_e.delete(0,tk.END)
	beta9p_e.delete(0,tk.END)
	beta10p_e.delete(0,tk.END)
	beta11p_e.delete(0,tk.END)
	beta12p_e.delete(0,tk.END)
	beta13p_e.delete(0,tk.END)
	beta14p_e.delete(0,tk.END)

	beta0p_e.insert(tk.END,Fuchs.beta[0][0])
	beta1p_e.insert(tk.END,Fuchs.beta[1][0])
	beta2p_e.insert(tk.END,Fuchs.beta[2][0])
	beta3p_e.insert(tk.END,Fuchs.beta[3][0])
	beta4p_e.insert(tk.END,Fuchs.beta[4][0])
	beta5p_e.insert(tk.END,Fuchs.beta[5][0])
	beta6p_e.insert(tk.END,Fuchs.beta[6][0])
	beta7p_e.insert(tk.END,Fuchs.beta[7][0])
	beta8p_e.insert(tk.END,Fuchs.beta[8][0])
	beta9p_e.insert(tk.END,Fuchs.beta[9][0])
	beta10p_e.insert(tk.END,Fuchs.beta[10][0])
	beta11p_e.insert(tk.END,Fuchs.beta[11][0])
	beta12p_e.insert(tk.END,Fuchs.beta[12][0])
	beta13p_e.insert(tk.END,Fuchs.beta[13][0])
	beta14p_e.insert(tk.END,Fuchs.beta[14][0])


def	calculateFuchs():
	setPP()
	IF=0
	if bln.get():
		IF=1
	else:
		IF=0
	Fuchs.Fuchs(IF)
	setbetaNeg()
	setbetaPos()

def	calculateCG():
	setPP()
	Fuchs.CG()
	setbetaNeg()
	setbetaPos()

def	electron():
	setPP()
	Fuchs.electron()
	setbetaNeg()

def timedevelop():
	Fuchs.C[0]=float(C0_e.get())
	Fuchs.C[1]=float(C1_e.get())
	Fuchs.C[2]=float(C2_e.get())
	Fuchs.C[3]=float(C3_e.get())
	Fuchs.C[4]=float(C4_e.get())
	Fuchs.C[5]=float(C5_e.get())
	Fuchs.C[6]=float(C6_e.get())
	Fuchs.C[7]=float(C7_e.get())
	Fuchs.C[8]=float(C8_e.get())
	Fuchs.C[9]=float(C9_e.get())
	Fuchs.C[10]=float(C10_e.get())
	Fuchs.C[11]=float(C11_e.get())
	Fuchs.C[12]=float(C12_e.get())
	Fuchs.C[13]=float(C13_e.get())
	Fuchs.C[14]=float(C14_e.get())

	Fuchs.beta[0][0]=float(beta0p_e.get())
	Fuchs.beta[1][0]=float(beta1p_e.get())
	Fuchs.beta[2][0]=float(beta2p_e.get())
	Fuchs.beta[3][0]=float(beta3p_e.get())
	Fuchs.beta[4][0]=float(beta4p_e.get())
	Fuchs.beta[5][0]=float(beta5p_e.get())
	Fuchs.beta[6][0]=float(beta6p_e.get())
	Fuchs.beta[7][0]=float(beta7p_e.get())
	Fuchs.beta[8][0]=float(beta8p_e.get())
	Fuchs.beta[9][0]=float(beta9p_e.get())
	Fuchs.beta[10][0]=float(beta10p_e.get())
	Fuchs.beta[11][0]=float(beta11p_e.get())
	Fuchs.beta[12][0]=float(beta12p_e.get())
	Fuchs.beta[13][0]=float(beta13p_e.get())
	Fuchs.beta[14][0]=float(beta14p_e.get())

	Fuchs.beta[0][1]=float(beta0m_e.get())
	Fuchs.beta[1][1]=float(beta1m_e.get())
	Fuchs.beta[2][1]=float(beta2m_e.get())
	Fuchs.beta[3][1]=float(beta3m_e.get())
	Fuchs.beta[4][1]=float(beta4m_e.get())
	Fuchs.beta[5][1]=float(beta5m_e.get())
	Fuchs.beta[6][1]=float(beta6m_e.get())
	Fuchs.beta[7][1]=float(beta7m_e.get())
	Fuchs.beta[8][1]=float(beta8m_e.get())
	Fuchs.beta[9][1]=float(beta9m_e.get())
	Fuchs.beta[10][1]=float(beta10m_e.get())
	Fuchs.beta[11][1]=float(beta11m_e.get())
	Fuchs.beta[12][1]=float(beta12m_e.get())
	Fuchs.beta[13][1]=float(beta13m_e.get())
	Fuchs.beta[14][1]=float(beta14m_e.get())
	setPP()
	Fuchs.Time_develop()


def calcDiff():
	setPP()
	Dp_e.delete(0,tk.END)
	Dp_e.insert(tk.END,Fuchs.calcDiff())



input_FM=tk.Button(window,text="Electron",width=5,height=1,command=lambda:electron())
input_FM.place(x=500, y=340)


time_development=tk.Button(window,text="Solve dn/dt",width=30,height=2,command=lambda:timedevelop())
time_development.place(x=600, y=420)


input_C0=tk.Button(window,text="Set C0",width=12,height=2,command=lambda:setC0())
input_C0.place(x=345, y=370)


input_Fuchs=tk.Button(window,text="Set Beta Fuchs",width=12,height=2,command=lambda:calculateFuchs())
input_Fuchs.place(x=480, y=370)


input_CG=tk.Button(window,text="Set Beta CG",width=12,height=2,command=lambda:calculateCG())
input_CG.place(x=615, y=370)

input_NC=tk.Button(window,text="Set normal air conditions",width=30,height=2,command=lambda:NormalCond())
input_NC.place(x=10, y=370)

diffusion_calcu=tk.Button(window,text="Calculate diffusion coefficient",width=30,height=2,command=lambda:calcDiff())
diffusion_calcu.place(x=10, y=420)

time_development=tk.Button(window,text="Solve dn/dt",width=30,height=2,command=lambda:timedevelop())
time_development.place(x=600, y=420)


bln=tk.BooleanVar()
imageforce=tk.Checkbutton(variable=bln, text="Image force")
imageforce.place(x=750, y=370)


window.mainloop()






