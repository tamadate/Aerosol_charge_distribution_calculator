import tkinter as tk
import numpy as np
from decimal import Decimal, ROUND_HALF_UP
import collisionRate
import functions

collisionRate=collisionRate.collisionRate()

#Generate window
window=tk.Tk()
window.title("Charge distribution calculator")
window.geometry("750x470")

#Generate physical property boxes
Base=40
A=20
x1=10
x2=150
x3=250
Boxsize=10

label1=tk.Label(text="   Physical properties")
label1.place(x=x1,y=Base-20)

labelNames=["Pressure","Temperature","Mass of gas","Gas viscosity","Particle diameter","Particle density","Particle diffusion","Particle permittivity","Time step","Total time","Ion mobility (-)","Ion mobility (+)","ion mass (-)","ion mass (+)","ion concentration (-)","ion concentration (+)"]
labelUnits=["atm","K","g/mol","Pas","um","kg/m3","cm2/s","F/m","s","s","cm2/Vs","cm2/Vs","g/mol","g/mol","cc-1","cc-1"]
initials=[1,298.15,28.8,1.8e-5,0.1,1000,7.08e-6,1e5,1e-3,300,1.9,1.4,50,100,1e10,1e10]
ppEntry=[]
for i in np.arange(np.size(labelNames)):
	label=tk.Label(text=labelNames[i])
	label.place(x=x1,y=Base+A*i)
	ppEntry=np.append(ppEntry,tk.Entry(width=Boxsize))
	ppEntry[i].place(x=x2,y=Base+A*i)
	label=tk.Label(text=labelUnits[i])
	label.place(x=x3,y=Base+A*i)


#Generate initial condition and beta boxes
x1=350
x2=400
x3=500
x4=600

label1=tk.Label(text="N0/Nt [-]")
label1.place(x=x2,y=Base-20)
label1=tk.Label(text="Beta(-) [m3/s]")
label1.place(x=x3,y=Base-20)
label1=tk.Label(text="Beta(+) [m3/s]")
label1.place(x=x4,y=Base-20)

initials=[1,298.15,28.8,1.8e-5,0.1,1000,7.08e-6,1e5,1e-6,0.1,1.9,1.4,50,100,1e8,1e8]
CEntry=[]
betaEntry_p=[]
betaEntry_m=[]

for i in np.arange(15):
	label=tk.Label(text="z="+str(i-7))
	label.place(x=x1,y=Base+A*i)
	CEntry=np.append(CEntry,tk.Entry(width=Boxsize))
	CEntry[i]=tk.Entry(width=Boxsize)
	CEntry[i].place(x=x2,y=Base+A*i)
	betaEntry_m=np.append(betaEntry_m,tk.Entry(width=Boxsize))
	betaEntry_m[i]=tk.Entry(width=Boxsize)
	betaEntry_m[i].place(x=x3,y=Base+A*i)
	betaEntry_p=np.append(betaEntry_p,tk.Entry(width=Boxsize))
	betaEntry_p[i]=tk.Entry(width=Boxsize)
	betaEntry_p[i].place(x=x4,y=Base+A*i)


#Get numbers from input boxes
def setPP():
	collisionRate.PPs.p=float(ppEntry[0].get())*101300
	collisionRate.PPs.T=float(ppEntry[1].get())
	collisionRate.PPs.Mgas=float(ppEntry[2].get())
	collisionRate.PPs.myu=float(ppEntry[3].get())
	collisionRate.PPs.dp=float(ppEntry[4].get())*1e-6
	collisionRate.PPs.rhop=float(ppEntry[5].get())
	collisionRate.PPs.Dp=float(ppEntry[6].get())*1e-4
	collisionRate.dt=float(ppEntry[8].get())
	collisionRate.tmax=float(ppEntry[9].get())
	collisionRate.PPs.Kion[1]=float(ppEntry[10].get())*1e-4
	collisionRate.PPs.Kion[0]=float(ppEntry[11].get())*1e-4
	collisionRate.PPs.Mion[1]=float(ppEntry[12].get())
	collisionRate.PPs.Mion[0]=float(ppEntry[13].get())
	collisionRate.PPs.Cion[1]=float(ppEntry[14].get())*1e6
	collisionRate.PPs.Cion[0]=float(ppEntry[15].get())*1e6
	if bln.get():
		collisionRate.PPs.ep=float(ppEntry[7].get())
	else:
		collisionRate.PPs.ep=0

#Put normal conditions into the boxes
def	NormalCond():
	for i in np.arange(np.size(ppEntry)):
		ppEntry[i].delete(0,tk.END)
		ppEntry[i].insert(tk.END,'{:.5g}'.format(initials[i]))


#Put initial conditions
def	setC0():
	for i in np.arange(15):
		CEntry[i].delete(0,tk.END)
		CEntry[i].insert(tk.END,0)
	CEntry[7].delete(0,tk.END)
	CEntry[7].insert(tk.END,1)

#Put beta from calculation in Fuchs.py
def	setbeta():
	for i in np.arange(15):
		betaEntry_m[i].delete(0,tk.END)
		betaEntry_m[i].insert(tk.END,'{:.5g}'.format(collisionRate.PPs.beta[i][1]))
		betaEntry_p[i].delete(0,tk.END)
		betaEntry_p[i].insert(tk.END,'{:.5g}'.format(collisionRate.PPs.beta[i][0]))

def	calculateFuchs():
	setPP()
	collisionRate.computeFuchsBeta()
	setbeta()


def	calculateCG():
	setPP()
	collisionRate.computeLDBeta()
	setbeta()


def timedevelop():
	for i in np.arange(15):
		collisionRate.C[i]=float(CEntry[i].get())
		collisionRate.PPs.beta[i][0]=float(betaEntry_p[i].get())
		collisionRate.PPs.beta[i][1]=float(betaEntry_m[i].get())
	setPP()
	collisionRate.Time_develop()

def calcDiff():
	setPP()
	ppEntry[6].delete(0,tk.END)
	ppEntry[6].insert(tk.END,'{:.5g}'.format(functions.Diff(collisionRate.PPs.myu,collisionRate.PPs.dp,collisionRate.PPs.T)*1e4))


#Generate buttons
input_C0=tk.Button(window,text="Set N0/Nt",width=12,height=2,command=lambda:setC0())
input_C0.place(x=345, y=370)

input_Fuchs=tk.Button(window,text="Set Beta Fuchs",width=12,height=2,command=lambda:calculateFuchs())
input_Fuchs.place(x=480, y=370)

input_CG=tk.Button(window,text="Set Beta LD",width=12,height=2,command=lambda:calculateCG())
input_CG.place(x=615, y=370)

input_NC=tk.Button(window,text="Set normal air conditions",width=30,height=2,command=lambda:NormalCond())
input_NC.place(x=10, y=370)

diffusion_calcu=tk.Button(window,text="Calculate diffusion coefficient",width=30,height=2,command=lambda:calcDiff())
diffusion_calcu.place(x=10, y=420)

time_development=tk.Button(window,text="Solve dn/dt",width=30,height=2,command=lambda:timedevelop())
time_development.place(x=400, y=420)

bln=tk.BooleanVar()
imageforce=tk.Checkbutton(variable=bln, text="Image force")
imageforce.place(x=615, y=350)



window.mainloop()






