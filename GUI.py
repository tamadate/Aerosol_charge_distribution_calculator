import Tkinter as tk
import numpy as np
import pandas as pd
import math
import V_Dp
import TF


def	longDMA():
	Qa_e.delete(0,tk.END)
	Qs_e.delete(0,tk.END)
	rin_e.delete(0,tk.END)
	rout_e.delete(0,tk.END)
	L_e.delete(0,tk.END)
	Vmin_e.delete(0,tk.END)
	Vmax_e.delete(0,tk.END)
	center_dp_e.delete(0,tk.END)
	file_e.delete(0,tk.END)

	Qa_e.insert(tk.END,1.5)
	Qs_e.insert(tk.END,15.0)
	rin_e.insert(tk.END,9.37)
	rout_e.insert(tk.END,19.61)
	L_e.insert(tk.END,443.69)
	Vmin_e.insert(tk.END,10)
	Vmax_e.insert(tk.END,10000)
	center_dp_e.insert(tk.END,3.0)
	file_e.insert(tk.END,"a.dat")


def nanoDMA():
	Qa_e.delete(0,tk.END)
	Qs_e.delete(0,tk.END)
	rin_e.delete(0,tk.END)
	rout_e.delete(0,tk.END)
	L_e.delete(0,tk.END)
	Vmin_e.delete(0,tk.END)
	Vmax_e.delete(0,tk.END)
	center_dp_e.delete(0,tk.END)
	file_e.delete(0,tk.END)

	Qa_e.insert(tk.END,1.5)
	Qs_e.insert(tk.END,15.0)
	rin_e.insert(tk.END,9.37)
	rout_e.insert(tk.END,19.05)
	L_e.insert(tk.END,49.87)
	Vmin_e.insert(tk.END,10)
	Vmax_e.insert(tk.END,10000)
	center_dp_e.insert(tk.END,3.0)
	file_e.insert(tk.END,"a.dat")





window=tk.Tk()
window.title("Charge distribution calculater")
window.geometry("350x450")

p_l=tk.Label(text="Pressure [Pa]")
p_e=tk.Entry(width=5)
p_l.place(x=10,y=30)
p_e.place(x=100,y=30)


rin_e=tk.Entry(width=15)
rout_l=tk.Label(text="DMA outer radius")
rout_e=tk.Entry(width=15)
Qa_l=tk.Label(text="Aerosol flow rate")
Qa_e=tk.Entry(width=15)
Qs_l=tk.Label(text="Seath flow rate")
Qs_e=tk.Entry(width=15)
Vmin_l=tk.Label(text="Minimum voltage")
Vmin_e=tk.Entry(width=15)
Vmax_l=tk.Label(text="Maxmun voltage")
Vmax_e=tk.Entry(width=15)
center_dp_l=tk.Label(text="Diameter")
center_dp_e=tk.Entry(width=15)
file_l=tk.Label(text="Output file")
file_e=tk.Entry(width=15)


L_unit=tk.Label(text="mm")
rin_unit=tk.Label(text="mm")
rout_unit=tk.Label(text="mm")
Qa_unit=tk.Label(text="L/min")
Qs_unit=tk.Label(text="L/min")
Vmin_unit=tk.Label(text="V")
Vmax_unit=tk.Label(text="V")
center_dp_unit=tk.Label(text="nm")


L_l.place(x=10,y=10)

L_unit.place(x=8*15+10,y=30)
rin_l.place(x=10,y=60)
rin_e.place(x=10,y=80)
rin_unit.place(x=8*15+10,y=80)
rout_l.place(x=10,y=110)
rout_e.place(x=10,y=130)
rout_unit.place(x=8*15+10,y=130)
Qa_l.place(x=10,y=160)
Qa_e.place(x=10,y=180)
Qa_unit.place(x=8*15+10,y=180)
Qs_l.place(x=10,y=210)
Qs_e.place(x=10,y=230)
Qs_unit.place(x=8*15+10,y=230)

Vmin_l.place(x=10+180,y=10)
Vmin_e.place(x=10+180,y=30)
Vmin_unit.place(x=8*15+10+180,y=30)
Vmax_l.place(x=10+180,y=60)
Vmax_e.place(x=10+180,y=80)
Vmax_unit.place(x=8*15+10+180,y=80)
center_dp_l.place(x=10+180,y=110)
center_dp_e.place(x=10+180,y=130)
center_dp_unit.place(x=8*15+10+180,y=130)
file_l.place(x=10+180,y=160)
file_e.place(x=10+180,y=180)


compute_range=tk.Button(window,text="Total transfer\nfunction",width=15,height=3,command=lambda:TF.peak_FWHM(float(Qs_e.get()),float(Qa_e.get()),float(L_e.get()),float(rin_e.get()),float(rout_e.get()),float(Vmin_e.get()),float(Vmax_e.get()),float(center_dp_e.get()),file_e.get()))
compute_transfer=tk.Button(window,text="Specific transfer\nfunction",width=15,height=3,command=lambda:TF.OMEGA(float(Qs_e.get()),float(Qa_e.get()),float(L_e.get()),float(rin_e.get()),float(rout_e.get()),float(center_dp_e.get()),file_e.get()))

input_nanoDMA=tk.Button(window,text="Nano-DMA",width=15,height=3,command=lambda:nanoDMA())
input_longDMA=tk.Button(window,text="Long-DMA",width=15,height=3,command=lambda:longDMA())
input_slenderDMA=tk.Button(window,text="Slender-DMA",width=15,height=3,command=lambda:SlenderDMA())

compute_range.place(x=10+180, y=220)
compute_transfer.place(x=10+180, y=280)
input_nanoDMA.place(x=10, y=270)
input_longDMA.place(x=10, y=320)
input_slenderDMA.place(x=10, y=370)


window.mainloop()






