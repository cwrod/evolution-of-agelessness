import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

run_panels = ["A","B","C","D","E"]
plot_panels = ["A","B","C","D","E"]


for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]



if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp,extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	b_vals = []
	a_vals = []
	r_vals = list(np.linspace(0.001,bfunc.rmax,num=1000))
	for r_val in r_vals:
		rep1.val = r_val
		b_vals.append(bfunc.ret_val())
		a_vals.append(mat_age.ret_val())
	
	
	with open("output/figures/data/Fig1A.p", "wb") as f:
		pickle.dump([r_vals, b_vals, a_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Fig1A.p", "rb") as f:
		r_vals, b_vals, a_vals = pickle.load(f)

	fig,ax = plt.subplots(figsize=(3,2),dpi=300)
	ax.plot(r_vals,b_vals,color=fss.COLORS[0],linewidth=2)
	ax.set_xlabel("r")
	ax.set_ylim(-0.5,10)
	ax.set_ylabel("b(r)",color=fss.COLORS[0])
	
	ax2=ax.twinx()
	ax2.plot(r_vals, a_vals,color=fss.COLORS[1],linewidth=2)
	ax2.set_ylabel("a(r)",color=fss.COLORS[1])
	#ax2.set_ylim(0.05,1)
	plt.savefig("output/figures/plots/Fig1A.svg",transparent=True)
	plt.show()


if "B" in run_panels or "C" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	r_vals = [0.1,0.2,0.3]
	surv_vals = []
	rep_vals = []
	age_vals = []
	for r_val in r_vals:
		rep1.val = r_val
		surv_list = []
		rep_list = []
		age_list = list(np.linspace(0,80,num=1000))
		for age in age_list:
			surv_list.append(surv_func.ret_val_at_time(age))
			rep_list.append(rep_func.ret_val_at_time(age))
		surv_vals.append(surv_list)
		rep_vals.append(rep_list)
		age_vals.append(age_list)
	
	with open("output/figures/data/Fig1BC.p", "wb") as f:
		pickle.dump([r_vals, age_vals, rep_vals, surv_vals], f)

if "B" in plot_panels:
	with open("output/figures/data/Fig1BC.p", "rb") as f:
		r_vals, age_vals, rep_vals, surv_vals = pickle.load(f)
		
	plt.figure(figsize=(3,2), dpi=300)
	plt.plot(age_vals[0],surv_vals[0],linewidth=2,label="r = {}".format(r_vals[0]))
	plt.plot(age_vals[1],surv_vals[1],linewidth=2,label="r = {}".format(r_vals[1]))
	plt.plot(age_vals[2],surv_vals[2],linewidth=2,label="r = {}".format(r_vals[2]))
	plt.xlabel("x")
	plt.ylabel("l(x)")
	plt.legend()
	plt.savefig("output/figures/plots/Fig1B.svg",transparent=True)
	plt.show()

if "C" in plot_panels:
	with open("output/figures/data/Fig1BC.p", "rb") as f:
		r_vals, age_vals, rep_vals, surv_vals = pickle.load(f)
		
	plt.figure(figsize=(3,2), dpi=300)
	plt.plot(age_vals[0], rep_vals[0], linewidth=2, label="r = {}".format(r_vals[0]))
	plt.plot(age_vals[1], rep_vals[1], linewidth=2, label="r = {}".format(r_vals[1]))
	plt.plot(age_vals[2], rep_vals[2], linewidth=2, label="r = {}".format(r_vals[2]))
	plt.xlabel("x")
	plt.ylabel("m(x)")
	plt.legend()
	plt.savefig("output/figures/plots/Fig1C.svg",transparent=True)
	plt.show()


if "D" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.01,bfunc.rmax,num=100))
	g_vals = []
	
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Fig1D.p", "wb") as f:
		pickle.dump([r_vals,g_vals], f)
	
if "D" in plot_panels:
	with open("output/figures/data/Fig1D.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)

	fig = plt.figure(figsize=(3, 2), dpi=300)
	ax = fig.add_subplot()
	ax.plot(r_vals, g_vals, linewidth=2)
	ax.set_ylim(None,0.16)
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]

	ax.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
				xy=(arg_max, max_val), xycoords="data",
				xytext=(-30, 7), textcoords="offset points", zorder=2.5)
	plt.xlabel("r")
	plt.ylabel("g")
	plt.savefig("output/figures/plots/Fig1D.svg",transparent=True)
	plt.show()

if "E" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.01,bfunc.rmax,num=100))
	g_vals = []
	
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Fig1E.p", "wb") as f:
		pickle.dump([r_vals, g_vals], f)

if "E" in plot_panels:
	with open("output/figures/data/Fig1E.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)
	
	fig = plt.figure(figsize=(3, 2), dpi=300)
	ax = fig.add_subplot()
	ax.plot(r_vals, g_vals, linewidth=2)
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]
	
	ax.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
	            xy=(arg_max, max_val), xycoords="data",
	            xytext=(-40, 7), textcoords="offset points", zorder=2.5)
	plt.ylim(None,16)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.savefig("output/figures/plots/Fig1E.svg", transparent=True)
	plt.show()

