import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

run_panels = []
plot_panels = ["A","B","C","D"]


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
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	rmax_vals = [0.1, 0.5, 0.9]
	r_vals_list = []
	b_vals_list = []
	for rmax_val in rmax_vals:
		b_vals = []
		bfunc.rmax = rmax_val
		r_vals = list(np.linspace(0.001, bfunc.rmax, num=100))
		for r_val in r_vals:
			rep1.val = r_val
			b_vals.append(bfunc.ret_val())
		r_vals_list.append(r_vals)
		b_vals_list.append(b_vals)
		
	with open("output/figures/data/Fig3A.p", "wb") as f:
		pickle.dump([rmax_vals, r_vals_list, b_vals_list], f)

if "A" in plot_panels:
	with open("output/figures/data/Fig3A.p", "rb") as f:
		rmax_vals, r_vals_list, b_vals_list = pickle.load(f)
	
	fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
	for i in range(len(r_vals_list)):
		ax.plot(r_vals_list[i], b_vals_list[i], color=fss.COLORS[i], linewidth=2, label=r"$r_{max}$ = " + str(rmax_vals[i]))
		
	ax.set_xlabel("r")
	ax.set_ylim(-0.02, 0.5)
	ax.set_ylabel(r"$\beta * b(r)$")
	
	plt.legend()
	plt.savefig("output/figures/plots/Fig3A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
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
	
	rmax_vals = [0.1, 0.5, 0.9]
	r_vals_list = []
	g_vals_list = []
	for rmax_val in rmax_vals:
		r_vals = list(np.linspace(0.01, rmax_val, num=100))
		g_vals = []
		bfunc.rmax = rmax_val
		for r_val in r_vals:
			rep1.val = r_val
			g_vals.append(growth.ret_val())
		r_vals_list.append(r_vals)
		g_vals_list.append(g_vals)
		
	with open("output/figures/data/Fig3B.p", "wb") as f:
		pickle.dump([rmax_vals, r_vals_list, g_vals_list], f)

if "B" in plot_panels:
	with open("output/figures/data/Fig3B.p", "rb") as f:
		rmax_vals, r_vals_list, g_vals_list = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	for i in range(len(r_vals_list)):
		plt.plot(r_vals_list[i], g_vals_list[i], linewidth=2, label=r"$r_{max}$ = " + str(rmax_vals[i]))
		
		max_val = max(g_vals_list[i])
		arg_max = r_vals_list[i][g_vals_list[i].index(max_val)]
		plt.scatter(arg_max, max_val, color="black", zorder=2.5)
		arg_max = round(arg_max, 2)
		max_val = round(max_val, 3)
		offset = [(-10, 7), (-40, 7), (-30, 7)][i]
		plt.annotate("(" + str(arg_max) + "," + str(max_val) + ")",
		             xy=(arg_max, max_val), xycoords="data",
		             xytext=offset, textcoords="offset points", zorder=2.5)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.ylim(None, 30)
	plt.legend()
	plt.savefig("output/figures/plots/Fig3B.svg", transparent=True)
	plt.show()

if "C" in run_panels or "D" in run_panels:
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
	
	rmax_vals = list(np.linspace(0.02, 1.0, num=100))
	opt_rs = []
	AMRs = []
	for rmax_val in rmax_vals:
		bfunc.rmax = rmax_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, rmax_val]])
		opt_rs.append(max_rep_point[0])
		rep1.val = max_rep_point[0]
		AMRs.append(gomp.get_AMR())
		rep1.val = bfunc.rmax
	
	with open("output/figures/data/Fig3CD.p", "wb") as f:
		pickle.dump([rmax_vals, opt_rs, AMRs], f)

if "C" in plot_panels:
	with open("output/figures/data/Fig3CD.p", "rb") as f:
		rmax_vals, opt_rs, AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot([0,1],[0,1], linestyle="dashed", color="black", linewidth=0.9)
	plt.plot(rmax_vals, opt_rs, linewidth=2)
	plt.xlabel(r"$r_{max}$")
	plt.ylabel(r"$r_{opt}$")
	plt.savefig("output/figures/plots/Fig3C.svg", transparent=True)
	plt.show()

if "D" in plot_panels:
	with open("output/figures/data/Fig3CD.p", "rb") as f:
		rmax_vals, opt_rs, AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(rmax_vals, AMRs, linewidth=2)
	plt.xlabel(r"$r_{max}$")
	plt.ylabel(r"$\beta * b\left(r_{opt}\right)$")
	plt.savefig("output/figures/plots/Fig3D.svg", transparent=True)
	plt.show()
