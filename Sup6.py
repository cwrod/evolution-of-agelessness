import math
import sys
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf

#run_panels = ["A","B","C","D"]
run_panels = []
plot_panels = ["A","B","C","D"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]
	

if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp1 = mf.GompertzRepAMR(0.0006, bfunc1)
	gomp2 = mf.GompertzCons(0.0006, 0.2)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	AMR2_vals = [0,0.2,0.4]
	r_vals = list(np.linspace(0.01, bfunc1.rmax, num=100))
	g_vals = []
	for AMR2_val in AMR2_vals:
		gomp2.B = AMR2_val
		g_val_list = []
		for r_val in r_vals:
			rep1.val = r_val
			g_val_list.append(growth.ret_val())
		g_vals.append(g_val_list)

	with open("output/figures/data/Sup6A.p", "wb") as f:
		pickle.dump([AMR2_vals, r_vals, g_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup6A.p", "rb") as f:
		AMR2_vals, r_vals, g_vals = pickle.load(f)

	fig = plt.figure(figsize=(3.5, 2.5), dpi=300)
	ax = fig.add_subplot()
	for i,AMR2_val in enumerate(AMR2_vals):
		ax.plot(r_vals, g_vals[i], linewidth=2, label=r"$AMR_2="+str(AMR2_val)+"$")
		max_val = max(g_vals[i])
		arg_max = r_vals[g_vals[i].index(max_val)]

		ax.scatter(arg_max, max_val, color="black", zorder=2.5)
		arg_max = round(arg_max, 2)
		max_val = round(max_val, 1)
		xy_poses = [(-15,6),(5,5),(4,2)]
		ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
					xy=(arg_max, max_val), xycoords="data", fontsize=8,
					xytext=xy_poses[i], textcoords="offset points", zorder=2.5)
	plt.ylim(0,17)
	plt.xlim(0,0.6)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.legend(loc="upper left",prop={"size":7})
	plt.savefig("output/figures/plots/Sup6A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp1 = mf.GompertzRepAMR(0.0006, bfunc1)
	gomp2 = mf.GompertzCons(0.0006, 0.2)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	AMR2_vals = list(np.linspace(0.0, 0.5, num=50))
	AMRs = []
	for AMR2_val in AMR2_vals:
		gomp2.B = AMR2_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc1.rmax]])
		rep1.val = max_rep_point[0]
		AMRs.append(gomp1.get_AMR())
	
	with open("output/figures/data/Sup6B.p", "wb") as f:
		pickle.dump([AMR2_vals, AMRs], f)

if "B" in plot_panels:
	with open("output/figures/data/Sup6B.p", "rb") as f:
		AMR2_vals, AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(AMR2_vals, AMRs, linewidth=2)
	plt.xlabel(r"$AMR_2$")
	plt.ylabel(r"$\beta × b\left(r_{opt}\right)$")
	plt.savefig("output/figures/plots/Sup6B.svg", transparent=True)
	plt.show()

if "C" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScalingDrain(rep1, 0.0)
	gomp1 = mf.GompertzRepAMR(0.0006, bfunc1)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	d_vals = [0.0,0.2,0.4]
	r_vals = list(np.linspace(0.01, bfunc1.rmax, num=100))
	g_vals = []
	for d_val in d_vals:
		hfunc.d = d_val
		g_val_list = []
		for r_val in r_vals:
			rep1.val = r_val
			g_val_list.append(growth.ret_val())
		g_vals.append(g_val_list)
	
	with open("output/figures/data/Sup6C.p", "wb") as f:
		pickle.dump([d_vals, r_vals, g_vals], f)

if "C" in plot_panels:
	with open("output/figures/data/Sup6C.p", "rb") as f:
		d_vals, r_vals, g_vals = pickle.load(f)
	
	fig = plt.figure(figsize=(3.5, 2.5), dpi=300)
	ax = fig.add_subplot()
	for i, d_val in enumerate(d_vals):
		ax.plot(r_vals, g_vals[i], linewidth=2, label=r"$d=" + str(d_val) + "$")
		max_val = max(g_vals[i])
		arg_max = r_vals[g_vals[i].index(max_val)]
		
		ax.scatter(arg_max, max_val, color="black", zorder=2.5)
		arg_max = round(arg_max, 2)
		max_val = round(max_val, 1)
		xy_poses = [(-15, 7), (-15, 7), (-15, 7)]
		ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
		            xy=(arg_max, max_val), xycoords="data", fontsize=8,
		            xytext=xy_poses[i], textcoords="offset points", zorder=2.5)
	plt.ylim(0, 17)
	plt.xlim(0, 0.6)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.legend(loc="upper left", prop={"size": 7})
	plt.savefig("output/figures/plots/Sup6C.svg", transparent=True)
	plt.show()

if "D" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScalingDrain(rep1, 0.0)
	gomp1 = mf.GompertzRepAMR(0.0006, bfunc1)
	gomp2 = mf.GompertzCons(0.0006, 0.2)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	d_vals = list(np.linspace(0.0, 0.4, num=50))
	AMRs = []
	for d_val in d_vals:
		hfunc.d = d_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc1.rmax]])
		rep1.val = max_rep_point[0]
		AMRs.append(gomp1.get_AMR())
	
	with open("output/figures/data/Sup6D.p", "wb") as f:
		pickle.dump([d_vals, AMRs], f)

if "D" in plot_panels:
	with open("output/figures/data/Sup6D.p", "rb") as f:
		d_vals, AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(d_vals, AMRs, linewidth=2)
	plt.xlabel(r"$d$")
	plt.ylabel(r"$\beta × b\left(r_{opt}\right)$")
	plt.savefig("output/figures/plots/Sup6D.svg", transparent=True)
	plt.show()
