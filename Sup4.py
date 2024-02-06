import math
import sys
import scipy.optimize
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

run_panels = []
plot_panels = ["A","B"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]
	

if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.0)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc.rmax]])
	
	gamma_vals = list(np.linspace(0, max_growth-0.005, num=20))

	R_AMRs = []
	for gamma_val in gamma_vals:
		extr.gamma = gamma_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc.rmax]])
		rep1.val = max_rep_point[0]
		R_AMRs.append(gomp.get_AMR())
	
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.0)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	K_AMRs = []
	for gamma_val in gamma_vals:
		extr.gamma = gamma_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc.rmax]])
		rep1.val = max_rep_point[0]
		K_AMRs.append(gomp.get_AMR())

	with open("output/figures/data/Sup4A.p", "wb") as f:
		pickle.dump([gamma_vals, R_AMRs, K_AMRs], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup4A.p", "rb") as f:
		gamma_vals, R_AMRs, K_AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3),dpi=300)
	plt.plot(gamma_vals, K_AMRs, linewidth=2, color=fss.COLORS[0], label=r"$K$-selection")
	plt.plot(gamma_vals, R_AMRs, linewidth=2, color=fss.COLORS[1], label=r"$r$-selection")
	plt.xlabel(r"$\gamma$")
	plt.ylabel(r"$\beta * b\left(r_{opt}\right)$")
	plt.legend()
	plt.savefig("output/figures/plots/Sup4A.svg", transparent=True)
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
	growth = mf.RSelPop(rep_func, surv_func)
	
	rep1.val = bfunc.rmax
	min_kappa_possible = 1 / mf.integrate(lambda x: surv_func.ret_val_at_time(x) * rep_func.ret_val_at_time(x),
	                                      0, surv_func.max_age())
	
	kappa_vals = list(np.linspace(min_kappa_possible + 0.001, 1.0, num=20))

	R_AMRs = []
	for kappa_val in kappa_vals:
		prv.scale = kappa_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc.rmax]])
		rep1.val = max_rep_point[0]
		R_AMRs.append(gomp.get_AMR())
	
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
	
	K_AMRs = []
	for kappa_val in kappa_vals:
		prv.scale = kappa_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, bfunc.rmax]])
		rep1.val = max_rep_point[0]
		K_AMRs.append(gomp.get_AMR())

	with open("output/figures/data/Sup4B.p", "wb") as f:
		pickle.dump([kappa_vals, R_AMRs, K_AMRs], f)

if "B" in plot_panels:
	with open("output/figures/data/Sup4B.p", "rb") as f:
		kappa_vals, R_AMRs, K_AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3),dpi=300)
	plt.plot(kappa_vals, K_AMRs, linewidth=2, color=fss.COLORS[0], label=r"$K$-selection")
	plt.plot(kappa_vals, R_AMRs, linewidth=2, color=fss.COLORS[1], label=r"$r$-selection")
	plt.xlabel(r"$\kappa$")
	plt.ylabel(r"$\beta * b\left(r_{opt}\right)$")
	plt.legend()
	plt.savefig("output/figures/plots/Sup4B.svg", transparent=True)
	plt.show()
