import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf

run_panels = []
plot_panels = ["A","B","C"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]

# Mat age
if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(40, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)

	rmax_vals = list(np.linspace(0.02, 0.95, num=100))
	opt_rs = []
	AMRs = []
	for rmax_val in rmax_vals:
		bfunc.rmax = rmax_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, rmax_val]])
		opt_rs.append(max_rep_point[0])
		rep1.val = max_rep_point[0]
		AMRs.append(gomp.get_AMR())
	
	with open("output/figures/data/Sup3A.p", "wb") as f:
		pickle.dump([rmax_vals, opt_rs, AMRs], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup3A.p", "rb") as f:
		rmax_vals, opt_rs, AMRs = pickle.load(f)
		
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(rmax_vals, AMRs, linewidth=2)
	plt.xlabel(r"$r_{max}$")
	plt.ylabel(r"$\beta * b\left(r_{opt}\right)$")
	plt.savefig("output/figures/plots/Sup3A.svg", transparent=True)
	plt.show()


if "B" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.ScalingTOff(rep1, 0.5, 0.8)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.01)
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

	with open("output/figures/data/Sup3B.p", "wb") as f:
		pickle.dump([rmax_vals, opt_rs, AMRs], f)

if "B" in plot_panels:
	with open("output/figures/data/Sup3B.p", "rb") as f:
		rmax_vals, opt_rs, AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(rmax_vals, AMRs, linewidth=2)
	plt.xlabel(r"$r_{max}$")
	plt.ylabel(r"$\beta * b\left(r_{opt}\right)$")
	plt.savefig("output/figures/plots/Sup3B.svg", transparent=True)
	plt.show()


# r-selection with high beta
if "C" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.6)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 4)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)

	rmax_vals = list(np.linspace(0.02, 1.0, num=20))
	opt_rs = []
	AMRs = []
	for rmax_val in rmax_vals:
		bfunc.rmax = rmax_val
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, rmax_val]])
		opt_rs.append(max_rep_point[0])
		rep1.val = max_rep_point[0]
		AMRs.append(gomp.get_AMR())
	
	with open("output/figures/data/Sup3C.p", "wb") as f:
		pickle.dump([rmax_vals, opt_rs, AMRs], f)

if "C" in plot_panels:
	with open("output/figures/data/Sup3C.p", "rb") as f:
		rmax_vals, opt_rs, AMRs = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(rmax_vals, AMRs, linewidth=2)
	plt.xlabel(r"$r_{max}$")
	plt.ylabel(r"$\beta * b\left(r_{opt}\right)$")
	plt.savefig("output/figures/plots/Sup3C.svg", transparent=True)
	plt.show()
