import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

run_panels = ["A"]
plot_panels = ["A"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]

if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.01)
	hfunc = mf.ReproductiveScaling(rep1)
	arm_doll = mf.ArmitageDoll(bfunc, 1e6, 6)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([arm_doll], mat_age, prv)
	surv_func = mf.SurvFunc([arm_doll, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	opt_rs = []
	opt_brs = []
	opt_scales = []
	C_sizes = list(np.logspace(1, 25, num=100))

	for C_size in C_sizes:
		arm_doll.C_size = C_size
		max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, 0.8]])
		opt_rs.append(max_rep_point[0])
		rep1.val = max_rep_point[0]
		opt_brs.append(bfunc.ret_val())
		opt_scales.append(arm_doll.risk_scale())
	
	with open("output/figures/data/Sup2A.p", "wb") as f:
		pickle.dump([C_sizes, opt_rs, opt_brs, opt_scales], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup2A.p", "rb") as f:
		C_sizes, opt_rs, opt_brs, opt_scales = pickle.load(f)
	
	fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
	line1 = ax.plot(C_sizes, opt_scales, linewidth=2)
	ax.set_xlabel("C")
	ax.set_xscale("log")
	ax.set_ylabel(r"$\frac{C × b\left(r_{opt}\right)^k}{(k-1)!}$")
	plt.savefig("output/figures/plots/Sup2A.svg", transparent=True)
	plt.show()
	
	#fig, ax = plt.subplots(figsize=(4, 3), dpi=300)
	#line1 = ax.plot(C_sizes, opt_brs, linewidth=2)
	#ax.set_xlabel("C")
	#ax.set_xscale("log")
	#ax.set_ylabel(r"$b(r_{opt})$")
	#plt.savefig("output/figures/plots/Sup7A.svg", transparent=True)
	#plt.show()
