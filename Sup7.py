import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

#run_panels = ["A","B","C","D","E","F"]
run_panels = []
#plot_panels = ["A","B","C","D","E","F"]
plot_panels = ["E"]


for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]

if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.6, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp,extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.58, bfunc.rmax, num=100))
	g_vals = []
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Sup7A.p", "wb") as f:
		pickle.dump([r_vals, g_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup7A.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(r_vals, g_vals, linewidth=2)
	
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]
	plt.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	offset = (-70, 5)
	plt.annotate("(" + str(arg_max) + "," + str(max_val) + ")",
				 xy=(arg_max, max_val), xycoords="data",
				 xytext=offset, textcoords="offset points", zorder=2.5)
	plt.ylim(None, 29.3)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.title(r"$r_{max}=0.6$")
	plt.locator_params(nbins=6)
	plt.savefig("output/figures/plots/Sup7A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.3, 0.2)
	bfunc2 = mf.InverseTOff(rep2, 0.3, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1,rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.295, 0.3, num=20))
	r2_vals = list(np.linspace(0.295, 0.3, num=20))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	
	with open("output/figures/data/Sup7B.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "B" in plot_panels:
	with open("output/figures/data/Sup7B.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	g_vals = np.transpose(g_vals)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.imshow(g_vals, cmap='coolwarm', extent=(min(r1_vals)/max(r1_vals), 1, min(r2_vals) / max(r2_vals), 1), origin='lower',
	           aspect='auto')
	plt.colorbar(label=r"$R_0$")
	plt.title(r"$r_{max,1}=0.3$    $r_{max,2}=0.3$")
	plt.xlabel(r"$r_1$ / $r_{max,1}$")
	plt.ylabel(r"$r_2$ / $r_{max,2}$")
	plt.locator_params(nbins=5)
	plt.savefig("output/figures/plots/Sup7B.svg", transparent=True)
	plt.show()

if "C" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.05, 0.2)
	bfunc2 = mf.InverseTOff(rep2, 0.55, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.0495, 0.05, num=20))
	r2_vals = list(np.linspace(0.53, 0.55, num=20))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	with open("output/figures/data/Sup7C.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "C" in plot_panels:
	with open("output/figures/data/Sup7C.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
		
	g_vals = np.transpose(g_vals)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.imshow(g_vals, cmap='coolwarm', extent=(min(r1_vals)/max(r1_vals), 1, min(r2_vals) / max(r2_vals), 1), origin='lower',
	           aspect='auto')
	plt.colorbar(label=r"$R_0$")
	plt.title(r"$r_{max,1}=0.05$      $r_{max,2}=0.55$")
	plt.xlabel(r"$r_1$ / $r_{max,1}$")
	plt.ylabel(r"$r_2$ / $r_{max,2}$")
	plt.locator_params(nbins=5)
	plt.savefig("output/figures/plots/Sup7C.svg", transparent=True)
	plt.show()

if "D" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.7, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.68, bfunc.rmax, num=100))
	g_vals = []
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Sup7D.p", "wb") as f:
		pickle.dump([r_vals, g_vals], f)

if "D" in plot_panels:
	with open("output/figures/data/Sup7D.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(r_vals, g_vals, linewidth=2)
	
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]
	plt.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	offset = (-20, 10)
	plt.annotate("(" + str(arg_max) + "," + str(max_val) + ")",
	             xy=(arg_max, max_val), xycoords="data",
	             xytext=offset, textcoords="offset points", zorder=2.5)
	plt.ylim(None, 20.4)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.title(r"$r_{max}=0.7$")
	plt.locator_params(nbins=6)
	plt.savefig("output/figures/plots/Sup7D.svg", transparent=True)
	plt.show()

if "E" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.35, 0.2)
	bfunc2 = mf.InverseTOff(rep2, 0.35, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.34, 0.35, num=20))
	r2_vals = list(np.linspace(0.34, 0.35, num=20))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	with open("output/figures/data/Sup7E.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "E" in plot_panels:
	with open("output/figures/data/Sup7E.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	g_vals = np.transpose(g_vals)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.imshow(g_vals, cmap='coolwarm', extent=(min(r1_vals) / max(r1_vals), 1, min(r2_vals) / max(r2_vals), 1),
	           origin='lower',
	           aspect='auto')
	plt.colorbar(label=r"$R_0$")
	plt.title(r"$r_{max,1}=0.35$    $r_{max,2}=0.35$")
	plt.xlabel(r"$r_1$ / $r_{max,1}$")
	plt.ylabel(r"$r_2$ / $r_{max,2}$")
	
	#plt.locator_params(nbins=6)
	plt.savefig("output/figures/plots/Sup7E.svg", transparent=True)
	plt.show()

if "F" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.05, 0.2)
	bfunc2 = mf.InverseTOff(rep2, 0.65, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.048, 0.05, num=20))
	r2_vals = list(np.linspace(0.6, 0.65, num=20))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	with open("output/figures/data/Sup7F.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "F" in plot_panels:
	with open("output/figures/data/Sup7F.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	g_vals = np.transpose(g_vals)
	
	plt.figure(figsize=(4, 3), dpi=300)
	plt.imshow(g_vals, cmap='coolwarm', extent=(min(r1_vals) / max(r1_vals), 1, min(r2_vals) / max(r2_vals), 1),
	           origin='lower',
	           aspect='auto')
	plt.colorbar(label=r"$R_0$")
	plt.title(r"$r_{max,1}=0.05$    $r_{max,2}=0.65$")
	plt.xlabel(r"$r_1$ / $r_{max,1}$")
	plt.ylabel(r"$r_2$ / $r_{max,2}$")
	plt.locator_params(nbins=6)
	plt.savefig("output/figures/plots/Sup7F.svg", transparent=True)
	plt.show()
