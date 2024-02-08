import math
import sys
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf

run_panels = ["A","B","C"]
plot_panels = ["A","B","C"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]


if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.4, 0.2)
	bfunc2 = mf.InverseTOff(rep2, 0.4, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1,rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.0, bfunc1.rmax, num=20))
	r2_vals = list(np.linspace(0.0, bfunc2.rmax, num=20))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	r1_vals, r2_vals = np.meshgrid(r1_vals, r2_vals)
	with open("output/figures/data/Sup5A.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup5A.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, dpi=300)
	surf = ax.plot_surface(r1_vals, r2_vals, g_vals, cmap="coolwarm")
	ax.zaxis.set_rotate_label(False)
	ax.set_ylabel(r"$r_1$")
	ax.set_xlabel(r"$r_2$")
	ax.set_zlabel(r"$R_0$")
	ax.view_init(60, -105)
	plt.savefig("output/figures/plots/Sup5A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.SharpTOff(rep1, 0.49, 0.7)
	bfunc2 = mf.SharpTOff(rep2, 0.49, 0.7)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.1)
	prv = mf.PeakReproductiveValue(hfunc, 50)
	mat_age = mf.ConstMatAge(10)
	rep_func = mf.ConstRepro(mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.0, 0.49, num=50))
	r2_vals = list(np.linspace(0.0, 0.49, num=50))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	r1_vals, r2_vals = np.meshgrid(r1_vals, r2_vals)
	with open("output/figures/data/Sup5B.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "B" in plot_panels:
	with open("output/figures/data/Sup5B.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, dpi=300)
	surf = ax.plot_surface(r1_vals, r2_vals, g_vals, cmap="coolwarm")
	ax.zaxis.set_rotate_label(False)
	ax.set_ylabel(r"$r_1$")
	ax.set_xlabel(r"$r_2$")
	ax.set_zlabel(r"$R_0$")
	ax.view_init(60, -105)
	plt.savefig("output/figures/plots/Sup5B.svg", transparent=True)
	plt.show()

if "C" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.4, 0.02)
	bfunc2 = mf.InverseTOff(rep2, 0.4, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r1_vals = list(np.linspace(0.0, 0.4, num=20))
	r2_vals = list(np.linspace(0.0, 0.4, num=20))
	g_vals = np.zeros((len(r1_vals), len(r2_vals)))
	for i, r1_val in enumerate(r1_vals):
		for j, r2_val in enumerate(r2_vals):
			rep1.val = r1_val
			rep2.val = r2_val
			g_vals[i][j] = growth.ret_val()
	
	r1_vals, r2_vals = np.meshgrid(r1_vals, r2_vals)
	with open("output/figures/data/Sup5C.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "C" in plot_panels:
	with open("output/figures/data/Sup5C.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, dpi=300)
	surf = ax.plot_surface(r1_vals, r2_vals, g_vals, cmap="coolwarm")
	ax.zaxis.set_rotate_label(False)
	ax.set_ylabel(r"$r_1$")
	ax.set_xlabel(r"$r_2$")
	ax.set_zlabel(r"$R_0$")
	ax.view_init(60, -105)
	plt.savefig("output/figures/plots/Sup5C.svg", transparent=True)
	plt.show()
