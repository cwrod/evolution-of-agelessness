import math
import sys

import scipy.optimize
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
	
	
	
gamma_vals = [0.0, 0.03, 0.05]
if "A" in run_panels:
	gamma_val = gamma_vals[0]
	rmax_vals = list(np.linspace(0.1, 0.9, num=10))
	beta_vals = list(np.linspace(0.1, 1.0, num=10))

	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp,extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	extr.gamma = gamma_val
	opt_r_mat = np.zeros((len(rmax_vals), len(beta_vals)))
	AMR_mat = np.zeros((len(rmax_vals), len(beta_vals)))
	
	for i in range(len(rmax_vals)):
		for j in range(len(beta_vals)):
			bfunc.rmax = rmax_vals[i]
			bfunc.beta = beta_vals[j]
			max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, rmax_vals[i]]])
			opt_r_mat[i, j] = max_rep_point[0]
			rep1.val = max_rep_point[0]
			AMR_mat[i, j] = gomp.get_AMR()
	
	with open("output/figures/data/Fig4A.p", "wb") as f:
		pickle.dump([gamma_val, rmax_vals, beta_vals, AMR_mat], f)

if "A" in plot_panels:
	with open("output/figures/data/Fig4A.p", "rb") as f:
		gamma_val, rmax_vals, beta_vals, AMR_mat = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	cur_mat = AMR_mat.transpose()
	centers = [min(rmax_vals), max(rmax_vals), min(beta_vals), max(beta_vals)]
	dx, = np.diff(centers[:2]) / (cur_mat.shape[1] - 1)
	dy, = -np.diff(centers[2:]) / (cur_mat.shape[0] - 1)
	extent = [centers[0] - dx / 2, centers[1] + dx / 2, centers[2] + dy / 2, centers[3] - dy / 2]
	pos = plt.imshow(cur_mat, cmap='hot', interpolation='nearest',
	                 extent=extent, origin="lower",
	                 aspect="auto", vmin=0, vmax=0.1)
	for x in range(len(rmax_vals)):
		for y in range(len(beta_vals)):
			if AMR_mat[x][y] == 0:
				plt.annotate('*', xycoords='data', xy=(rmax_vals[x], beta_vals[y]),
				             ha="center", va="center", color="white")
	plt.ylabel(r"$\beta$")
	plt.xlabel(r"$r_{max}$")
	cbar = plt.colorbar()
	cbar.set_label(r"$\beta × b\left(r_{opt}\right)$")
	plt.title(r"$\gamma$ = " + str(gamma_val))
	plt.savefig("output/figures/plots/Fig4A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
	gamma_val = gamma_vals[1]
	rmax_vals = list(np.linspace(0.1, 0.9, num=10))
	beta_vals = list(np.linspace(0.1, 1.0, num=10))
	
	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp,extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	extr.gamma = gamma_val
	opt_r_mat = np.zeros((len(rmax_vals), len(beta_vals)))
	AMR_mat = np.zeros((len(rmax_vals), len(beta_vals)))
	
	for i in range(len(rmax_vals)):
		for j in range(len(beta_vals)):
			bfunc.rmax = rmax_vals[i]
			bfunc.beta = beta_vals[j]
			max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, rmax_vals[i]]])
			opt_r_mat[i, j] = max_rep_point[0]
			rep1.val = max_rep_point[0]
			AMR_mat[i, j] = gomp.get_AMR()
	
	with open("output/figures/data/Fig4B.p", "wb") as f:
		pickle.dump([gamma_val, rmax_vals, beta_vals, AMR_mat], f)

if "B" in plot_panels:
	with open("output/figures/data/Fig4B.p", "rb") as f:
		gamma_val, rmax_vals, beta_vals, AMR_mat = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	cur_mat = AMR_mat.transpose()
	centers = [min(rmax_vals), max(rmax_vals), min(beta_vals), max(beta_vals)]
	dx, = np.diff(centers[:2]) / (cur_mat.shape[1] - 1)
	dy, = -np.diff(centers[2:]) / (cur_mat.shape[0] - 1)
	extent = [centers[0] - dx / 2, centers[1] + dx / 2, centers[2] + dy / 2, centers[3] - dy / 2]
	pos = plt.imshow(cur_mat, cmap='hot', interpolation='nearest',
	                 extent=extent, origin="lower",
	                 aspect="auto", vmin=0.0, vmax=0.1)
	for x in range(len(rmax_vals)):
		for y in range(len(beta_vals)):
			if AMR_mat[x][y] == 0:
				plt.annotate('*', xycoords='data', xy=(rmax_vals[x], beta_vals[y]),
				             ha="center", va="center", color="white")
	plt.ylabel(r"$\beta$")
	plt.xlabel(r"$r_{max}$")
	cbar = plt.colorbar()
	cbar.set_label(r"$\beta × b\left(r_{opt}\right)$")
	plt.title(r"$\gamma$ = " + str(gamma_val))
	plt.savefig("output/figures/plots/Fig4B.svg", transparent=True)
	plt.show()

if "C" in run_panels:
	gamma_val = gamma_vals[2]
	rmax_vals = list(np.linspace(0.1, 0.9, num=10))
	beta_vals = list(np.linspace(0.1, 1.0, num=10))

	rep1 = mf.Repair(0.32)
	bfunc = mf.InverseTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	extr.gamma = gamma_val
	opt_r_mat = np.zeros((len(rmax_vals), len(beta_vals)))
	AMR_mat = np.zeros((len(rmax_vals), len(beta_vals)))
	
	for i in range(len(rmax_vals)):
		for j in range(len(beta_vals)):
			bfunc.rmax = rmax_vals[i]
			bfunc.beta = beta_vals[j]
			max_growth, max_rep_point = mf.optimize_growth(growth, [rep1], [[0, rmax_vals[i]]])
			opt_r_mat[i, j] = max_rep_point[0]
			rep1.val = max_rep_point[0]
			AMR_mat[i, j] = gomp.get_AMR()
	
	with open("output/figures/data/Fig4C.p", "wb") as f:
		pickle.dump([gamma_val, rmax_vals, beta_vals, AMR_mat], f)

if "C" in plot_panels:
	with open("output/figures/data/Fig4C.p", "rb") as f:
		gamma_val, rmax_vals, beta_vals, AMR_mat = pickle.load(f)
	
	plt.figure(figsize=(4, 3), dpi=300)
	cur_mat = AMR_mat.transpose()
	centers = [min(rmax_vals), max(rmax_vals), min(beta_vals), max(beta_vals)]
	dx, = np.diff(centers[:2]) / (cur_mat.shape[1] - 1)
	dy, = -np.diff(centers[2:]) / (cur_mat.shape[0] - 1)
	extent = [centers[0] - dx / 2, centers[1] + dx / 2, centers[2] + dy / 2, centers[3] - dy / 2]
	pos = plt.imshow(cur_mat, cmap='hot', interpolation='nearest',
	                 extent=extent, origin="lower",
	                 aspect="auto", vmin=0.0, vmax=0.1)
	for x in range(len(rmax_vals)):
		for y in range(len(beta_vals)):
			if AMR_mat[x][y] == 0:
				plt.annotate('*', xycoords='data', xy=(rmax_vals[x], beta_vals[y]),
				             ha="center", va="center", color="white")
	plt.ylabel(r"$\beta$")
	plt.xlabel(r"$r_{max}$")
	cbar = plt.colorbar()
	cbar.set_label(r"$\beta × b\left(r_{opt}\right)$")
	plt.title(r"$\gamma$ = " + str(gamma_val))
	plt.savefig("output/figures/plots/Fig4C.svg", transparent=True)
	plt.show()
