import math
import sys
import matplotlib
import matplotlib.colors
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pickle
import operator
import model_funcs as mf
import os
import pandas as pd
import datetime

run_panels = ["A","B"]
plot_panels = ["A","B"]


MODE = "DEBUG"

for arg in sys.argv[1:]:
	if arg == "-t":
		MODE = "TEST"
	elif arg == "-d":
		MODE = "DEBUG"
	else:
		run_panels = [x for x in arg.strip()]
		plot_panels = [x for x in arg.strip()]

if MODE == "DEBUG":
	mf.MODE = "DEBUG"
else:
	mf.MODE = "TEST"

if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.25, 0.2)
	bfunc2 = mf.InverseTOff(rep2, 0.25, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
	gomp1 = mf.GompertzRepAMR(0.0003, bfunc1)
	gomp2 = mf.GompertzRepAMR(0.0003, bfunc2)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp1, gomp2], mat_age, prv)
	surv_func = mf.SurvFunc([gomp1, gomp2, extr])
	growth = mf.KSelPop(rep_func, surv_func)

	beta1_vals = np.linspace(0.2, 1.5, num=10)
	beta2_vals = np.linspace(0.2, 1.5, num=10)
	AMR1_vals = np.zeros((len(beta1_vals), len(beta2_vals)))
	AMR2_vals = np.zeros((len(beta1_vals), len(beta2_vals)))
	for i, beta1_val in enumerate(beta1_vals):
		bfunc1.beta = beta1_val
		for j, beta2_val in enumerate(beta2_vals):
			bfunc2.beta = beta2_val
			res = mf.optimize_growth(growth, [rep1, rep2], [[0, bfunc1.rmax], [0, bfunc2.rmax]])
			
			rep1.val = res[1][0]
			rep2.val = res[1][1]
			AMR1_vals[i][j] = gomp1.get_AMR()
			AMR2_vals[i][j] = gomp2.get_AMR()
	
	with open("output/figures/data/Fig5A.p", "wb") as f:
		pickle.dump([beta1_vals, beta2_vals, AMR1_vals, AMR2_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Fig5A.p", "rb") as f:
		beta1_vals, beta2_vals, AMR1_vals, AMR2_vals = pickle.load(f)
	AMR1_vals = AMR1_vals.transpose()
	AMR2_vals = AMR2_vals.transpose()
	
	
	def getCustomSymbol1(path_index=1):
		if path_index == 1:
			verts = [
				(-0.5, -0.5),
				(0.5, -0.5),
				(0.5, 0.5),
				(-0.5, -0.5), ]
		else:
			verts = [
				(-0.5, -0.5),
				(-0.5, 0.5),
				(0.5, 0.5),
				(-0.5, -0.5), ]
		codes = [matplotlib.path.Path.MOVETO,
		         matplotlib.path.Path.LINETO,
		         matplotlib.path.Path.LINETO,
		         matplotlib.path.Path.CLOSEPOLY,
		         ]
		pathCS1 = matplotlib.path.Path(verts, codes)
		return pathCS1, verts
	
	
	def plot_mat(matrix, path_index=1, alpha=1.0, vmin=0., vmax=1.):
		nx, ny = matrix.shape
		X, Y, values = zip(*[(i, j, matrix[i, j]) for i in range(nx) for j in range(ny)])
		marker, verts = getCustomSymbol1(path_index=path_index)
		X = [beta1_vals[x] for x in X]
		Y = [beta2_vals[y] for y in Y]
		ax.scatter(X, Y, s=650,
		           marker=marker,
		           c=values,
		           cmap='hot',
		           alpha=alpha,
		           vmin=vmin, vmax=vmax)
		
		for i, val in enumerate(values):
			if val < 1e-5:
				if path_index == 1:
					dx = 1.5
					dy = -1.5
				else:
					dx = -1.5
					dy = 1.5
				star_marker = matplotlib.markers.MarkerStyle('*', fillstyle='none')
				star_marker._transform = star_marker.get_transform().translate(dx, dy)
				ax.scatter(X[i], Y[i], s=10, marker=star_marker, c="white")
		return
	
	
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	A = np.random.uniform(20, 50, 30).reshape([6, 5])
	B = np.random.uniform(40, 70, 30).reshape([6, 5])
	vmin = np.min([AMR1_vals, AMR2_vals])
	vmax = np.max([AMR1_vals, AMR2_vals]) + 0.001
	plot_mat(path_index=2, vmin=vmin, vmax=vmax, matrix=AMR1_vals)
	plot_mat(path_index=1, vmin=vmin, vmax=vmax, matrix=AMR2_vals)
	plt.xlim([0.12, 1.58])
	plt.ylim([0.12, 1.58])
	sm = plt.cm.ScalarMappable(cmap='hot', norm=plt.Normalize(vmin=vmin, vmax=vmax))
	sm._A = []
	cbar = plt.colorbar(sm)
	cbar.set_label(r"$\beta * b(r_{opt})$",size=16)
	plt.ylabel(r"$\beta_1$",size=16)
	plt.xlabel(r"$\beta_2$",size=16)
	cbar.ax.tick_params(labelsize=12)
	ax.tick_params(axis='both', which='major', labelsize=12)
	ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
	ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
	plt.savefig("output/figures/plots/Fig5A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.InverseTOff(rep1, 0.4, 0.2)
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
	
	rmax1_vals = np.linspace(0.05, 0.2, num=10)
	rmax2_vals = np.linspace(0.05, 0.2, num=10)
	AMR1_vals = np.zeros((len(rmax1_vals), len(rmax2_vals)))
	AMR2_vals = np.zeros((len(rmax1_vals), len(rmax2_vals)))
	for i, rmax1_val in enumerate(rmax1_vals):
		bfunc1.rmax = rmax1_val
		for j, rmax2_val in enumerate(rmax2_vals):
			bfunc2.rmax = rmax2_val
			res = mf.optimize_growth(growth, [rep1, rep2], [[0, bfunc1.rmax], [0, bfunc2.rmax]])
			
			rep1.val = res[1][0]
			rep2.val = res[1][1]
			AMR1_vals[i][j] = gomp1.get_AMR()
			AMR2_vals[i][j] = gomp2.get_AMR()
	
	with open("output/figures/data/Fig5B.p", "wb") as f:
		pickle.dump([rmax1_vals, rmax2_vals, AMR1_vals, AMR2_vals], f)

if "B" in plot_panels:
	with open("output/figures/data/Fig5B.p", "rb") as f:
		rmax1_vals, rmax2_vals, AMR1_vals, AMR2_vals = pickle.load(f)
	AMR1_vals = AMR1_vals.transpose()
	AMR2_vals = AMR2_vals.transpose()
	
	
	def getCustomSymbol1(path_index=1):
		if path_index == 1:
			verts = [
				(-0.5, -0.5),
				(0.5, -0.5),
				(0.5, 0.5),
				(-0.5, -0.5), ]
		else:
			verts = [
				(-0.5, -0.5),
				(-0.5, 0.5),
				(0.5, 0.5),
				(-0.5, -0.5), ]
		codes = [matplotlib.path.Path.MOVETO,
		         matplotlib.path.Path.LINETO,
		         matplotlib.path.Path.LINETO,
		         matplotlib.path.Path.CLOSEPOLY,
		         ]
		pathCS1 = matplotlib.path.Path(verts, codes)
		return pathCS1, verts
	
	
	def plot_mat(matrix, path_index=1, alpha=1.0, vmin=0., vmax=1.):
		nx, ny = matrix.shape
		X, Y, values = zip(*[(i, j, matrix[i, j]) for i in range(nx) for j in range(ny)])
		marker, verts = getCustomSymbol1(path_index=path_index)
		X = [rmax1_vals[x] for x in X]
		Y = [rmax2_vals[y] for y in Y]
		ax.scatter(X, Y, s=650,
		           marker=marker,
		           c=values,
		           cmap='hot',
		           alpha=alpha,
		           vmin=vmin, vmax=vmax)
		
		for i,val in enumerate(values):
			if val < 1e-5:
				print(i)
				if path_index == 1:
					dx = 1.5
					dy = -1.5
				else:
					dx = -1.5
					dy = 1.5
				star_marker = matplotlib.markers.MarkerStyle('*', fillstyle='none')
				star_marker._transform = star_marker.get_transform().translate(dx,dy)
				ax.scatter(X[i], Y[i], s=10, marker=star_marker,c="white")
		return
	
	
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	A = np.random.uniform(20, 50, 30).reshape([6, 5])
	B = np.random.uniform(40, 70, 30).reshape([6, 5])
	vmin = np.min([AMR1_vals, AMR2_vals])
	vmax = np.max([AMR1_vals, AMR2_vals]) + 0.001
	plot_mat(path_index=2, vmin=vmin, vmax=vmax, matrix=AMR1_vals)
	plot_mat(path_index=1, vmin=vmin, vmax=vmax, matrix=AMR2_vals)
	plt.xlim([0.0401, 0.209])
	plt.ylim([0.0401, 0.209])
	sm = plt.cm.ScalarMappable(cmap='hot', norm=plt.Normalize(vmin=vmin, vmax=vmax))
	sm._A = []
	cbar = plt.colorbar(sm)
	cbar.set_label(r"$\beta * b(r_{opt})$",size=16)
	plt.ylabel(r"$r_{max,1}$",size=16)
	plt.xlabel(r"$r_{max,2}$",size=16)
	cbar.ax.tick_params(labelsize=12)
	ax.tick_params(axis='both', which='major', labelsize=12)
	ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
	ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
	plt.savefig("output/figures/plots/Fig5B.svg", transparent=True)
	plt.show()

