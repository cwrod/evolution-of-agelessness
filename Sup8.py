import math
import sys
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss
import pandas as pd
import datetime
import os

#run_panels = ["A","B","C"]
run_panels = []
#plot_panels = ["A","B","C"]
plot_panels = ["C"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]

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
			rep1.val = bfunc1.rmax
			res = mf.optimize_growth(growth, [rep2], [[0, bfunc2.rmax]])
			
			rep2.val = res[1][0]
			AMR1_vals[i][j] = gomp1.get_AMR()
			AMR2_vals[i][j] = gomp2.get_AMR()
	
	with open("output/figures/data/Sup8A.p", "wb") as f:
		pickle.dump([beta1_vals, beta2_vals, AMR1_vals, AMR2_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup8A.p", "rb") as f:
		beta1_vals, beta2_vals, AMR1_vals, AMR2_vals = pickle.load(f)
	AMR1_vals = AMR1_vals.transpose()
	AMR2_vals = AMR2_vals.transpose()
	
	with open("output/figures/data/Fig5A.p", "rb") as f:
		fig5_beta1_vals, fig5_beta2_vals, fig5_AMR1_vals, fig5_AMR2_vals = pickle.load(f)
	fig5_AMR1_vals = fig5_AMR1_vals.transpose()
	fig5_AMR2_vals = fig5_AMR2_vals.transpose()
	
	
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
				
				fig5_i = list(beta1_vals).index(list(X)[i])
				fig5_j = list(beta2_vals).index(list(Y)[i])
				if path_index == 1:
					is_ageless = (fig5_AMR1_vals[fig5_i, fig5_j] < 1e-5)
				else:
					is_ageless = (fig5_AMR2_vals[fig5_i, fig5_j] < 1e-5)

				star_color = "cyan" if not is_ageless else "white"
				ax.scatter(X[i], Y[i], s=10, marker=star_marker, c=star_color)
	
	
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	vmin = np.min([AMR1_vals, AMR2_vals])
	vmax = np.max([AMR1_vals, AMR2_vals]) + 0.001
	plot_mat(path_index=2, vmin=vmin, vmax=vmax, matrix=AMR1_vals)
	plot_mat(path_index=1, vmin=vmin, vmax=vmax, matrix=AMR2_vals)
	plt.xlim([0.12, 1.58])
	plt.ylim([0.12, 1.58])
	sm = plt.cm.ScalarMappable(cmap='hot', norm=plt.Normalize(vmin=vmin, vmax=vmax))
	sm._A = []
	cbar = plt.colorbar(sm)
	cbar.set_label(r"$\beta × b(r_{opt})$", size=16)
	plt.ylabel(r"$\beta_1$", size=16)
	plt.xlabel(r"$\beta_2$", size=16)
	cbar.ax.tick_params(labelsize=12)
	ax.tick_params(axis='both', which='major', labelsize=12)
	ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
	ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
	plt.savefig("output/figures/plots/Sup8A.svg", transparent=True)
	plt.show()

if "B" in run_panels:
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
	
	rmax1_vals = np.linspace(0.05, 0.2, num=10)
	rmax2_vals = np.linspace(0.05, 0.2, num=10)
	AMR1_vals = np.zeros((len(rmax1_vals), len(rmax2_vals)))
	AMR2_vals = np.zeros((len(rmax1_vals), len(rmax2_vals)))
	for i, rmax1_val in enumerate(rmax1_vals):
		bfunc1.rmax = rmax1_val
		for j, rmax2_val in enumerate(rmax2_vals):
			bfunc2.rmax = rmax2_val
			
			rep1.val = bfunc1.rmax
			res = mf.optimize_growth(growth, [rep2], [[0, bfunc2.rmax]])
			
			rep2.val = res[1][0]
			AMR1_vals[i][j] = gomp1.get_AMR()
			AMR2_vals[i][j] = gomp2.get_AMR()
	
	with open("output/figures/data/Sup8B.p", "wb") as f:
		pickle.dump([rmax1_vals, rmax2_vals, AMR1_vals, AMR2_vals], f)

if "B" in plot_panels:
	with open("output/figures/data/Sup8B.p", "rb") as f:
		rmax1_vals, rmax2_vals, AMR1_vals, AMR2_vals = pickle.load(f)
	AMR1_vals = AMR1_vals.transpose()
	AMR2_vals = AMR2_vals.transpose()
	
	# Load data from Fig5B for comparison
	with open("output/figures/data/Fig5B.p", "rb") as f:
		fig5_rmax1_vals, fig5_rmax2_vals, fig5_AMR1_vals, fig5_AMR2_vals = pickle.load(f)
	fig5_AMR1_vals = fig5_AMR1_vals.transpose()
	fig5_AMR2_vals = fig5_AMR2_vals.transpose()
	
	
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
				
				fig5_i = list(rmax1_vals).index(list(X)[i])
				fig5_j = list(rmax2_vals).index(list(Y)[i])
				if path_index == 1:
					is_ageless = (fig5_AMR1_vals[fig5_i, fig5_j] < 1e-5)
				else:
					is_ageless = (fig5_AMR2_vals[fig5_i, fig5_j] < 1e-5)
				
				star_color = "cyan" if not is_ageless else "white"
				ax.scatter(X[i], Y[i], s=10, marker=star_marker, c=star_color)
	
	
	fig = plt.figure(dpi=300)
	ax = fig.add_subplot(111)
	vmin = np.min([AMR1_vals, AMR2_vals])
	vmax = np.max([AMR1_vals, AMR2_vals]) + 0.001
	plot_mat(path_index=2, vmin=vmin, vmax=vmax, matrix=AMR1_vals)
	plot_mat(path_index=1, vmin=vmin, vmax=vmax, matrix=AMR2_vals)
	plt.xlim([0.04, 0.21])
	plt.ylim([0.04, 0.21])
	sm = plt.cm.ScalarMappable(cmap='hot', norm=plt.Normalize(vmin=vmin, vmax=vmax))
	sm._A = []
	cbar = plt.colorbar(sm)
	cbar.set_label(r"$\beta × b(r_{opt})$", size=16)
	plt.ylabel(r"$r_{max,1}$", size=16)
	plt.xlabel(r"$r_{max,2}$", size=16)
	cbar.ax.tick_params(labelsize=12)
	ax.tick_params(axis='both', which='major', labelsize=12)
	ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
	ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
	
	plt.savefig("output/figures/plots/Sup8B.svg", transparent=True)
	plt.show()


def bootstrap_ci(lst, n_bootstraps=10000, alpha=0.05):
	# Convert list to numpy array
	lst = np.array(lst)
	
	# Calculate observed proportion of True values
	obs_prop = np.mean(lst)
	
	# Create array to store bootstrap estimates
	bootstrap_props = np.zeros(n_bootstraps)
	
	# Perform bootstrapping
	for i in range(n_bootstraps):
		# Sample with replacement from the list
		resampled = np.random.choice(lst, size=len(lst), replace=True)
		# Calculate proportion of True values in resampled data
		bootstrap_props[i] = np.mean(resampled)
	
	# Calculate lower and upper bounds of confidence interval
	lower_bound = np.percentile(bootstrap_props, 100 * alpha / 2)
	upper_bound = np.percentile(bootstrap_props, 100 * (1 - alpha / 2))
	
	# Return observed proportion and confidence interval bounds
	return obs_prop, (lower_bound, upper_bound)


if "C" in run_panels:
	# To actually run the analysis, run MonteCarloPlei.py and MonteCarlo.py
	filepath = "output/figures/data/"
	
	plei_ageless_means = []
	plei_ageless_cis = []
	noplei_ageless_means = []
	noplei_ageless_cis = []

	for damclasses in [1, 2, 3, 4]:
		# Pleiotropy
		filename = f"MonteCarloPlei_{damclasses}.csv"
		modified_time = os.path.getmtime(filepath + filename)
		modified_datetime = datetime.datetime.fromtimestamp(modified_time)
		print(f"{filename} last modified on {modified_datetime.strftime('%Y-%m-%d %I:%M:%S %p')}")
		df = pd.read_csv(filepath + filename)
		
		ageless_classes = []
		for ageless_cat in range(damclasses):
			ageless_counts = list(df[f"ageless{ageless_cat + 1}"])
			ageless_classes.append(ageless_counts)
		
		one_ageless_list = []
		all_ageless_list = []
		for i in range(len(ageless_classes[0])):
			one_ageless_list.append(True if any([lst[i] for lst in ageless_classes]) else False)
			all_ageless_list.append(True if all([lst[i] for lst in ageless_classes]) else False)
		
		one_ageless_mean, one_ageless_ci = bootstrap_ci(one_ageless_list)
		all_ageless_mean, all_ageless_ci = bootstrap_ci(all_ageless_list)
		
		
		plei_ageless_means.append(all_ageless_mean)
		plei_ageless_cis.append(all_ageless_ci)
		
		# No pleiotropy
		filename = f"MonteCarlo_{damclasses}.csv"
		modified_time = os.path.getmtime(filepath + filename)
		modified_datetime = datetime.datetime.fromtimestamp(modified_time)
		print(f"{filename} last modified on {modified_datetime.strftime('%Y-%m-%d %I:%M:%S %p')}")
		df = pd.read_csv(filepath + filename)
		
		ageless_classes = []
		for ageless_cat in range(damclasses):
			ageless_counts = list(df[f"ageless{ageless_cat + 1}"])
			ageless_classes.append(ageless_counts)
		
		one_ageless_list = []
		all_ageless_list = []
		for i in range(len(ageless_classes[0])):
			one_ageless_list.append(True if any([lst[i] for lst in ageless_classes]) else False)
			all_ageless_list.append(True if all([lst[i] for lst in ageless_classes]) else False)
		
		one_ageless_mean, one_ageless_ci = bootstrap_ci(one_ageless_list)
		all_ageless_mean, all_ageless_ci = bootstrap_ci(all_ageless_list)
		
		
		noplei_ageless_means.append(all_ageless_mean)
		noplei_ageless_cis.append(all_ageless_ci)
	
	with open("output/figures/data/Sup8C.p", "wb") as f:
		pickle.dump([plei_ageless_means, plei_ageless_cis, noplei_ageless_means, noplei_ageless_cis], f)

if "C" in plot_panels:
	with open("output/figures/data/Sup8C.p", "rb") as f:
		plei_ageless_means, plei_ageless_cis, noplei_ageless_means, noplei_ageless_cis = pickle.load(f)
	dam_classes = [1, 2, 3, 4]
	
	fig = plt.figure(figsize=(6, 3), dpi=300)
	plt.plot(dam_classes, plei_ageless_means, marker="o", label=r"Pleiotropy")
	plt.fill_between(dam_classes, [plei_ageless_cis[x][0] for x in range(len(plei_ageless_means))],
	                 [plei_ageless_cis[x][1] for x in range(len(plei_ageless_means))], alpha=0.2)
	plt.plot(dam_classes, noplei_ageless_means, marker="o", label=r"No pleiotropy")
	plt.fill_between(dam_classes, [noplei_ageless_cis[x][0] for x in range(len(noplei_ageless_means))],
	                 [noplei_ageless_cis[x][1] for x in range(len(noplei_ageless_means))], alpha=0.2)
	plt.xticks(dam_classes)
	plt.ylabel("Fraction of Parameter Space\n"+r"Favoring All $r_{opt}=r_{max}$")
	plt.xlabel("Damage Classes")
	plt.legend(fontsize=9)
	plt.savefig("output/figures/plots/Sup8C.svg", transparent=True)
	plt.show()
