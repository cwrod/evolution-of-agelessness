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


if "A" in run_panels:
	# To actually run the analysis, run MonteCarlo.py
	filepath = "output/figures/data/"
	
	one_ageless_means = []
	one_ageless_cis = []
	all_ageless_means = []
	all_ageless_cis = []
	
	for damclasses in [1, 2, 3, 4]:
		filename = f"MonteCarlo_{damclasses}.csv"
		modified_time = os.path.getmtime(filepath + filename)
		modified_datetime = datetime.datetime.fromtimestamp(modified_time)
		print(f"{filename} last modified on {modified_datetime.strftime('%Y-%m-%d %I:%M:%S %p')}")
		df = pd.read_csv(filepath + filename)
		
		ageless_classes = []
		for ageless_cat in range(damclasses):
			ageless_counts = list(df[f"ageless{ageless_cat + 1}"])
			ageless_classes.append(ageless_counts)
			print(ageless_counts)
		
		one_ageless_list = []
		all_ageless_list = []
		for i in range(len(ageless_classes[0])):
			one_ageless_list.append(True if any([lst[i] for lst in ageless_classes]) else False)
			all_ageless_list.append(True if all([lst[i] for lst in ageless_classes]) else False)
		print(one_ageless_list)
		print(all_ageless_list)
		print()
		
		one_ageless_mean, one_ageless_ci = bootstrap_ci(one_ageless_list)
		all_ageless_mean, all_ageless_ci = bootstrap_ci(all_ageless_list)
		
		print((one_ageless_mean, one_ageless_ci))
		print((all_ageless_mean, all_ageless_ci))
		
		one_ageless_means.append(one_ageless_mean)
		one_ageless_cis.append(one_ageless_ci)
		
		all_ageless_means.append(all_ageless_mean)
		all_ageless_cis.append(all_ageless_ci)
	
	with open("output/figures/data/Fig6A.p", "wb") as f:
		pickle.dump([one_ageless_means, one_ageless_cis, all_ageless_means, all_ageless_cis], f)

if "A" in plot_panels:
	with open("output/figures/data/Fig6A.p", "rb") as f:
		one_ageless_means, one_ageless_cis, all_ageless_means, all_ageless_cis = pickle.load(f)
	dam_classes = [1, 2, 3, 4]
	
	fig = plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(dam_classes, one_ageless_means, marker="o", label=r"Some $r_{opt}=r_{max}$")
	plt.fill_between(dam_classes, [one_ageless_cis[x][0] for x in range(len(one_ageless_means))],
	                 [one_ageless_cis[x][1] for x in range(len(one_ageless_means))], alpha=0.2)
	plt.plot(dam_classes, all_ageless_means, marker="o", label=r"All $r_{opt}=r_{max}$")
	plt.fill_between(dam_classes, [all_ageless_cis[x][0] for x in range(len(all_ageless_means))],
	                 [all_ageless_cis[x][1] for x in range(len(all_ageless_means))], alpha=0.2)
	plt.xticks(dam_classes)
	plt.ylabel("Fraction of Parameter Space")
	plt.xlabel("Damage Classes")
	plt.legend(fontsize=9)
	plt.title("Adding")
	plt.ylim(-0.05,1.05)
	plt.savefig("output/figures/plots/Fig6A.svg", transparent=True)
	plt.show()
	
if "B" in run_panels:
	# To actually run the analysis, run MonteCarloSplit.py
	filepath = "output/figures/data/"
	
	one_ageless_means = []
	one_ageless_cis = []
	all_ageless_means = []
	all_ageless_cis = []
	
	for damclasses in [1, 2, 3, 4]:
		filename = f"MonteCarloSplit_{damclasses}.csv"
		modified_time = os.path.getmtime(filepath + filename)
		modified_datetime = datetime.datetime.fromtimestamp(modified_time)
		print(f"{filename} last modified on {modified_datetime.strftime('%Y-%m-%d %I:%M:%S %p')}")
		df = pd.read_csv(filepath + filename)
		
		ageless_classes = []
		for ageless_cat in range(damclasses):
			ageless_counts = list(df[f"ageless{ageless_cat + 1}"])
			ageless_classes.append(ageless_counts)
			print(ageless_counts)
		
		one_ageless_list = []
		all_ageless_list = []
		for i in range(len(ageless_classes[0])):
			one_ageless_list.append(True if any([lst[i] for lst in ageless_classes]) else False)
			all_ageless_list.append(True if all([lst[i] for lst in ageless_classes]) else False)
		
		one_ageless_mean, one_ageless_ci = bootstrap_ci(one_ageless_list)
		all_ageless_mean, all_ageless_ci = bootstrap_ci(all_ageless_list)
		
		print((one_ageless_mean, one_ageless_ci))
		print((all_ageless_mean, all_ageless_ci))
		
		one_ageless_means.append(one_ageless_mean)
		one_ageless_cis.append(one_ageless_ci)
		
		all_ageless_means.append(all_ageless_mean)
		all_ageless_cis.append(all_ageless_ci)
	
	with open("output/figures/data/Fig6B.p", "wb") as f:
		pickle.dump([one_ageless_means, one_ageless_cis, all_ageless_means, all_ageless_cis], f)

if "B" in plot_panels:
	with open("output/figures/data/Fig6B.p", "rb") as f:
		one_ageless_means, one_ageless_cis, all_ageless_means, all_ageless_cis = pickle.load(f)
	dam_classes = [1, 2, 3, 4]
	
	fig = plt.figure(figsize=(4, 3), dpi=300)
	plt.plot(dam_classes, one_ageless_means, marker="o", label=r"Some $r_{opt}=r_{max}$")
	plt.fill_between(dam_classes, [one_ageless_cis[x][0] for x in range(len(one_ageless_means))],
	                 [one_ageless_cis[x][1] for x in range(len(one_ageless_means))], alpha=0.2)
	plt.plot(dam_classes, all_ageless_means, marker="o", label=r"All $r_{opt}=r_{max}$")
	plt.fill_between(dam_classes, [all_ageless_cis[x][0] for x in range(len(all_ageless_means))],
	                 [all_ageless_cis[x][1] for x in range(len(all_ageless_means))], alpha=0.2)
	plt.xticks(dam_classes)
	plt.ylabel("Fraction of Parameter Space")
	plt.xlabel("Damage Classes")
	plt.legend(fontsize=9)
	plt.title("Splitting")
	plt.ylim(-0.05,1.05)
	plt.savefig("output/figures/plots/Fig6B.svg", transparent=True)
	plt.show()
