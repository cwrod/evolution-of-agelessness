import fcntl
import math
import sys
import matplotlib
import matplotlib.colors
from matplotlib import pyplot as plt
import numpy as np
import pickle
import model_funcs as mf



dam_classes = 2
should_overwrite = False
trials = 1000

for i,arg in enumerate(sys.argv[1:]):
	if arg == "-c":
		dam_classes = int(sys.argv[i+2])
	elif arg == "-n":
		trials = int(sys.argv[i+2])
	elif arg == "-f":
		should_overwrite = True


if should_overwrite:
	with open("output/figures/data/MonteCarloSplit_"+str(dam_classes)+".csv","w") as f:
		header = []
		for class_num in range(dam_classes):
			header.append(f"rmax{class_num + 1}")
			header.append(f"beta{class_num + 1}")
		header.append("opt_growth")
		for class_num in range(dam_classes):
			header.append(f"ropt{class_num + 1}")
			header.append(f"ageless{class_num + 1}")
		f.write(",".join(header)+"\n")
	sys.exit()

rep1 = mf.Repair(0.32)
rep_vars = [mf.Repair(0.32) for x in range(dam_classes)]
bfunc_vars = [mf.InverseTOff(rep_vars[x], 0.5, 0.2) for x in range(dam_classes)]
hfunc = mf.ReproductiveScalingMulti(rep_vars)
gomp_vars = [mf.GompertzRepAMR(0.0006/dam_classes, bfunc_vars[x]) for x in range(dam_classes)]
extr = mf.Extrinsic(0.03)
prv = mf.PeakReproductiveValue(hfunc, 1)
mat_age = mf.MatAge(8, prv)
rep_func = mf.DecayingRepro(gomp_vars, mat_age, prv)
surv_func = mf.SurvFunc(gomp_vars + [extr])
growth = mf.KSelPop(rep_func, surv_func)

rmax_total_vals = np.random.uniform(1e-5,0.2,size=trials)
beta_vals = np.random.uniform(1e-5,0.3,size=trials)
perc_vals = []
if dam_classes > 1:
	rmax_split_points = np.sort(np.random.uniform(size=(trials,dam_classes-1)),axis=1)
	for row in rmax_split_points:
		row = list(row)
		row.insert(0,0)
		row.append(1)
		perc_vals.append(list(np.diff(row)))
else:
	perc_vals = [1]*trials
	
rmax_vals = np.asarray(perc_vals)*rmax_total_vals[:,np.newaxis]

for trial in range(trials):
	entry = []
	print(trial)
	print(f"rmax total: {rmax_total_vals[trial]}")
	for class_num in range(dam_classes):
		print(
			f"rmax{class_num + 1} = {rmax_vals[trial][class_num]}, beta{class_num + 1} = {beta_vals[trial]}, ",
			end="")
		bfunc_vars[class_num].rmax = rmax_vals[trial][class_num]
		bfunc_vars[class_num].beta = beta_vals[trial]
		entry.append(rmax_vals[trial][class_num])
		entry.append(beta_vals[trial])
	print()
	
	rep_ranges = [[0, bfunc_vars[x].rmax] for x in range(dam_classes)]
	res = mf.optimize_growth(growth, rep_vars, rep_ranges)
	print(res)
	
	entry.append(res[0])
	for class_num in range(dam_classes):
		entry.append(res[1][class_num])
		entry.append(res[1][class_num] == bfunc_vars[class_num].rmax)
	print(entry)
	with open("output/figures/data/MonteCarloSplit_" + str(dam_classes) + ".csv", "a") as f:
		fcntl.flock(f.fileno(), fcntl.LOCK_EX)
		f.write(",".join(list(map(str, entry))) + "\n")
		fcntl.flock(f.fileno(), fcntl.LOCK_UN)
