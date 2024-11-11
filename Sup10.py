import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

#run_panels = ["A","B","C","D","E","F","G","H"]
run_panels = []
#plot_panels = ["A","B","C","D","E","F","G","H"]
plot_panels = ["D","H"]


for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]



if "A" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.StepTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp,extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	b_vals = []
	a_vals = []
	r_vals = list(np.linspace(0.001,bfunc.rmax,num=1000))
	for r_val in r_vals:
		rep1.val = r_val
		b_vals.append(bfunc.ret_val())
		a_vals.append(mat_age.ret_val())
	
	
	with open("output/figures/data/Sup10A.p", "wb") as f:
		pickle.dump([r_vals, b_vals, a_vals], f)

if "A" in plot_panels:
	with open("output/figures/data/Sup10A.p", "rb") as f:
		r_vals, b_vals, a_vals = pickle.load(f)

	fig,ax = plt.subplots(figsize=(3,2),dpi=300)
	ax.plot(r_vals,b_vals,color=fss.COLORS[0],linewidth=2)
	ax.set_xlabel("r")
	ax.set_ylim(-0.5,10)
	ax.set_ylabel("b(r)",color=fss.COLORS[0])
	
	ax2=ax.twinx()
	ax2.plot(r_vals, a_vals,color=fss.COLORS[1],linewidth=2)
	ax2.set_ylabel("a(r)",color=fss.COLORS[1])
	#ax2.set_ylim(0.05,1)
	plt.savefig("output/figures/plots/Sup10A.svg",transparent=True)
	plt.show()


if "B" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.StepTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.01,bfunc.rmax,num=100))
	g_vals = []
	
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Sup10B.p", "wb") as f:
		pickle.dump([r_vals,g_vals], f)
	
if "B" in plot_panels:
	with open("output/figures/data/Sup10B.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)

	fig = plt.figure(figsize=(3, 2), dpi=300)
	ax = fig.add_subplot()
	ax.plot(r_vals, g_vals, linewidth=2)
	ax.set_ylim(None,0.16)
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]

	ax.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
				xy=(arg_max, max_val), xycoords="data",
				xytext=(-30, 7), textcoords="offset points", zorder=2.5)
	plt.xlabel("r")
	plt.ylabel("g")
	plt.savefig("output/figures/plots/Sup10B.svg",transparent=True)
	plt.show()

if "C" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.StepTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.01,bfunc.rmax,num=100))
	g_vals = []
	
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Sup10C.p", "wb") as f:
		pickle.dump([r_vals, g_vals], f)

if "C" in plot_panels:
	with open("output/figures/data/Sup10C.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)
	
	fig = plt.figure(figsize=(3, 2), dpi=300)
	ax = fig.add_subplot()
	ax.plot(r_vals, g_vals, linewidth=2)
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]
	
	ax.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
	            xy=(arg_max, max_val), xycoords="data",
	            xytext=(-60, 7), textcoords="offset points", zorder=2.5)
	plt.ylim(None,16)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.savefig("output/figures/plots/Sup10C.svg", transparent=True)
	plt.show()

if "D" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.StepTOff(rep1, 0.4, 0.2)
	bfunc2 = mf.StepTOff(rep2, 0.4, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
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
	with open("output/figures/data/Sup10D.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "D" in plot_panels:
	with open("output/figures/data/Sup10D.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, dpi=300, figsize=(3,3))
	surf = ax.plot_surface(r1_vals, r2_vals, g_vals, cmap="coolwarm")
	ax.zaxis.set_rotate_label(False)
	ax.set_ylabel(r"$r_1$")
	ax.set_xlabel(r"$r_2$")
	ax.set_zlabel(r"$R_0$")
	ax.view_init(60, -105)
	plt.savefig("output/figures/plots/Sup10D.svg", transparent=True, bbox_inches='tight')
	plt.show()

if "E" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.SmoothStepTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	b_vals = []
	a_vals = []
	r_vals = list(np.linspace(0.001, bfunc.rmax, num=1000))
	for r_val in r_vals:
		rep1.val = r_val
		b_vals.append(bfunc.ret_val())
		a_vals.append(mat_age.ret_val())
	
	with open("output/figures/data/Sup10E.p", "wb") as f:
		pickle.dump([r_vals, b_vals, a_vals], f)

if "E" in plot_panels:
	with open("output/figures/data/Sup10E.p", "rb") as f:
		r_vals, b_vals, a_vals = pickle.load(f)
	
	fig, ax = plt.subplots(figsize=(3, 2), dpi=300)
	ax.plot(r_vals, b_vals, color=fss.COLORS[0], linewidth=2)
	ax.set_xlabel("r")
	ax.set_ylim(-0.025, 0.5)
	ax.set_ylabel("b(r)", color=fss.COLORS[0])
	
	ax2 = ax.twinx()
	ax2.plot(r_vals, a_vals, color=fss.COLORS[1], linewidth=2)
	ax2.set_ylabel("a(r)", color=fss.COLORS[1])
	# ax2.set_ylim(0.05,1)
	plt.savefig("output/figures/plots/Sup10E.svg", transparent=True)
	plt.show()

if "F" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.SmoothStepTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 0.25)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.RSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.01, bfunc.rmax, num=100))
	g_vals = []
	
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Sup10F.p", "wb") as f:
		pickle.dump([r_vals, g_vals], f)

if "F" in plot_panels:
	with open("output/figures/data/Sup10F.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)
	
	fig = plt.figure(figsize=(3, 2), dpi=300)
	ax = fig.add_subplot()
	ax.plot(r_vals, g_vals, linewidth=2)
	ax.set_ylim(None, 0.06)
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]
	
	ax.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
	            xy=(arg_max, max_val), xycoords="data",
	            xytext=(-30, 7), textcoords="offset points", zorder=2.5)
	plt.xlabel("r")
	plt.ylabel("g")
	plt.savefig("output/figures/plots/Sup10F.svg", transparent=True)
	plt.show()

if "G" in run_panels:
	rep1 = mf.Repair(0.32)
	bfunc = mf.SmoothStepTOff(rep1, 0.5, 0.2)
	hfunc = mf.ReproductiveScaling(rep1)
	gomp = mf.GompertzRepAMR(0.0006, bfunc)
	extr = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp], mat_age, prv)
	surv_func = mf.SurvFunc([gomp, extr])
	growth = mf.KSelPop(rep_func, surv_func)
	
	r_vals = list(np.linspace(0.01, bfunc.rmax, num=100))
	g_vals = []
	
	for r_val in r_vals:
		rep1.val = r_val
		g_vals.append(growth.ret_val())
	
	with open("output/figures/data/Sup10G.p", "wb") as f:
		pickle.dump([r_vals, g_vals], f)

if "G" in plot_panels:
	with open("output/figures/data/Sup10G.p", "rb") as f:
		r_vals, g_vals = pickle.load(f)
	
	fig = plt.figure(figsize=(3, 2), dpi=300)
	ax = fig.add_subplot()
	ax.plot(r_vals, g_vals, linewidth=2)
	max_val = max(g_vals)
	arg_max = r_vals[g_vals.index(max_val)]
	
	ax.scatter(arg_max, max_val, color="black", zorder=2.5)
	arg_max = round(arg_max, 2)
	max_val = round(max_val, 3)
	ax.annotate("(" + str(arg_max) + ", " + str(max_val) + ")",
	            xy=(arg_max, max_val), xycoords="data",
	            xytext=(-60, 7), textcoords="offset points", zorder=2.5)
	plt.ylim(None, 16)
	plt.xlabel("r")
	plt.ylabel(r"$R_0$")
	plt.savefig("output/figures/plots/Sup10G.svg", transparent=True)
	plt.show()

if "H" in run_panels:
	rep1 = mf.Repair(0.32)
	rep2 = mf.Repair(0.32)
	bfunc1 = mf.SmoothStepTOff(rep1, 0.4, 0.2)
	bfunc2 = mf.SmoothStepTOff(rep2, 0.4, 0.2)
	hfunc = mf.ReproductiveScalingMulti([rep1, rep2])
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
	with open("output/figures/data/Sup10H.p", "wb") as f:
		pickle.dump([r1_vals, r2_vals, g_vals], f)

if "H" in plot_panels:
	with open("output/figures/data/Sup10H.p", "rb") as f:
		r1_vals, r2_vals, g_vals = pickle.load(f)
	
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, dpi=300, figsize=(3,3))
	surf = ax.plot_surface(r1_vals, r2_vals, g_vals, cmap="coolwarm")
	ax.zaxis.set_rotate_label(False)
	ax.set_ylabel(r"$r_1$")
	ax.set_xlabel(r"$r_2$")
	ax.set_zlabel(r"$R_0$")
	ax.view_init(60, -105)
	plt.savefig("output/figures/plots/Sup10H.svg", transparent=True, bbox_inches='tight')
	plt.show()
