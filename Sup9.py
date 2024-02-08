import math
import sys
from matplotlib import pyplot as plt
import numpy as np
import pickle
import operator
import model_funcs as mf
import FigStyleSchemes as fss

run_panels = ["A","B","C","D","E","F","G","H"]
plot_panels = ["A","B","C","D","E","F","G","H"]

for arg in sys.argv[1:]:
	run_panels = [x for x in arg.strip()]
	plot_panels = [x for x in arg.strip()]


# Plots: beta, rmax, gamma, drain, mort
def gen_panel_data(growth_FR_single, rep_var_single, dam_sing_ret_func,
		beta_sing_assign_func, rmax_sing_assign_func, gamma_sing_assign_func, drain_assign_func, mort_assign_func,
		beta_vals, rmax_vals, gamma_vals, drain_vals, mort_vals,
        beta_sing_default, rmax_sing_default, gamma_sing_default,
		output_file):
	plot_data = {}
	
	beta_sing_assign_func(beta_sing_default)
	rmax_sing_assign_func(rmax_sing_default)
	gamma_sing_assign_func(gamma_sing_default)
	drain_assign_func(0)
	mort_assign_func(0)
	dam_vals = []
	for beta_val in beta_vals:
		beta_sing_assign_func(beta_val)
		max_growth, max_rep_point = mf.optimize_growth(growth_FR_single, [rep_var_single], [[0, rmax_sing_default]])
		rep_var_single.val = max_rep_point[0]
		dam_vals.append(dam_sing_ret_func())
	plot_data["beta_sing_dams"] = dam_vals
	
	beta_sing_assign_func(beta_sing_default)
	rmax_sing_assign_func(rmax_sing_default)
	gamma_sing_assign_func(gamma_sing_default)
	drain_assign_func(0)
	mort_assign_func(0)
	dam_vals = []
	for rmax_val in rmax_vals:
		rmax_sing_assign_func(rmax_val)
		max_growth, max_rep_point = mf.optimize_growth(growth_FR_single, [rep_var_single], [[0, rmax_val]])
		rep_var_single.val = max_rep_point[0]
		dam_vals.append(dam_sing_ret_func())
	plot_data["rmax_sing_dams"] = dam_vals

	beta_sing_assign_func(beta_sing_default)
	rmax_sing_assign_func(rmax_sing_default)
	gamma_sing_assign_func(gamma_sing_default)
	drain_assign_func(0)
	mort_assign_func(0)
	dam_vals = []
	for gamma_val in gamma_vals:
		gamma_sing_assign_func(gamma_val)
		max_growth, max_rep_point = mf.optimize_growth(growth_FR_single, [rep_var_single], [[0, rmax_sing_default]])
		rep_var_single.val = max_rep_point[0]
		dam_vals.append(dam_sing_ret_func())
	plot_data["gamma_sing_dams"] = dam_vals
	
	beta_sing_assign_func(beta_sing_default)
	rmax_sing_assign_func(rmax_sing_default)
	gamma_sing_assign_func(gamma_sing_default)
	drain_assign_func(0)
	mort_assign_func(0)
	dam_vals = []
	for drain_val in drain_vals:
		drain_assign_func(drain_val)
		max_growth, max_rep_point = mf.optimize_growth(growth_FR_single, [rep_var_single], [[0, rmax_sing_default]])
		rep_var_single.val = max_rep_point[0]
		dam_vals.append(dam_sing_ret_func())
	plot_data["drain_dams"] = dam_vals
	
	beta_sing_assign_func(beta_sing_default)
	rmax_sing_assign_func(rmax_sing_default)
	gamma_sing_assign_func(gamma_sing_default)
	drain_assign_func(0)
	mort_assign_func(0)
	dam_vals = []
	for mort_val in mort_vals:
		mort_assign_func(mort_val)
		max_growth, max_rep_point = mf.optimize_growth(growth_FR_single, [rep_var_single], [[0, rmax_sing_default]])
		rep_var_single.val = max_rep_point[0]
		dam_vals.append(dam_sing_ret_func())
	plot_data["mort_dams"] = dam_vals
	
	
	plot_val_ranges = {"beta_vals":beta_vals,"rmax_vals":rmax_vals,"gamma_vals":gamma_vals,
	                   "drain_vals" : drain_vals, "mort_vals" : mort_vals,}
	with open(output_file, "wb") as f:
		pickle.dump([plot_data, plot_val_ranges], f)


def gen_panel_plot(data_file, plot_file, fig_title):
	with open(data_file, "rb") as f:
		plot_data, plot_val_ranges = pickle.load(f)

	fig, axes = plt.subplots(1, 5, figsize=(15, 3), dpi=300)
	beta_vals = plot_val_ranges["beta_vals"]
	dam_vals = plot_data["beta_sing_dams"]
	ax = axes[0]
	ax.plot(beta_vals, dam_vals)
	ax.set_ylabel(r"$\beta * b\left(r_{opt}\right)$")
	ax.set_xlabel(r"$\beta$")
	ax.set_xscale("log")
	
	rmax_vals = plot_val_ranges["rmax_vals"]
	dam_vals = plot_data["rmax_sing_dams"]
	ax = axes[1]
	ax.plot(rmax_vals, dam_vals)
	ax.set_xlabel(r"$r_{max}$")
	ax.set_ylabel(r"$\beta * b\left(r_{opt}\right)$")
	
	gamma_vals = plot_val_ranges["gamma_vals"]
	dam_vals = plot_data["gamma_sing_dams"]
	ax = axes[2]
	ax.plot(gamma_vals, dam_vals)
	ax.set_xlabel(r"$\gamma$")
	ax.set_ylabel(r"$\beta * b\left(r_{opt}\right)$")
	
	drain_vals = plot_val_ranges["drain_vals"]
	dam_vals = plot_data["drain_dams"]
	ax = axes[3]
	ax.plot(drain_vals, dam_vals)
	ax.set_xlabel(r"$d$")
	ax.set_ylabel(r"$\beta * b\left(r_{opt}\right)$")
	
	mort_vals = plot_val_ranges["mort_vals"]
	dam_vals = plot_data["mort_dams"]
	ax = axes[4]
	ax.plot(mort_vals, dam_vals)
	ax.set_xlabel(r"$AMR_2$")
	ax.set_ylabel(r"$\beta * b\left(r_{opt}\right)$")
	
	
	fig.suptitle(fig_title)
	plt.savefig(plot_file, transparent=True)
	plt.show()


if "A" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.InverseTOff(rep_var_single, 0.5, 0.2)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single,0.0)
	gomp_single = mf.GompertzRepAMR(0.0006, bfunc_single)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_AMR,
			lambda x : setattr(bfunc_single, "beta", x), lambda x : setattr(bfunc_single, "rmax", x),
			lambda x: setattr(extr_single, "gamma", x), lambda x : setattr(hfunc_single, "d", x), lambda x: setattr(gomp_single_cons, "B", x),
			list(np.logspace(-3,1,num=50)), list(np.linspace(0.001,0.9,num=50)),
			list(np.linspace(0,0.05,num=50)), list(np.linspace(0,0.4,num=50)), list(np.linspace(0,0.03,num=50)),
			0.2, 0.5, 0.03, "output/figures/data/Sup9A.p")
	
if "A" in plot_panels:
	fig_title = "Gompertz - Repair on AMR"
	gen_panel_plot("output/figures/data/Sup9A.p", "output/figures/plots/Sup9A.svg", fig_title)

if "B" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.InverseTOff(rep_var_single, 0.5, 1e-6)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.GompertzRepIMR(bfunc_single, 0.2)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_IMR,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-8, -1, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               1e-6, 0.5, 0.03,
	               "output/figures/data/Sup9B.p")

if "B" in plot_panels:
	fig_title = "Gompertz - Repair on IMR"
	gen_panel_plot("output/figures/data/Sup9B.p", "output/figures/plots/Sup9B.svg", fig_title)

if "C" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.InverseTOff(rep_var_single, 0.5, 5)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.WeibullRepPower(4e-10, bfunc_single)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_power,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-1, 10, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               5, 0.5, 0.03,
	               "output/figures/data/Sup9C.p")

if "C" in plot_panels:
	fig_title = "Weibull - Repair on Exponent"
	gen_panel_plot("output/figures/data/Sup9C.p", "output/figures/plots/Sup9C.svg", fig_title)

if "D" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.InverseTOff(rep_var_single, 0.5, 4e-10)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.WeibullRepCoeff(bfunc_single, 5)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_coeff,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-15, -5, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               4e-10, 0.5, 0.03,
	               "output/figures/data/Sup9D.p")

if "D" in plot_panels:
	fig_title = "Weibull - Repair on Coefficient"
	gen_panel_plot("output/figures/data/Sup9D.p", "output/figures/plots/Sup9D.svg", fig_title)

if "E" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.SqrtTOff(rep_var_single, 0.5, 0.2)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.GompertzRepAMR(0.0006, bfunc_single)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_FR_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_FR_single, rep_var_single, gomp_single.get_AMR,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-3, 1, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               0.2, 0.5, 0.01,
	               "output/figures/data/Sup9E.p")

if "E" in plot_panels:
	fig_title = "Gompertz - Repair on AMR - Sqrt"
	gen_panel_plot("output/figures/data/Sup9E.p", "output/figures/plots/Sup9E.svg", fig_title)

if "F" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.SquareTOff(rep_var_single, 0.5, 0.2)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.GompertzRepAMR(0.0006, bfunc_single)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.01)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_AMR,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-3, 1, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               0.2, 0.5, 0.01,
	               "output/figures/data/Sup9F.p")

if "F" in plot_panels:
	fig_title = "Gompertz - Repair on AMR - Square"
	gen_panel_plot("output/figures/data/Sup9F.p", "output/figures/plots/Sup9F.svg", fig_title)

if "G" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.ExpTOff(rep_var_single, 0.5, 0.2)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.GompertzRepAMR(0.0006, bfunc_single)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.MatAge(8, prv)
	rep_func = mf.DecayingRepro([gomp_single], mat_age, prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_AMR,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-3, 1, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               0.2, 0.5, 0.03,
	               "output/figures/data/Sup9G.p")

if "G" in plot_panels:
	fig_title = "Gompertz - Repair on AMR - Exp"
	gen_panel_plot("output/figures/data/Sup9G.p", "output/figures/plots/Sup9G.svg", fig_title)

if "H" in run_panels:
	rep_var_single = mf.Repair(0.32)
	bfunc_single = mf.InverseTOff(rep_var_single, 0.5, 0.2)
	hfunc_single = mf.ReproductiveScalingDrain(rep_var_single, 0.0)
	gomp_single = mf.GompertzRepAMR(0.0006, bfunc_single)
	gomp_single_cons = mf.GompertzCons(0.0006, 0.0)
	extr_single = mf.Extrinsic(0.03)
	prv = mf.PeakReproductiveValue(hfunc_single, 1)
	mat_age = mf.ConstMatAge(0)
	rep_func = mf.ConstRepro(mat_age,prv)
	surv_func = mf.SurvFunc([gomp_single, gomp_single_cons, extr_single])
	growth_K_single = mf.KSelPop(rep_func, surv_func)
	
	gen_panel_data(growth_K_single, rep_var_single, gomp_single.get_AMR,
	               lambda x: setattr(bfunc_single, "beta", x), lambda x: setattr(bfunc_single, "rmax", x),
	               lambda x: setattr(extr_single, "gamma", x), lambda x: setattr(hfunc_single, "d", x),
	               lambda x: setattr(gomp_single_cons, "B", x),
	               list(np.logspace(-3, 1, num=50)), list(np.linspace(0.001, 0.9, num=50)),
	               list(np.linspace(0, 0.05, num=50)), list(np.linspace(0, 0.4, num=50)),
	               list(np.linspace(0, 0.03, num=50)),
	               0.2, 0.5, 0.03,
	               "output/figures/data/Sup9H.p")

if "H" in plot_panels:
	fig_title = "Gompertz - Repair on AMR - Constant Reproduction"
	gen_panel_plot("output/figures/data/Sup9H.p", "output/figures/plots/Sup9H.svg", fig_title)

