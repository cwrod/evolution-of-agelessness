# Evolution of aging/agelessness under the Disposable Soma Theory

This is the code that produces the figures for the paper "Agelessness is Possible Under the Disposable Soma Theory but System Complexity Makes it Unlikely."

## Installation

Run the following code to install the repository and all the dependencies.

```
git clone https://github.com/cwrod/evolution-of-agelessness
python3 -m venv
venv/bin/pip install -r requirements.txt
```

## Reproducing Figures

You can reproduce figures by calling FigX.py (or SupX.py) and providing the panels you want to reproduce. For instance, to reproduce panels A, B, and C from Figure 1, you could type the following.

```
venv/bin/python3 Fig1.py ABC
```

If you just call FigX.py, it will default to reproducing all panels. The program saves the results from the analysis under `output/figures/data/` and can retrieve pickled files from this folder to save analysis time. If you set `run_panels = []` in the first few lines of the python files, it will default to using cached results.

The sampled data and results for the Monte Carlo simulation is saved in the repository. To rerun this result, you can run the MonteCarlo.py, MonteCarloPlei.py, and MonteCarloSplit.py files. You can use the -c flag to adjust the number of damage classes and the -n flag to adjust the number of samples taken. Also, the program currently defaults to appending results to the CSV file. To overwrite the previous results, you can use the -f flag. So to run a Monte Carlo simulation 1000 times with 3 damage classes and overwriting the previous results, you would run the following.

```
venv/bin/pyton3 MonteCarlo.py -c 3 -n 1000 -f
```

You can modify equation forms and default parameter values in model_funcs.py. 
