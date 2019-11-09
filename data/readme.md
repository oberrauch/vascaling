# Data directory

- `eq_runs/`: The directory contains the results produced by the `../code/equilibrium_run.py` script, performing equilibrium runs using constant and random massbalance models.
- `run_comparison.csv`: Result of the `../code/run_comparison.py` script, which runs and compares the VAS model and the OGGM flowline model for the Upper Grindelwald Glacier.
- `start_area_results/`: The directory contains Matlab workspace dumps, from the iterative start area seeking process. Each file corresponds to one iteration step and contains yearly values of glacier length, surface area and volume, as well as the corresponding time scales (response times) for length and area change ($\tau_L$ and $\tau_A$).
- `mat_files/`: Directory contains other workspace dumps from the original Matlab model.
- `rgi_ref_glaciers.csv`: The table contains all the *reference* glaciers used by the OGGM, i.e., all glaciers with a substantially long massbalance records.
- `run_oggm.csv`: Result of the `../code/run_oggm.py` script, which sets up and runs OGGM flowline model from start to finish (one of the first scripts.)
- `run_vas.csv`: Result of the `../code/run_vas.py` script, which sets up and runs VAS model from start to finish (one of the first scripts.)