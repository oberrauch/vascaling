# Code



- `commitment_run_oggm.py`: Run the OGGM flowline model for all glaciers of the Rofental Basin under a random climate based on the HistAlp climate between 1984 and 2014. The random model is initialized with the same seed (12), the get reproducible results. The run output is stored with the suffix `_commitment` in the working directory (`../working_directories/commitment_run_oggm/`) as netCDF file.

- `commitment_run_vas.py`: Run the volume/area scaling model for all glaciers of the Rofental Basin under a random climate based on the HistAlp climate between 1984
  and 2014. The random model is initialized with the same seed (12), the get reproducible results. The run output is stored with the suffix `_commitment` in the working directory (`../working_directories/commitment_run_vas/`) as netCDF file.

- compare.py

- comparison_function.py

- eq_tmp.py

- equilibrium_run.py

- eval_eq_runs.py

- matlab_comparison.py

- mb_calibration.py

- mb_calibration_vas.py

- mb_calibration_vas_RGI6_HISTALP.py

- model_ben_matlab

- prepo_template.py

- ret_t_star_list.py

- `run_ben.py`: Set up a full run with the volume/area scaling model from start to finish. This includes: a) initialisation and calibration, b) the mass balance model, c) the 'dynamic' model. 

  The results are stored under `run_ben.csv` in the `../data/` directory. The plots show the temporal evolution of glacier length, surface area and volume and are store under `length/area/volume_ben.png` in the `../plots/` directory. 

- `run_comparison.py`: @WEITERMACHEN

- `run_oggm.py`: Set up a full run with the volume/area scaling model from start to finish. This includes: a) initialisation and calibration, b) the mass balance model, c) the 'dynamic' model. 

  The results are stored under `run_oggm.csv` in the `../data/` directory. The plots show the temporal evolution of glacier length, surface area and volume and are store under `length/area/volume_oggm.png` in the `../plots/` directory. 

- seek_start_area.py

- set_up_run.py

- start_area.py

- tmp.py

- unused_snippets.py