## Setting up an analysis using DeepDiveR – a worked example

You can find instructions for installing the DeepDive programme [here](https://github.com/DeepDive-project/deepdive) and the DeepDiveR library [here.](https://github.com/DeepDive-project/DeepDiveR)

Once the programme has been installed, the R script `carnivora_runner.R` takes input occurrence data from `Carnivora_data.csv`.
The script then defines time bins to be used in the analysis and the present diversity of the clade, i.e. the number of living species in the clade. This will be used to inform the model.
Additional settings can be specified such as the number of replicates (age randomizations), the number of simulations, and information about area connectivity.
The output of this script is formatted input data and a config file.

Once the configuration and input files are created, the full DeepDive analysis, inclusive of simulation, model training and empirical predictions, can be carried out through a single command line entered in a Terminal (MacOS and Linux) or Command prompt (Windows) window executing the Python script [`run_dd_config.py`](https://github.com/DeepDive-project/deepdive/run_dd_config.py) as follows:

```
python run_dd_config.py your_path/config_file.ini
```

You can additionally specify a working directory where all output files are saved and the number of CPUs used for the parallelized simulations, which will overwrite the corresponding settings in the configuration file, using the flags `-wd` and `-cpu`, e.g. 
```
python run_dd_config.py your_path/config_file.ini -wd your_working_directory -cpu 64
``` 

This script will create a "simulations" folder containing the training and test sets, and a "trained_models" folder containing the trained models and plots of the training history. This folder will additionally include plots comparing the empirical and simulated fossil features (e.g. number of occurrences through time and per area, number of localities, fraction of singletons, and sampled diversity), CSV files with the predicted diversity trajectories for the test set and for the empirical dataset, and a plot of the estimated diversity trajectory.
