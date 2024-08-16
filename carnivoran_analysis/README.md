# carnivore_analysis

The R script carnivora_runner.R takes input occurrence data from carnivora_data_cleaned.xlsx, defines time bins, present diversity of the clade and the number of replicates to produce formatted input data and a configuration file.

Once the configuration and input files are created, the full DeepDive analysis, inclusive of simulation, model training and empirical predictions, can be carried out through a single command line entered in a Terminal (MacOS and Linux) or Command prompt (Windows) window using the Python script run_dd_config.py:

```
python run_dd_config.py your_path/config_file.ini
```

You can additionally specify a working directory where all output files are saved and the number of CPUs used for the parallelized simulations, which will overwrite the corresponding settings in the configuration file, using the flags 
```
-wd your\_working\_directory} 
```
and e.g. -cpu 64. 

This script will create a "simulations" folder containing the training and test sets, and a "trained_models" folder containing the trained models and plots of the training history. This folder will additionally include plots comparing the empirical and simulated fossil features (e.g. number of occurrences through time and per area, number of localities, fraction of singletons, and sampled diversity), CSV files with the predicted diversity trajectories for the test set and for the empirical dataset, and a plot of the estimated diversity trajectory.
