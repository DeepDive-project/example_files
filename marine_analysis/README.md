# DeepDive 


---
## Running DeepDive
Running DeepDive will typically involve four steps:

1. Preparing the input data
2. Simulating training datasets
3. training a DD model
4. running diversity predictions on the empirical data

The first step is performed in R, while steps 2-4 are performed by launching Python scripts after adjusting few settings as described below. The full pipeline has currently only been tested on Unix systems (MacOS or Linux). 

The example below, shows how to run these steps for the example of the marine empirical dataset. Note that the output of the intermediate steps (i.e. pre-formatted input files, simulated datasets, trained models) are also provided in compressed folders (which need to be unzipped to be used in DeepDive).


### 1. Prepare input data
The input data can be formatted for use in DeepDive using the R script `1.data_pipeline_marine.R`, where settings can be adjusted to change time binning of the data, age assignment method, the number of age assignment replicates and the taxonomic level of the analysis in the lines:

```
setwd("your_path/example_files") # <- change this to your path to the 'example_files' directory
bins_scale = "stages" # time bins set to stages, epochs or equal_bins (latter gives 100 equal bins)
age_method = "random_by_loc" # in each replicate randomize the age of occurrences by locality to reflect their stratigraphic age range
replicates = 10
begin_bins = 259.51
end_bins = 192.9
taxonomic_level = "Genus"
```

In this example, the only setting that must be changed is the working directory using (`setwd()`).
Running the script in the R console will generate the age randomisations and save the occurrence datasets, locality counts and time bins required as input for DeepDive in '*.csv' format. All outputs will be saved in a directory named `deepdive_data` and placed in `example_files/marine_analysis/marine_deepdive_data/`. These files will be used in the prediction script (step 4).


### 2. Simulate training data
The number of simulations and number of CPUs used to generate them can be specified in the runner file `2.simulate_marine.py`, editing the following lines: 

```
n_CPUS = 3  # number of CPUs used for simulations
n_training_simulations = 5  # simulations per CPU (the total should be ~10,000)
n_test_simulations = 10  # simulations per CPU (the total should be e.g. 100 or 1000)
training_seed = 1234 # rnd seed for training set
test_seed = 4321 # rnd seed for test set
outname = "marine"
```

After adjusting the settings as desired the script can be launched from a terminal window. 

First, activate the virtual environment:

`source path_to_virtual_env/dd_env/bin/activate`

Then browse to the marine_analysis folder

`cd your_path/example_files/marine_analysis`

Finally launch the script (note that you might have to use `python3` or `py` instead of `python` depending on your OS and settings). 

`python 2.simulate_marine.py`

This will run the simulations and generate training and test datasets including features (the pre-processed fossil data) and labels (the true diversity trajectories). These files are saved in the compressed Python format `*.npy` and can be used in the next steps to train a DeepDive model. 

Note that 1000 simulated datasets are already provided in the `example_files/simulation` folder, and can be used to test the following steps (however to train properly a model you will need a larger training set). 


### 3. Train a model
The script `3.deepdive_model_training.py` shows how to train a set of 12 models with different architecture (number of layers and nodes) based on input data generated in step 2.  
You can edit the following lines at the beginning of the script to adjust some of the settings and specify: 

1. the directory in which the input files are found, the name of the input files (features and labels in `*.npy` format)  
2. the path and folder where you want to save the trained models  
3. whether you want to parallelize the training of all models (this will need at least 12 CPUs, but more will be used if available):  

```
feature_file = "marine_test_features.npy"
label_file = "marine_test_labels.npy"
wd = "./simulations"
model_wd = "./model_marine"
parallelize_all_models = False # if True all models are trained in parallel

```

The script can be launched as shown in step 2 (note that you might have to use `python3` or `py` instead of `python` depending on your OS and settings). 

```
python 3.deepdive_model_training.py `
```

The output includes the trained models in Tensorflow format (one directory for each model), Python files with training history and settings (in pickle format, `*.pkl`), and PDF files with the plot of the training and validation accuracies computed at each epoch ('step') of the training. 

Note that pre-trained models (with training based on a large training set) are provided to help running the following steps. 


### 4. Predict diversity curves
Lastly, you can use the trained models to make diversity predictions for the empirical datasets prepared in step 1. This is done using the script `4.predict_marine.py`. 
The script will also run the trained models using the test set, to assess the expected prediction error (mean squared error on the test set) and plot the simulated features against the empirical ones (to check that the simulations were indeed able to cover the range of patterns observed in the real data)
Before running the script, you can edit the following lines to adjust some of the settings and specify: 

1. the number of randomized datasets as created in step 1 (where each replicate will have randomized fossil occurrence ages to account for their age uncertainties)
2. the number of predictions performed on each dataset
3. the correct paths and file names for the empirical data and test set simulated data:

```
n_predictions = 5  # number of predictions per input file
replicates = 10  # number of randomized input files used in data_pipeline.R

data_wd = './marine_deepdive_data' # path to empirical data
testset_wd = "./simulations" # path to the simulated test set data
model_wd = "./model_marine" # path to the trained models
output_wd = "./output" # path to directory where the putput will be saved

feature_file = "marine_test_features.npy" # testset files to calculate prediction error
label_file = "marine_test_labels.npy"
res_file = "marine_results" # name used to save the outputs. 
```

The script is launched as shown in steps 2 and 3 and generates three output files: a PDF file with the estimated diversity trajectory for the empirical clade, a CSV table with the training, validation and test MSE for each model, and a PNG file with plots summarizing the simulated and empirical features.




