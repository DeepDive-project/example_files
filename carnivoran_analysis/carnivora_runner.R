library(devtools)
# Option 1: install the package from GitHub
devtools::install_github("DeepDive-project/DeepDiveR")
# Option 2: load it from a directory after donwloading it from
# https://github.com/DeepDive-project/DeepDiveR
deepdiver_path <- "path_to_DeepDiveR"
setwd(deepdiver_path)
load_all(".")
library(DeepDiveR)


# Load an occurrence table here
# This typically includes one row for each occurrence and 5 columns:
# taxon name, geographic area, min age, max age, and locality
path_dat <- "your_path"  # path to the input data 
setwd(path_dat)
dat <- read.csv("carnivora_data.csv")

# Specify a vector of time bins (they do not need to be equally spaced)
bins <- c(max(dat$MaxAge), 65, 64, 63, 61.6, 60, 59.2, 58.13333, 57.06667, 56, 
          54.975, 53.95, 52.925, 51.9, 50.875, 49.85, 48.825, 47.8, 46.85714, 
          45.91429, 44.97143, 44.02857, 43.08571, 42.14286, 41.2, 40.03667, 
          38.87333, 37.71, 36.7575, 35.805, 34.8525, 33.9, 32.88667, 31.87333, 
          30.86, 29.84667, 28.83333, 27.82, 26.862, 25.904, 24.946, 23.988, 
          23.03, 22.16667, 21.30333, 20.44, 19.3225, 18.205, 17.0875, 15.97, 
          14.895, 13.82, 12.725, 11.63, 10.534, 9.438, 8.342, 7.246, 6.2895, 
          5.333, 4.4665, 3.6, 2.58, 1.8, 0.774, 0.129, 0)

dd_file_name <- "carnivora_deepdive_input.csv"

# Prepare input file for deepdive, set output_file name to save
# with the argument 'r=10' we are randomly resampling occurrence ages 
# 10 times from their stratigraphic ranges
prep_dd_input(dat = dat, bins = bins, r = 10, output_file = dd_file_name)


# Create a config for the DeepDive analysis
# If applicable you can specify the number of living taxa which will be used 
# by the model to calibrate the predicted diversity trajectories
config <- create_config(
      path_wd = path_dat,  
      bins = bins,
      sim_name="carnivora",
      n_areas = length(unique(dat$Area)),
      simulations_file = "simulations_carnivora", 
      models_file = "trained_models_carnivora", 
      present_diversity = 313,
      empirical_input_file = dd_file_name
)


# To have regions becoming available over time, create area_ages object
# Areas are expected to be in alphabetical order
# Here for example we assume that dispersal to South America only becomes  
# possible between 11 and 7 Ma
area_ages <- rbind(c(max(bins), max(bins)),  # Africa 
                   c(max(bins), max(bins)),  # Asia
                   c(max(bins), max(bins)),  # Europe
                   c(max(bins), max(bins)),  # North America
                   c(11.608, 7.3))           # South America 


# Regions disappearing instead of connecting can be made via adding the argument
# label = "end"
areas_matrix(area_ages, n_areas = length(unique(dat$Area)), config)

# For the purpose of runnign a quick test analysis we use low number of training 
# and test simulations 
# Note that you will need at least 10,000 training simulations for a real analysis
set_value(attribute_name = "n_training_simulations", 
          value=100, 
          module="simulations", config)

set_value(attribute_name = "n_test_simulations", 
        value=10, 
        module="simulations", config)


# write the configuration file
config$write("carnivora.ini")
