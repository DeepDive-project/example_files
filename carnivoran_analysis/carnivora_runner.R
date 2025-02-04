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
# For now, we will use the carnivoran occurrences included in the package
data(carnivora)
dat <- carnivora

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
# with the argument 'r=100' we are randomly resampling occurrence ages 
# 10 times from their stratigraphic ranges
prep_dd_input(dat = dat, bins = bins, r = 100, output_file = dd_file_name)


# Create a config for the DeepDive analysis
# If applicable you can specify the number of living taxa which will be used 
# by the model to calibrate the predicted diversity trajectories
config <- create_config(
      bins = bins,
      name="carnivora",
      n_regions = length(unique(dat$Region)),
      present_diversity = 313,
      data_file = dd_file_name
)


# To have regions becoming available over time, create area_ages object
# Areas are expected to be in alphabetical order
# Here for example we assume that dispersal to South America only becomes  
# possible between 11 and 7 Ma
region_ages <- rbind(c("Africa", max(bins), max(bins)),  
                     c("Asia", max(bins), max(bins)),  
                     c("Europe", max(bins), max(bins)), 
                     c("North America", max(bins), max(bins)),  
                     c("South America", 11.608, 7.3))            
region_ages <- as.data.frame(region_ages)
# Label columns
colnames(region_ages) <- c("Region", "MaxAge", "MinAge")

# Regions disappearing instead of connecting can be made via setting the label argument
# to: label = "end"
regions_matrix(config, region_ages, presence = TRUE)

# add models 
add_model(config=config, lstm_nodes = c(64, 32), dense_nodes = c(64, 32), model_name = "1")
add_model(config=config, lstm_nodes = c(128, 64, 32), dense_nodes = c(64, 32), model_name = "2")
add_model(config=config, lstm_nodes = c(256, 128, 64), dense_nodes = c(64, 32), model_name = "3")
add_model(config=config, lstm_nodes = c(512, 128), dense_nodes = c(64, 32), model_name = "4")


# write the configuration file
config$write("carnivora.ini")
