library(devtools)
load_all(".")
library(DeepDiveR)
library(readxl)
library(data.table)


# load an occurrence table here
path_dat <- "your_path"  # path to folder containing cleaned data
dat <- read_xlsx(paste0(path_dat, "/carnivora_data_cleaned.xlsx"), sheet="Records")

# Specify a vector of time bins
bins <- c(max(dat$MaxAge), 65, 64, 63, 61.6, 60, 59.2, 58.13333, 57.06667, 56, 
          54.975, 53.95, 52.925, 51.9, 50.875, 49.85, 48.825, 47.8, 46.85714, 
          45.91429, 44.97143, 44.02857, 43.08571, 42.14286, 41.2, 40.03667, 
          38.87333, 37.71, 36.7575, 35.805, 34.8525, 33.9, 32.88667, 31.87333, 
          30.86, 29.84667, 28.83333, 27.82, 26.862, 25.904, 24.946, 23.988, 
          23.03, 22.16667, 21.30333, 20.44, 19.3225, 18.205, 17.0875, 15.97, 
          14.895, 13.82, 12.725, 11.63, 10.534, 9.438, 8.342, 7.246, 6.2895, 
          5.333, 4.4665, 3.6, 2.58, 1.8, 0.774, 0.129, 0.0117, 0)

dd_input_file_name <- "deepdive_input.csv"
dd_dataset <- paste(path_dat, dd_input_file_name, sep="/")

# Prepare input file for deepdive, set output_file name to save
prep_dd_input(dat = dat, bins = bins, r = 100, 
              age_m = "random_by_loc", output_file = dd_dataset)


# create a config for a full analysis
config <- create_config(
  path_wd = path_dat,  
  bins = bins,
  sim_name="carnivora",
  n_areas = length(unique(dat$Area)),
  simulations_file = "simulations_carnivora", 
  models_file = "trained_models_carnivora", 
  add_test = T, 
  present_diversity = 313,
  autotune=TRUE,
  empirical_input_file = dd_input_file_name
)

# edit number of living taxa, present_diversity is 313
set_value(attribute_name = "extant_sp", value=c(31, 3100), module="simulations", config)

# edit total number of simulated taxa (minimum and maximum) 
set_value(attribute_name = "total_sp", value=c(618, 2000), module="simulations", config)

# edit carrying capacity in equilibrium simulations
set_value(attribute_name="dd_K", value=c(31, 3100), module="simulations", config)

set_value(attribute_name = "pr_extant_clade", value=1, module="simulations", config)


# To have regions becoming available over time, create area_ages object
# should be written in alphabetical order
area_ages <- rbind(c(61.6, 59.2),  # Africa 
                   c(max(bins), max(bins)),  # Asia
                   c(59.2, 56),  # Europe
                   c(max(bins), max(bins)),  # North America
                   c(11.608, 7.3))  # South America 


# if entered as area_ages = NULL, also provide bins as an argument.
# regions disappearing instead of connecting can be made via adding the argument
# label = "end"
areas_matrix(area_ages, n_areas = length(unique(dat$Area)), config)

# write the configuration file
config$write(paste(path_dat, "carnivora.ini", sep="/"))
