library(devtools)
load_all(".")
library(DeepDiveR)
library(readxl)
library(data.table)


# load an occurrence table here
# your_dat <- read.csv("your_file_path")
path_dat <- "C:/Users/CooperR/Documents/GitHub/carnivore_analysis"
dat <- read_xlsx(paste0(path_dat, "/carnivora_data_cleaned.xlsx"), sheet="Records")

dat <- dat[-which(dat$Maximum_Age =="Unknown"),]
dat <- dat[-which(dat$Latitude == "New_Record"),]
dat <- dat[-which(dat$Species == "Genus_Only"),]
dat$Complete_name <- paste(dat$Genus,dat$Species)

dat <- data.frame(dat$Complete_name, dat$Continent, 
                  as.numeric(dat$Minimum_Age), as.numeric(dat$Maximum_Age), 
                  dat$Latitude, dat$Longitude)
colnames(dat) <- c('Taxon', 'Area', 'MinAge', 'MaxAge', 
                   'Latitude', 'Longitude')

# assign localities from coordinates
dat$Latitude <- round(as.numeric(dat$Latitude), digits=1)
dat$Longitude <- round(as.numeric(dat$Longitude), digits=1)
setDT(dat)[, Locality := .GRP, by =.(Latitude, Longitude, MinAge, MaxAge)]

dat <- dat[!duplicated(dat), ]  # remove duplicated rows


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

# Check the number of extant taxa
summary_dat <- read_xlsx(paste0(path_dat, "/carnivora_data_cleaned.xlsx"), sheet="Species summary")
summary_dat <- summary_dat[-which(summary_dat$Minimum_age =="Unknown"),]
summary_dat <- summary_dat[-which(summary_dat$Status == "Extinct"),]
summary_dat$Complete_name <- paste(summary_dat$Genus,summary_dat$Species)
summary_dat <- summary_dat[1:313,]  ## remove row with NA
summary_dat <- data.frame(summary_dat$Complete_name, summary_dat$Status)
colnames(summary_dat) <- c("Taxon", "Status")
summary_dat <- summary_dat[!duplicated(summary_dat), ]
length(summary_dat$Status == "Extant")


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


# From the data, find when areas must have been occupied by (MaxAge for oldest
# fossil sampled per continent, in this case).
area_tables <- split(dat, f = dat$Area)  # Split data by area
Africa <- max(area_tables$Africa$MaxAge)
Asia <- max(area_tables$Asia$MaxAge)
Europe <- max(area_tables$Europe$MaxAge)
N_America <- max(area_tables$NorthAmerica$MaxAge)
S_America <- max(area_tables$SouthAmerica$MaxAge)

# to have regions becoming available over time, create area_ages object
# should be written in alphabetical order, give a label and sort them
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
