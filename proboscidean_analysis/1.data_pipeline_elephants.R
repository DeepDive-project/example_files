# Data pipeline Elephants
library(dplyr)
library(readxl)
library(data.table)
library(tidyr)
library(stringr)

setwd("your_path/deepdive_for_review") # <- change this to your path 
source("utilities/data_pipeline_utilities.R")
input_occs_file <- "Proboscidea_occurrences_def_Nov.2020.xlsx"
name <- "elephant_analysis/"
create_folders(name)


# Settings
bins_scale = "equal_bins"  # time bins set to stages, epochs or equal_bins (latter gives 100 equal bins)
age_method = "random_by_loc"  # can be set to random_by_loc or random, for median ages use alternative pipeline
replicates = 100  # number of age randomisation replicates
begin_bins = 66
end_bins = 3.6
age_range_threshold = NA  # e.g. = 2 occs with age range > 2 Myr are treated as low resolution, if NA assigns whether max_age-min_age is within the high or low resolution bins set
taxonomic_level = "Species"

# Read in the data
dat <- read_xlsx(input_occs_file)
dat <- data.frame(dat$species_complete_corrected, dat$species_complete_corrected,
                  dat$continent, dat$MIN_AGE_homog, dat$MAX_AGE_homog, dat$NAME)
colnames(dat) <- c('Complete_name', 'Genus_species', 'Area', 'MinAge', 'MaxAge', 
                   'Locality')

dat <- dat %>% separate(Genus_species, c('Genus', 'Species'), sep=" ")

# re-code geography
setDT(dat)[, Locality := .GRP, by = Locality]
dat$Area <- str_replace_all(dat$Area, c("Northern Africa" = "Africa",
                                        "Eastern Africa" = "Africa",
                                        "Western Africa" = "Africa",
                                        "Southern Africa" = "Africa",
                                        "Middle Africa" = "Africa",
                                        "Eastern Asia" = "Asia", 
                                        "South-Eastern Asia" = "Asia",
                                        "Western Asia" = "Asia",
                                        "Southern Asia" = "Asia",
                                        "South-Asia" = "Asia",
                                        "Central Asia" = "Asia",
                                        "Northern Europe" = "Europe", 
                                        "Eastern Europe" = "Europe", 
                                        "Western Europe" = "Europe", 
                                        "Southern Europe" = "Europe",
                                        "Northern America" = "North America",
                                        "Central America" = "North America"))

# check some stats
dat <- dat[!duplicated(dat), ]

# low res bins, stages
lr_bins <- time_bins(dat, scale="stages", begin = begin_bins, finish=0, merge_holo = T, 
                     geochart=paste0(getwd(), "/utilities/geochart.xlsx"))
lr_bins_0 <- -as.numeric(lr_bins$start, 0)
write.csv(diff(lr_bins_0), paste0(name, "deepdive_data/t_bins.csv"), row.names=FALSE)


# Build higher res time bins, with boundaries that align to lr_bins
# with length 1 My on average
# find start and end of lr bins with finer bins after 0.774
low_res_start <- c(as.numeric(lr_bins$start)[-length(lr_bins$start)])
low_res_end <- c(as.numeric(lr_bins$end)[-length(lr_bins$end)])
n <- round(low_res_start-low_res_end) # = 3
start <- c()
lr_bin_indices <- c()
for(i in 1:length(low_res_start)){
  high_res <- seq(low_res_start[i], low_res_end[i], length.out=n[i]+1)
  high_res <- high_res[-length(high_res)]
  start <- append(start, high_res)
  lr_bin_indices <- append(lr_bin_indices, rep(i, times=length(high_res)))
}
start <- append(start, c(0.129, 0.0117, 0.0082, 0.0042))
end <- append(start[-1], 0)
midpoint <- (start-end)/2 + end
bins <- data.frame(start, end, midpoint)
avg_bin_length <- mean(start-end)
range_bin_lengths <- range(start-end)
lr_bin_indices <- append(lr_bin_indices, tail(lr_bin_indices, n=1) + 1:length(c(0.129, 0.0117, 0.0082, 0.0042)))
bins_0 <- -c(bins$start, 0)
write.csv(diff(bins_0), paste0(name, "/t_bins.csv"), row.names = FALSE)
saveRDS(bins, paste0(name, "/bins.RDS"))
write.table(lr_bin_indices, paste0(name, "bin_indices.txt")) 


for(i in 1:nrow(dat)){
  if (dat$MinAge[i] < 0.004 & dat$Complete_name[i] == "Stegodon zhaotongensis"){
    dat$MinAge[i] <- 0.0117
  }   
  if (dat$MinAge[i] < 0.004 & dat$Complete_name[i] == "Elephas hysudrindicus"){
    dat$MinAge[i] <- 0.015
  } 
  if (dat$MinAge[i] < 0.004 & dat$Complete_name[i] == "Mammut americanum"){
    dat$MinAge[i] <- 0.011
  } 
}


# Get ages and append using either median, random or random_by_loc
# If using random ages repeat the random draws in the for loop
for(replicate in 1:replicates){
  write_dd_files(dat, res=replicate, taxon_level = taxonomic_level, bins=bins, age_m=age_method)
}
