# Data pipeline Permo-Triassic marine animals
library(dplyr)
library(readxl)
library(data.table)
library(tidyr)
library(stringr)

#------ Script settings:
setwd("your_path_to_example_files") # <- change this to your path to the 'example_files' directory
bins_scale = "stages" # time bins set to stages, epochs or equal_bins (latter gives 100 equal bins)
age_method = "random_by_loc" # in each replicate randomize the age of occurrences by locality to reflect their stratigraphic age range
replicates = 10 
begin_bins = 259.51
end_bins = 192.9
taxonomic_level = "Genus"
#------


source("utilities/data_pipeline_utilities.R")
input_occs_file <- "marine_analysis/marine_deepdive_data/NC_raw.csv" # raw data
name <- "marine_analysis/marine_deepdive_data/"
create_folders(name)

dat <- read.csv(input_occs_file)
dat <- dat[dat$max_ma <= begin_bins,]

# Re-code geography
setDT(dat)[, Loc_ID:=.GRP, by=.(lng, lat, max_ma, min_ma)]
dat <- data.frame(dat$genus, dat$genus, dat$genus, dat$reg, dat$min_ma, dat$max_ma, dat$Loc_ID)
colnames(dat) <- c('Complete_name', 'Genus', 'Species', 'Area', 'MinAge', 'MaxAge', 'Locality')
dat <- dat[!duplicated(dat), ]
dat$Complete_name <- NA
dat$Species <- NA

# Build low resolution time bins
bins <- time_bins(dat, scale=bins_scale, begin=begin_bins, finish=end_bins, use_q=T,
    geochart="utilities/geochart.xlsx")
bins_0 <- -c(bins$start, 192.9)
write.csv(diff(bins_0), paste0(name, "deepdive_data/t_bins.csv"), row.names=FALSE)


# Get ages and append using either median, random or random_by_loc
# If using random ages repeat the random draws in the for loop
for(replicate in 1:replicates){
    write_dd_files(dat, res=replicate, taxon_level=taxonomic_level, bins=bins, age_m=age_method)
}
