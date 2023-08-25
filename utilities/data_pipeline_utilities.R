# Functions for organising data in the correct format for deepdive
# Run functions in data_pipeline.R
library(dplyr)
library(readxl)
library(data.table)
library(tidyr)


create_folders <- function(name){
  dir.create(paste0(name, "deepdive_data"))
  folder_names <- c("Species_occurrences", "Genus_occurrences", "Locality", "Info")
  for (i in 1:length(folder_names)){
    dir.create(paste0(name, "deepdive_data/", folder_names[i]))
  }
}


time_bins <- function(dat, scale, begin, finish, n_bins=100, remove_bins=NULL, 
                      use_q=FALSE, merge_holo=FALSE, holo_together= FALSE, geochart = NULL){
  if(scale == "equal_bins"){  # add remove holocene and extend pleistocene to equal_bins
    bins <- build_bins(start = begin, end = finish, n = n_bins) 
  }
  if(scale == "epochs"){
    stages <- read_xlsx(geochart)
    bins <- data.frame(stages$series, stages$bottom, stages$top, stages$mid)
    colnames(bins) <- c('epoch', 'start', 'end', 'midpoint')
    bins <- bins[which(stages$bottom<=begin), ]
    bins <- bins[which(bins$end>=finish), ]
    epoch_tables <- split(bins, f = bins$epoch)
    new_bins <- data.frame()
    for(i in 1:length(epoch_tables)){
      epoch <- unique(epoch_tables[[i]][,1])
      start <- max(epoch_tables[[i]][,2])
      end <- min(epoch_tables[[i]][,3])
      midpoint <- (start-end)/2 + end
      e <- cbind(epoch, start, end, midpoint)
      new_bins <- rbind(new_bins, e)
    }
    bins <- new_bins[order(new_bins$start, decreasing = TRUE), ]
    if(!is.null(remove_bins)){
      bins <- bins[which(bins$stage!=remove_bins), ]
      # bins$end[length(bins$end)]<-0
    }
    if(use_q == T & finish==0){
      bins<-bins[1:((nrow(bins))-7),]
      bins <- rbind(bins, c("Quaternary", 2.58, 0, 1.29))
    }
  }
  if(scale == "stages"){
    stages <- read_xlsx(geochart)
    bins <- data.frame(stages$stage, stages$bottom, stages$top, stages$mid)
    colnames(bins) <- c('stage', 'start', 'end', 'midpoint')
    bins <- bins[which(stages$bottom<=begin), ]
    bins <- bins[which(bins$end>=finish), ]
    if(!is.null(remove_bins)){
      bins <- bins[which(bins$stage!=remove_bins), ]
      # bins$end[length(bins$end)]<-0
    }
    if(use_q == T & finish==0){
      bins<-bins[1:((nrow(bins))-7),]
      bins <- rbind(bins, c("Quaternary", 2.58, 0, 1.29))
    }
    if(merge_holo == T & finish<=0.0117){
      n_holo_bins_merge <- length(which(bins$start<=0.129))
      bins<-bins[1:((nrow(bins))-n_holo_bins_merge),]
      midpoint <- (0.129-finish)/2
      bins <- rbind(bins, c("Upper", 0.129, finish, midpoint))
    }
    if(holo_together == T & finish<=0.0117){
      n_holo_bins_merge <- length(which(bins$start<=0.0117))
      bins<-bins[1:((nrow(bins))-n_holo_bins_merge),]
      midpoint <- (0.0117-finish)/2
      bins <- rbind(bins, c("Holocene", 0.0117, finish, midpoint))#0.0645))
    }
  }    
  return(bins)
}



# Build time bins function (adapted from master's project code)
build_bins <- function (start, end, n = NULL) {
  # convert number of bins to bin length
  interval <- abs((start - end) / n)
  
  if (start > end) {
    interval <- -abs(interval)
  }
  
  # vector of bin boundary times
  bin_bound_seq <- seq(start, end, interval)
  
  # offset start and end times by 1 to build bins of required length
  bins <- tibble(
    start = bin_bound_seq[-length(bin_bound_seq)],
    end   = bin_bound_seq[-1]
  )
  
  bins$midpoint <- rowMeans(bins)
  bins
}


# Get ages function 
ages <- function(dat, method, area_tables){
  if (method == "median") {
    dat <- mutate(rowwise(dat), SampledAge = median(c(MinAge, MaxAge)))
  }
  if (method == "random") {
    SampledAge <- runif(length(dat$Complete_name), min = dat$MinAge, max = dat$MaxAge)
    dat <- cbind(dat, SampledAge)
  } 
  if (method == "random_by_loc") {
    dat$SampledAge <- dat$Locality
    for (i in 1:length(unique(dat$Locality))) {
      locate <- which(dat$Locality == i)
      loc <- dat[locate, ]
      Age <- runif(n = 1, min = loc$MinAge, max = loc$MaxAge)
      dat$SampledAge[dat$SampledAge == i] <- Age
    }
  }
  return(dat)
}


get_lr_hr_dat <- function(data=dat, res="high", age_range_threshold = NA){
  if(!is.na(age_range_threshold == TRUE)){
    print("Using age_range_threshold")
    age_range <- dat$MaxAge - dat$MinAge  # find age ranges
    
    # find low and high res data and their indices
    lr <- dat[which(age_range > age_range_threshold), ]
    lr_ind <- which(age_range > age_range_threshold)
    hr <- dat[which(age_range <= age_range_threshold), ]
    hr_ind <- which(age_range <= age_range_threshold)
    
    if(res == "low"){
      # make a low res data subset of same dimensions as original data
      r_dat <- dat
      for(i in r_dat) {
        r_dat[lr_ind, ] <- lr
        r_dat[hr_ind, 5:6] <- NA
      }
      colnames(r_dat) <- c('Complete_name', 'Genus', 'Species', 'Area', 'MinAge', 'MaxAge', 
                           'Locality')
      
    }
    
    if(res == "high"){
      # make a high res data subset of same dimensions as original data
      r_dat <- dat
      for(i in r_dat) {
        r_dat[lr_ind, 5:6] <- NA
        r_dat[hr_ind, ] <- hr
      }
      colnames(r_dat) <- c('Complete_name', 'Genus', 'Species', 'Area', 'MinAge', 'MaxAge', 
                           'Locality')
      
    }
  }
  else{
    print("Using assignment to lr or hr bins")
    hr <- dat_res[which(dat_res$resolution == "HR"),]
    lr <- dat_res[which(dat_res$resolution == "LR"),]
    hr_ind <- which(dat_res$resolution == "HR")
    lr_ind <- which(dat_res$resolution == "LR")
    
    if(res == "low"){
      # make a low res data subset of same dimensions as original data
      r_dat <- dat_res
      for(i in r_dat) {
        r_dat[lr_ind, ] <- lr
        r_dat[hr_ind, 5:6] <- NA
      }
      colnames(r_dat) <- c('Complete_name', 'Genus', 'Species', 'Area', 'MinAge', 'MaxAge', 
                           'Locality', 'Resolution')
      
    }
    
    if(res == "high"){
      # make a high res data subset of same dimensions as original data
      r_dat <- dat_res
      for(i in r_dat) {
        r_dat[lr_ind, 5:6] <- NA
        r_dat[hr_ind, ] <- hr
      }
      colnames(r_dat) <- c('Complete_name', 'Genus', 'Species', 'Area', 'MinAge', 'MaxAge', 
                           'Locality', 'Resolution')
      
    }
    
    
  }
  return(r_dat)
}


fraction_lr_hr_dat <- function(dat, res="low", r_dat=NA, age_range_threshold = NA){
  # check fraction of occurrences assigned to the high res and low res datasets
  if(!is.na(age_range_threshold == TRUE)){
    print("Using age_range_threshold")
    age_range <- dat$MaxAge-dat$MinAge 
    if(res == "high"){
      fraction_res <- nrow(dat[which(age_range <= age_range_threshold), ]) / (nrow(dat[which(age_range > age_range_threshold), ]) + nrow(dat[which(age_range <= age_range_threshold), ]))
    }
    if(res == "low"){
      fraction_res <- nrow(dat[which(age_range > age_range_threshold), ]) / (nrow(dat[which(age_range > age_range_threshold), ]) + nrow(dat[which(age_range <= age_range_threshold), ]))
    }
  }
  else{
    print("Using assignment to lr or hr bins")
    fraction_res <- length(r_dat$MinAge)/length(dat$MinAge)
  }
  return(fraction_res)
}


# TAXON X AREA X TIME TABLE FOR EITHER TAXONOMIC LEVEL 
# taxon time per area table, needs to be run per area
taxa_time_per_area <- function(dat, area_tables, bins, taxon_level, ind, name){
  for(i in 1:length(area_tables)){
    area_table <- area_tables[[i]]
    area_name <- area_table$Area[1]
    if (taxon_level == "Species"){
      area_taxa <- unique(area_table$Complete_name)
      global_taxa_list <- unique(dat$Complete_name)
    }
    if (taxon_level == "Genus"){
      area_taxa <- unique(area_table$Genus)
      global_taxa_list <- unique(dat$Genus)
    }
    all_taxa <- length(global_taxa_list)
    taxa_time_table <- data.frame(matrix(0, all_taxa, nrow(bins)))
    bins_0 <- -c(as.numeric(bins$start), 0)  # add zero to time bins (negative value time)
    for (j in 1:length(area_taxa)){
      if (taxon_level == "Species") {
        indices_occurrences <- which(area_table$Complete_name == area_taxa[j])
      }
      if (taxon_level == "Genus") {
        indices_occurrences <- which(area_table$Genus == area_taxa[j])
      }
      age_occs <- area_table[indices_occurrences,]$SampledAge
      h <- hist(x = -as.numeric(age_occs), breaks=bins_0, plot=F)
      if (taxon_level == "Species"){
        t_i <- which(global_taxa_list == area_taxa[j])
      }
      if (taxon_level == "Genus"){
        t_i <- which(global_taxa_list == area_taxa[j])
      }
      taxa_time_table[t_i,] <- h$counts
    }
    write.csv(taxa_time_table, paste0(name, "deepdive_data/", taxon_level, 
                                      "_occurrences/", ind, "_", area_name, 
                                      "_occs.csv"), row.names=FALSE)
  }
}


# Locality data set (localities per area and time bin)
generate_locality_dataset <- function(dat, bins, ind, name){
  list_areas <- unique(dat$Area)
  localities <- data.frame(matrix(0, length(list_areas), nrow(bins)))
  bins_0 <- -c(as.numeric(bins$start), 0) # add zero to bins (using negative values for time)
  for (i in seq_len(length(list_areas))){      
    indices_areas <- which(dat$Area == list_areas[i]) 
    locality_ids <- dat[indices_areas,]$Locality
    uni_loc_ids <- unique(locality_ids)
    no_loc_in_area <- c()
    for(j in uni_loc_ids){
      t <- dat$SampledAge[which(dat$Locality == j)]
      no_loc_in_area <- c(no_loc_in_area, unique(t))
    }
    h <- hist(x = -as.numeric(no_loc_in_area), breaks=bins_0, plot=F)
    localities[i,] <- h$counts
  }
  write.csv(localities, paste0(name, "deepdive_data/Locality/", ind, "_localities.csv"), row.names=FALSE)
}


# occurrences per continent and time bin
generate_occurrence_dataset <- function(dat, area_tables, bins, res, name){
  list_areas <- unique(dat$Area)
  n_occurrences <- data.frame(matrix(0, length(list_areas), nrow(bins)))
  bins_0 <- -c(as.numeric(bins$start), 0) # add zero to bins (using negative values for time)
  for (i in seq_len(length(list_areas))){
    indices_areas <- which(dat$Area == list_areas[i]) 
    total_occs_for_area <- length(indices_areas)
    area_dat <- dat[indices_areas,]
    age_occs <- area_dat$SampledAge
    # uni_age_occs <- unique(age_occs)
    h <- hist(x = -as.numeric(age_occs), breaks=bins_0, plot=F)
    n_occurrences[i,] <- h$counts
  }
  write.csv(n_occurrences, paste0(name, "deepdive_data/Info/", res, "_occurrences.csv"), row.names=FALSE)
}

write_dd_files <- function(data, res="hr", age_m="median", taxon_level, bins){
  # Get ages and append using median
  sampled_dat <- ages(data, method=age_m)
  # sampled_dat <- sampled_dat[-which(sampled_dat$SampledAge>259.1),]
  saveRDS(sampled_dat, paste0(name, "deepdive_data/Info/", res, "dat_with_ages.RDS"))
  
  area_tables <- split(sampled_dat, f = sampled_dat$Area)  # Split data by area 
  
  # Get species or genera level data
  taxa_time_per_area(sampled_dat, area_tables, bins=bins, taxon_level=taxon_level, res, name)
  
  # Get locality data
  generate_locality_dataset(sampled_dat, bins=bins, res, name)
  locs <- read.csv(paste0(name, "deepdive_data/Locality/", res, "_localities.csv"))
  
  # Get occurrence data
  #generate_occurrence_dataset(sampled_dat, area_tables, bins=bins, res, name)
  #occs <- read.csv(paste0(name, "deepdive_data/Info/", res, "_occurrences.csv"))
  
  # Occurrences per locality through time for each area
  #occs_per_loc <- occs/(locs+1)
  #write.csv(occs_per_loc, paste0(name, "deepdive_data/Info/", res, "_occurrences_per_locality.csv"), row.names=FALSE)
  
}
