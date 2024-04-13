import os
import sys
import deepdive as dd
import numpy as np
from datetime import datetime
import pandas as pd
import multiprocessing
today = datetime.now()
now = datetime.now().strftime('%Y%m%d')
try:
    os.mkdir('simulations/')
except:
    pass

output_path = "simulations/"
outname = "elephant_sim"

# init simulation environment
n_CPUS = 3  # number of CPUs used for simulations
n_training_simulations = 10  # simulations per CPU (the total should be ~10,000)
n_test_simulations = 10  # simulations per CPU (the total should be e.g. 100 or 1000)

training_seed = 1234
test_seed = 4321
max_clade_age = 100
n_time_bins = 100



#------
use_bins = True
use_area_constraints = True
if use_bins:
    time_bins = np.sort(np.array([66, 65, 64, 63, 61.6, 60, 59.2, 58.13333, 57.06667, 56, 54.975, 53.95, 52.925, 51.9,
                                  50.875, 49.85, 48.825, 47.8, 46.85714, 45.91429, 44.97143, 44.02857, 43.08571,
                                  42.14286, 41.2, 40.03667, 38.87333, 37.71, 36.7575, 35.805, 34.8525, 33.9, 32.88667,
                                  31.87333, 30.86, 29.84667, 28.83333, 27.82, 26.862, 25.904, 24.946, 23.988, 23.03,
                                  22.16667, 21.30333, 20.44, 19.3225, 18.205, 17.0875, 15.97, 14.895, 13.82, 12.725,
                                  11.63, 10.534, 9.438, 8.342, 7.246, 6.2895, 5.333, 4.4665, 3.6, 2.58, 1.8, 0.774,
                                  0.129, 0.0117, 0.0082, 0.0042, 0]))

    n_time_bins = len(time_bins) - 1


if use_area_constraints:
    # area constraints
    # data.frame (rows: n_areas, cols: [area_id, start, end]) if -1 start at the beginning or end at time 0
    # the dataframe can also be read from a file e.g. using pd.read_csv()
    n_areas = 5
    area_tbl = np.ones((n_areas, 3))
    area_tbl[:, 0] = np.arange(n_areas)  # set areas IDs
    area_tbl[:, 1:] = -1  # set all values to -1 (all areas exist throughout)
    eurasia_age = np.random.uniform(27, 33.9)  # Cantalapiedra + Eocene/Oligocene boundary Eurasia appears
    area_tbl[1, 1] = eurasia_age
    area_tbl[2, 1] = eurasia_age  # repeat as analysis is splitting Europe and Asia
    area_tbl[3, 1] = np.random.uniform(16, 20)  # Cantalapiedra, North America appears
    area_tbl[4, 1] = np.random.uniform(0.8, 5.3)  # Carrillo et al. and Cantalapiedra's data South America appears
    area_tbl = pd.DataFrame(area_tbl)

    time_bins_duration = np.diff(time_bins)
    mid_time_bins = []
    for i in range(0, n_time_bins):
        mid_time_i = time_bins[i] + 0.5 * time_bins_duration[i]
        mid_time_bins.append(mid_time_i)
    mid_time_bins = np.array(mid_time_bins)

#------


# create simulator object
def create_sim_obj(rseed):
    bd_sim = dd.bd_simulator(s_species=1,  # number of starting species
                             rangeSP=[200, 2000],  # min/max size data set
                             minEX_SP=0,  # minimum number of extinct lineages allowed
                             root_r=[50., max_clade_age],  # range root ages
                             minEXTANT_SP=3,  # min number of living species
                             maxEXTANT_SP=30,
                             rangeL=[0.01, 0.3],  # range of birth rates
                             rangeM=[0.01, 0.3],  # range of death rates
                             p_mass_extinction=0.01,  # probability of mass extinction per my
                             poiL=4,  # expected number of birth rate shifts
                             poiM=4,  # expected number of death rate shifts
                             seed=rseed,
                             scale=10,
                             p_equilibrium=0.5,
                             vectorize=True)  # if > 0 fixes the random seed to make simulations reproducible

    # create fossil simulator object
    fossil_sim = dd.fossil_simulator(n_areas=n_areas,
                                     n_bins=n_time_bins,  # number of time bins
                                     time_bins=time_bins,
                                     eta=[1, 1.75],  # area-sp stochasticity - previously 1, 5
                                     p_gap=[0.01, 0.95],  # probability of 0 preservation in a time bin
                                     dispersal_rate=None,
                                     max_dist=1,
                                     disp_rate_mean=[0, 1],
                                     disp_rate_variance=1,
                                     area_mean=50,  # G(a,b) distributed preservation rates across areas, previously 50
                                     area_variance=1,  # previously 5
                                     size_concentration_parameter=[0.1, 3],  # single value or array of length n_areas
                                     link_area_size_carrying_capacity=[1, 10],
                                     # positive, larger numbers = stronger link between area size and carrying capacity
                                     p_origination_a_slope_mean=2,  # mean slope of probability of origination area mean
                                     p_origination_a_slope_sd=0.5,  # std of the slopes
                                     sp_mean=[0.3, 1],  # G(a,b) distributed preservation rates across species
                                     sp_variance=2,
                                     slope=[-0.01, 0],  # change in log-sampling rate through time (log-linear)
                                     intercept=[0.1, 0.5],  # initial sampling rate
                                     sd_through_time=[0.001, 0.01],  # st dev in log-sampling rate through time
                                     sd_through_time_skyline=1,
                                     mean_n_epochs_skyline=4,
                                     fraction_skyline_sampling=0.5,
                                     maximum_localities_per_bin=300,
                                     singletons_frequency=[0, 0.5],
                                     seed=rseed)  # if > 0 fixes the random seed to make simulations reproducible
    return bd_sim, fossil_sim


def run_sim(rep):
    batch_features = []
    batch_labels = []

    # simulate training data
    bd_sim, fossil_sim = create_sim_obj(training_seed + rep)
    for i in range(n_training_simulations):
        if i % 1 == 0 and rep == 0:
            dd.print_update("%s of %s done" % (i+1, n_training_simulations))
        sp_x = bd_sim.run_simulation(print_res=False)

        #----
        if use_area_constraints:
            c1, c2 = dd.set_area_constraints(sp_x=sp_x,
                                             n_time_bins=n_time_bins,
                                             area_tbl=area_tbl,
                                             mid_time_bins=mid_time_bins)

            fossil_sim.set_carrying_capacity_multiplier(m_species_origin=c1, m_sp_area_time=c2)
        #----

        # using min_age and max_age ensures that the time bins always span the same amount of time
        sim = fossil_sim.run_simulation(sp_x, min_age=0, max_age=max_clade_age, return_sqs_data=False)

        #----
        sim_features = dd.extract_sim_features(sim)
        sim_y = sim['global_true_trajectory']
        batch_features.append(sim_features)
        batch_labels.append(sim_y)

    res = {'features': batch_features, 
           'labels': batch_labels}
    return res


def run_test_sim(rep):
    batch_features = []
    batch_labels = []

    # simulate training data
    sim_settings = []
    bd_sim, fossil_sim = create_sim_obj(test_seed + rep)
    for i in range(n_test_simulations):
        if i % 1 == 0 and rep == 0:
            dd.print_update("%s of %s done" % (i+1, n_test_simulations))
        sp_x = bd_sim.run_simulation(print_res=False)
        #----
        if use_area_constraints:
            c1, c2 = dd.set_area_constraints(sp_x=sp_x,
                                             n_time_bins=n_time_bins,
                                             area_tbl=area_tbl,
                                             mid_time_bins=mid_time_bins)

            fossil_sim.set_carrying_capacity_multiplier(m_species_origin=c1, m_sp_area_time=c2)
        #----

        # using min_age and max_age ensures that the time bins always span the same amount of time
        sim = fossil_sim.run_simulation(sp_x, min_age=0, max_age=max_clade_age, return_sqs_data=False)
        sim_features = dd.extract_sim_features(sim)
        sim_y = sim['global_true_trajectory']
        batch_features.append(sim_features)
        batch_labels.append(sim_y)
        
        keys = ['time_specific_rate', 'species_specific_rate', 'area_specific_rate',
                'a_var', 'n_bins', 'area_size', 'n_areas', 'n_species', 'n_sampled_species', 'tot_br_length',
                'n_occurrences', 'slope_pr', 'pr_at_origination', 'time_bins_duration', 'eta', 'p_gap',
                'area_size_concentration_prm', 'link_area_size_carrying_capacity',
                'slope_log_sampling', 'intercept_initial_sampling', 'sd_through_time', 'additional_info']

        s = {key: sim[key] for key in keys}
        sim_settings.append(s)

    res = {'features': batch_features, 
           'labels': batch_labels,
           'settings': sim_settings}
    return res


if __name__ == "__main__":
    ### SIMULATE MULTIPLE DATASETS ###
    # parallel simulations
    if n_training_simulations:

        list_args = list(np.arange(n_CPUS))
        print("\nSimulating training data...")
        pool = multiprocessing.Pool()
        res = pool.map(run_sim, list_args)
        pool.close()

        features = []
        labels = []

        for i in range(n_CPUS):
            features = features + res[i]['features']
            labels = labels + res[i]['labels']

        Xt = np.array(features)
        Yt = np.array(labels)
        print(Xt.shape, Yt.shape)
        # save simulations
        np.save(os.path.join(output_path, outname + "_features" + now + ".npy"), Xt)
        np.save(os.path.join(output_path, outname + "_labels" + now + ".npy"), Yt)

        print("Training features saved as: \n", os.path.join(output_path, outname + "_features" + now + ".npy"))
        print("Training labels saved as: \n", os.path.join(output_path, outname + "_labels" + now + ".npy"))

    ### TEST DATASETS ###
    if n_test_simulations:
        res = run_test_sim(0)

        features = res['features']
        labels = res['labels']

        Xt = np.array(features)
        Yt = np.array(labels)
        print(Xt.shape, Yt.shape)
        # save simulations
        np.save(os.path.join(output_path, outname + "test_features" + now + ".npy"), Xt)
        np.save(os.path.join(output_path, outname + "test_labels" + now + ".npy"), Yt)

        print("Test features saved as: \n", os.path.join(output_path, outname + "test_features" + now + ".npy"))
        print("Test labels saved as: \n", os.path.join(output_path, outname + "test_labels" + now + ".npy"))
        sim_settings = res['settings']
        dd.save_pkl(sim_settings, os.path.join(output_path, "test_sim_settings" + now + ".pkl"))
        
    print("\ndone.\n")
