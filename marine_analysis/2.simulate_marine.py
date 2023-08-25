import os
import sys
import deepdive as dd
import numpy as np
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing

#------ Simulation settings:
n_CPUS = 3  # number of CPUs used for simulations
n_training_simulations = 5  # simulations per CPU (the total should be ~10,000)
n_test_simulations = 10  # simulations per CPU (the total should be e.g. 100 or 1000)
training_seed = 1234
test_seed = 4321
outname = "marine_sim"
#------

try:
    os.mkdir('simulations/')
except:
    pass
try:
    os.mkdir('simulations/sqs_')
except:
    pass
output_path = "simulations/"

use_bins = True
use_area_constraints = False
n_areas = 6
# Shifting the time bins to finish at 0 ma, truncates simulations at the end of the study interval
# Remember to shift these back when you plot predictions!

time_bins = np.sort(np.array([259.51, 254.14, 251.902, 251.2, 247.2, 242, 237, 227, 208.5, 201.4, 199.5, 192.9]))
# boundaries of low resolution time bins
min_age_truncation = np.min(time_bins)
time_bins = time_bins - min_age_truncation
n_time_bins = len(time_bins) - 1
min_age = np.min(time_bins)
root_age = np.max(time_bins)


def plot_feat(features, indx):
    plt.step(-(time_bins + min_age_truncation), [features[0, indx]] + list(features[:, indx]))
    plt.gca().set_title("feature %s" % indx, fontweight="bold", fontsize=24)
    plt.show()


# create simulator object
def create_sim_obj(rseed):
    bd_sim = dd.bd_simulator(s_species=[100, 2000],  # number of starting species
                             rangeSP=[0, 30000],
                             # minEX_SP=0,  # minimum number of extinct lineages allowed
                             minEXTANT_SP=100,  # min number of living species
                             root_r=root_age,  # root age
                             rangeL=[0.05, 0.25],  # range of birth rates - range based on
                             rangeM=[0.05, 0.25],  # range of death rates
                             log_uniform_rates=False,
                             p_mass_extinction=[0, 2 / root_age],  # 50% of dataset expected to have at least 1 ME
                             p_equilibrium=0.5,
                             p_constant_bd=0,
                             p_mass_speciation=[0, 0.001],
                             # probability of mass extinction per my, 2 known mass extinctions in the 70 Ma range
                             poiL=4,  # expected number of birth rate shifts
                             poiM=4,  # expected number of death rate shifts
                             seed=rseed,   # if > 0 fixes the random seed to make simulations reproducible
                             scale=10,
                             vectorize=True)

    fossil_sim = dd.fossil_simulator(n_areas=n_areas,
                                     n_bins=n_time_bins,  # number of time bins
                                     time_bins=time_bins,
                                     eta=[1, 1.75],  # area-sp stochasticity - previously 1, 5
                                     p_gap=[0.01, 0.75],  # probability of 0 preservation in a time bin
                                     dispersal_rate=None,
                                     max_dist=1,
                                     disp_rate_mean=[0, 0.5],
                                     disp_rate_variance=1,
                                     area_mean=100,  # G(a,b) distributed preservation rates across areas, previously 50
                                     area_variance=2,  # previously 5
                                     size_concentration_parameter=[0.1, 3],  # single value or array of length n_areas
                                     link_area_size_carrying_capacity=[1, 10],
                                     # positive, larger numbers = stronger link between area size and carrying capacity
                                     p_origination_a_slope_mean=2,  # mean slope of probability of origination area mean
                                     p_origination_a_slope_sd=0.5,  # std of the slopes
                                     sp_mean=[0.2, 0.5],  # G(a,b) distributed preservation rates across species
                                     sp_variance=2,
                                     slope=[-0.01, 0],  # change in log-sampling rate through time (log-linear)
                                     intercept=[0.1, 0.5],  # initial sampling rate
                                     sd_through_time=[0.001, 0.01],  # st dev in log-sampling rate through time
                                     sd_through_time_skyline=3,
                                     mean_n_epochs_skyline=11,
                                     fraction_skyline_sampling=0.75,
                                     maximum_localities_per_bin=350,
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

        # using min_age and max_age ensures that the time bins always span the same amount of time
        sim = fossil_sim.run_simulation(sp_x, min_age=0, max_age=root_age)
        sim_y = sim['global_true_trajectory']
        sim_features = dd.extract_sim_features(sim)
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
        sim = fossil_sim.run_simulation(sp_x, min_age=0, max_age=root_age)
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
        np.save(os.path.join(output_path, outname + "_features" + ".npy"), Xt)
        np.save(os.path.join(output_path, outname + "_labels" + ".npy"), Yt)

        print("Training features saved as: \n", os.path.join(output_path, outname + "_features" + ".npy"))
        print("Training labels saved as: \n", os.path.join(output_path, outname + "_labels" + ".npy\n"))

    ### TEST DATASETS ###
    if n_test_simulations:
        res = run_test_sim(0)

        features = res['features']
        labels = res['labels']

        Xt = np.array(features)
        Yt = np.array(labels)
        # save simulations
        np.save(os.path.join(output_path, outname + "test_features" + ".npy"), Xt)
        np.save(os.path.join(output_path, outname + "test_labels" + ".npy"), Yt)

        print("Test features saved as: \n", os.path.join(output_path, outname + "_test_features" + ".npy"))
        print("Test labels saved as: \n", os.path.join(output_path, outname + "_test_labels" + ".npy"))
        sim_settings = res['settings']
        dd.save_pkl(sim_settings, os.path.join(output_path, outname + "_test_sim_settings" + ".pkl"))
        
    print("\ndone.\n")
