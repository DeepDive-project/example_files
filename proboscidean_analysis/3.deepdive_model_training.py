import os
import sys
import deepdive as dd
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from keras import layers
import pickle as pk
import matplotlib.pyplot as plt
import multiprocessing

#------ Script settings:
feature_file = "elephants_features.npy"
label_file = "elephants_labels.npy"
sim_wd = "./simulations_elephants"
model_wd = "./model_elephants"
parallelize_all_models = False
#------

try:
    os.mkdir(model_wd)
except:
    pass

def get_model_settings():
    lstm_nodes = [# [128,64, 64, 32],
                  #[128, 64, 32],
                  [64, 32],
                  [64]
                  #[16],
                  ]

    dense_nodes = [#[64, 32],
                  [8],
                  []
                   ]

    dropout_frac = [0.05]

    loss_f = ['mse'] 

    list_settings = []
    model_n = 0
    for l in lstm_nodes:
        for d in dense_nodes:
            for f in loss_f:
                for o in dropout_frac:
                    out = 'lstm%s_d%s_o%s_%s' % (len(l), len(d), o, f)
                    d_item = {
                        'model_n': model_n,
                        'lstm_nodes': l,
                        'dense_nodes': d,
                        'loss_f': f,
                        'dropout': o,
                        'model_name': out
                    }
                    list_settings.append(d_item)
                    model_n += 1

    return list_settings


def run_model_training(d):
    Xt = np.load(os.path.join(sim_wd, d['feature_file']))
    Yt = np.load(os.path.join(sim_wd, d['label_file']))
    infile_name = feature_file.split('.npy')[0]
    outname = infile_name + d['model_name']

    feature_rescaler = dd.FeatureRescaler(Xt)
    Xt_r = feature_rescaler.feature_rescale(Xt)
    Yt_r = dd.normalize_labels(Yt, rescaler=1, log=True)
    #Yt_r = dd.normalize_labels(Yt, rescaler=0, log=False)
    model = dd.build_rnn(Xt_r,
                         lstm_nodes=d['lstm_nodes'],
                         dense_nodes=d['dense_nodes'],
                         loss_f=d['loss_f'],
                         dropout_rate=d['dropout'])
    verbose = 0
    if d['model_n'] == 0:
        verbose = 1
    history = dd.fit_rnn(Xt_r, Yt_r, model, verbose=verbose, max_epochs=1000, patience=5,batch_size=100)
    print("\nSaving DeepDive model...\n")
    dd.save_rnn_model(model_wd, history, model, feature_rescaler, filename=outname)
    dd.plot_training_history(history, criterion='val_loss', wd=model_wd, show=False, filename=outname)
    print("done.")


if __name__ == "__main__":
    nametag= 'base file name' 

    list_settings = get_model_settings()
    
    if parallelize_all_models:
        for j in list_settings:
            j['feature_file'] = feature_file
            j['label_file'] = label_file

        # run all jobs in parallel
        pool = multiprocessing.Pool(len(list_settings))
        pool.map(run_model_training, list_settings)
        pool.close()
    else:
        for j in list_settings:
            j['feature_file'] = feature_file
            j['label_file'] = label_file
            run_model_training(j)
        

