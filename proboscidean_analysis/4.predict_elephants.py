import os
import sys
import glob
import deepdive as dd
import numpy as np
from datetime import datetime
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends import backend_pdf  # saves pdfs
import multiprocessing
import scipy.ndimage as nd
import copy
today = datetime.now()
now = datetime.now().strftime('%Y%m%d')

np.random.seed(123)

#------ Script settings:
n_predictions = 5  # number of predictions per input file
replicates = 10  # number of age randomisation replicates used in data_pipeline.R

data_wd = './elephant_deepdive_data'
testset_wd = "./simulations_elephants"
model_wd = "./model_elephants"
output_wd = "./output"

feature_file = "elephants_test_features.npy"
label_file = "elephants_test_labels.npy"
res_file = "elephants_results"
#------


try:
    os.mkdir(output_wd)
except:
    pass
level = "Species_occurrences"
# scaling_options: None, "1-mean", "first-bin"
scaling = None

# run predictions across all models
model_list = glob.glob(os.path.join(model_wd, "*rnn_model*"))

# create time bin indices from recent to old
time_bins = np.sort(np.array([66, 65, 64, 63, 61.6, 60, 59.2, 58.13333, 57.06667, 56, 54.975, 53.95, 52.925, 51.9,
                              50.875, 49.85, 48.825, 47.8, 46.85714, 45.91429, 44.97143, 44.02857, 43.08571, 42.14286,
                              41.2, 40.03667, 38.87333, 37.71, 36.7575, 35.805, 34.8525, 33.9, 32.88667, 31.87333,
                              30.86, 29.84667, 28.83333, 27.82, 26.862, 25.904, 24.946, 23.988, 23.03, 22.16667,
                              21.30333, 20.44, 19.3225, 18.205, 17.0875, 15.97, 14.895, 13.82, 12.725, 11.63, 10.534,
                              9.438, 8.342, 7.246, 6.2895, 5.333, 4.4665, 3.6, 2.58, 1.8, 0.774, 0.129, 0.0117, 0.0082,
                              0.0042, 0]))

n_time_bins = len(time_bins) - 1

# make and plot predictions:
fig = plt.figure(figsize=(12, 8))
predictions = []

def plot_all_models():
    fig = plt.figure(figsize=(12, 8))
    prediction_color = "b"
    alpha = 0.2
    predictions = []
    for model_i in model_list:
        filename = model_i.split(sep="rnn_model")[1]
        print("\nModel", filename)
        # load model trained using age uncertainty
        history, model, feature_rescaler = dd.load_rnn_model(model_wd, filename=filename)

        for replicate in range(1, replicates + 1):
            features, info = dd.prep_dd_input(data_wd,
                                              bin_duration_file='t_bins.csv',  # from old to recent, array of shape (t)
                                              locality_file='%s_localities.csv' % replicate,  # array of shape (a, t)
                                              locality_dir='Locality',
                                              taxon_dir=level,
                                              hr_time_bins=time_bins,  # array of shape (t)
                                              rescale_by_n_bins=True,
                                              no_age_u=True,
                                              replicate=replicate,
                                              debug=False)

            # from recent to old
            plot_time_axis = np.sort(time_bins)

            dd.print_update("Running replicate n. %s" % replicate)

            # from recent to old
            pred_div = dd.predict(features, model, feature_rescaler,
                                  n_predictions=n_predictions, dropout=True)

            pred = np.mean(np.exp(pred_div) - 1, axis=0)
            if scaling == "1-mean":
                den = np.mean(pred)
            elif scaling == "first-bin":
                den = pred[-1]
            else:
                den = 1
            pred /= den

            plt.step(-plot_time_axis,  # pred,
                     [pred[0]] + list(pred),
                     label="Mean prediction",
                     linewidth=2,
                     c=prediction_color,
                     alpha=0.05)

            predictions.append(pred)

        predictions = np.array(predictions)

        dd.add_geochrono(0, -4.8, max_ma=-66, min_ma=0)
        plt.ylim(bottom=-4.8, top=80)
        plt.xlim(-66, 0)
        plt.ylabel("Species diversity", fontsize=15)
        plt.xlabel("Time (Ma)", fontsize=15)
        fig.show()
        file_name = os.path.join(output_wd, res_file + ".pdf")
        ele_plot = matplotlib.backends.backend_pdf.PdfPages(file_name)
        ele_plot.savefig(fig)
        ele_plot.close()
        print("Plot saved as:", file_name)
        return features



if __name__=="__main__":
    # Run predictions and create plots
    features = plot_all_models()
    
    # Get stats for model training in a pandas dataframe
    print("\n\nCalculating MSA on test set...")
    #  Run model on test set to estimate accuracy
    ffile = os.path.join(testset_wd, feature_file)
    lfile = os.path.join(testset_wd, label_file)
    f = np.load(ffile)
    l = np.load(lfile)

    res = list()
    for i in model_list:
        filename = i.split(sep="rnn_model")[1]
        history, model, feature_rescaler = dd.load_rnn_model(model_wd, filename=filename)
        val_loss = np.min(history["val_loss"])  # check validation loss
        t_loss = history["loss"][np.argmin(history["val_loss"])]  # training loss
        epochs = np.argmin(history["val_loss"])  # number of epochs used to train
        try:
            Xt_r = feature_rescaler.feature_rescale(f)
        except:
            Xt_r = f * feature_rescaler
        Yt_r = dd.normalize_labels(l, rescaler=1, log=True)
        # run predictions
        y = np.array(model(Xt_r))
        print("Running model:", filename)
        mse = np.mean((y[:, :, 0] - Yt_r)**2)

        res.append([filename, epochs, t_loss, val_loss,mse])
    res = pd.DataFrame(res)
    res.columns = ["Model", "training_epochs", "training_MSE", "validation_MSE", "test_MSE"]
    res.to_csv( os.path.join(output_wd, res_file + ".csv" ),
                index=False)
    print("Output saved in:", os.path.join(output_wd, res_file + ".csv"))
    
    print("Plotting features")
    dd.plot_feature_hists(test_features=f, empirical_features=features,
        show=False, wd=output_wd, output_name=res_file + "_features")

#
#
# # Get stats for model training in a pandas dataframe
# res = list()
# for i in model_list:
#     filename = i.split(sep="rnn_model")[1]
#     history, model, feature_rescaler = dd.load_rnn_model(model_wd, filename=filename)
#     val_loss = np.min(history["val_loss"])  # check validation loss
#     t_loss = history["loss"][np.argmin(history["val_loss"])]  # training loss
#     epochs = np.argmin(history["val_loss"])  # number of epochs used to train
#     res.append([filename, val_loss, t_loss, epochs])
# res = pd.DataFrame(res)
#
#
# #  Run model on test set to estimate accuracy
# ffile = os.path.join(testset_wd, "elephants_features.npy")
# lfile = os.path.join(testset_wd, "elephants_labels.npy")
# f = np.load(ffile)
# l = np.load(lfile)
#
# for model_i in model_list:
#     filename = model_i.split(sep="rnn_model")[1]
#     history, model, feature_rescaler = dd.load_rnn_model(model_wd, filename)
#     Xt_r = feature_rescaler(f)
#     Yt_r = dd.normalize_labels(l, rescaler=1, log=True)
#     # run predictions
#     y = np.array(model(Xt_r))
#     print(filename)
#     print("MSE: ", np.mean((y[:, :, 0] - Yt_r)**2), "\n")
#
# dd.plot_feature_hists_ele(test_features=f, empirical_features=features, show=True)
