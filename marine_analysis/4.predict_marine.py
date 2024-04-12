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
import copy

np.random.seed(123)

#------ Script settings:
n_predictions = 5  # number of predictions per input file
replicates = 10  # number of age randomisation replicates used in data_pipeline.R

data_wd = './marine_deepdive_data'
testset_wd = "./simulations"
model_wd = "./model_marine"
output_wd = "./output"

feature_file = "marine_test_features.npy"
label_file = "marine_test_labels.npy"
res_file = "marine_results"
#------


try:
    os.mkdir(output_wd)
except:
    pass

#  Specify settings
level = "Genus_occurrences"
# scaling_options: None, "1-mean", "first-bin"
scaling = "1-mean"

# run predictions across all models
model_list = glob.glob(os.path.join(model_wd, "*rnn_model*"))

# create time bin indices from recent to old
time_bins = np.sort(np.array([259.51, 254.14, 251.902, 251.2, 247.2, 242, 237, 227, 208.5, 201.4, 199.5, 192.9]))

min_age = np.min(time_bins)
# Shifting the time bins to finish at 0 ma, truncates simulations at the end of the study interval
# Remember to shift these back when you plot predictions!
time_bins = time_bins - min_age
delta_t = np.diff(time_bins)

### PLOT BY MODEL
def plot_all_models():
    fig = plt.figure(figsize=(12, 8))
    prediction_color = "b"
    alpha = 0.2
    predictions = []

    for model_i in model_list:
        filename = model_i.split(sep="rnn_model")[1]
        print("\nLoading model:", filename)
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
            plot_time_axis = np.sort(time_bins) + min_age

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

    dd.add_geochrono(0, -0.15, max_ma=-259.51, min_ma=-192.9)
    plt.ylim(bottom=-0.15, top=2.5)
    dd.add_pt_events(height=2.5)
    plt.xlim(-259.51, -192.9)
    plt.ylabel("Relative deviation from mean genus diversity", fontsize=15)
    plt.xlabel("Time (Ma)", fontsize=15)
    fig.show()
    file_name = os.path.join(output_wd, res_file + ".pdf")
    marine_plot = matplotlib.backends.backend_pdf.PdfPages(file_name)
    marine_plot.savefig(fig)
    marine_plot.close()
    print("\nPlot saved as:", file_name)
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

        Xt_r = feature_rescaler.feature_rescale(f)
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
