#!/usr/bin/env python3
import numpy as np
import sys
import argparse
import pandas as pd
from mne_bids import BIDSPath, read_raw_bids
from glob import glob
import os.path as op
from mne.preprocessing.nirs import optical_density, beer_lambert_law
from mne_nirs.statistics import run_GLM
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.channels import get_short_channels, get_long_channels
from mne_nirs.utils._io import glm_to_tidy
from mne.utils import warn
import statsmodels.formula.api as smf
from mne_nirs.statistics import statsmodels_to_results, glm_region_of_interest
import os
import subprocess

__version__ = "v0.1.0"

def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        line = str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() != None:
            break
    if process.returncode != 0:
        raise Exception("Non zero return code: %d" % process.returncode)

parser = argparse.ArgumentParser(description='Quality Reports')
parser.add_argument('--input-datasets', default="/bids_dataset", type=str,
                    help='The directory with the input dataset formatted according to the BIDS standard.')
parser.add_argument('--output-location', default="/bids_dataset/derivatives/fnirs-apps-glm-pipeline",
                    type=str, help='The directory where the output files should be stored.')
parser.add_argument('--subject-label',
                    help='The label(s) of the participant(s) that should be '
                    'analyzed. The label corresponds to '
                    'sub-<subject-label> from the BIDS spec (so it does '
                    'not include "sub-"). If this parameter is not provided '
                    'all subjects should be analyzed. Multiple participants '
                    'can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('--task-label',
                    help='The label(s) of the tasks(s) that should be '
                    'analyzed. If this parameter is not provided '
                    'all tasks should be analyzed. Multiple tasks '
                    'can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('--short-regression', type=bool, default=True,
                    help='Include short channels as regressors.')
parser.add_argument('--export-drifts', type=bool, default=False,
                    help='Export the drift coefficients in csv.')
parser.add_argument('--export-shorts', type=bool, default=False,
                    help='Export the short channel coefficients in csv.')
parser.add_argument('-v', '--version', action='version',
                    version='BIDS-App Scalp Coupling Index version '
                    f'{__version__}')
args = parser.parse_args()


########################################
# Extract parameters
########################################


ids = []
# only for a subset of subjects
if args.subject_label:
    ids = args.subject_label
    print(f"Processing requested participants: {ids}")
# for all subjects
else:
    subject_dirs = glob(op.join(args.input_datasets, "sub-*"))
    ids = [subject_dir.split("-")[-1] for
           subject_dir in subject_dirs]
    print(f"No participants specified, processing: {ids}")


tasks = []
if args.task_label:
    tasks = args.task_label
    print(f"Processing requested tasks: {tasks}")
else:
    all_snirfs = glob(f"{args.input_datasets}/**/*_nirs.snirf", recursive=True)
    for a in all_snirfs:
        s = a.split("_task-")[1]
        s = s.split("_nirs.snirf")[0]
        tasks.append(s)
    tasks = np.unique(tasks)
    print(f"No tasks specified, processing {tasks}")


########################################
# Report Sections
########################################

def individual_analysis(bids_path, ID, srate=0.6, short=True):

    raw_intensity = read_raw_bids(bids_path=bids_path, verbose=True)

    # Convert signal to haemoglobin and resample
    raw_od = optical_density(raw_intensity)
    raw_haemo = beer_lambert_law(raw_od)
    raw_haemo.resample(srate, verbose=True)

    # Cut out just the short channels for creating a GLM repressor
    sht_chans = get_short_channels(raw_haemo)
    raw_haemo = get_long_channels(raw_haemo)

    if ~np.all(raw_haemo.annotations.duration ==
               raw_haemo.annotations.duration[0]):
        raise ValueError("Support is only available for experiments where "
                         "all durations are the same. "
                         "See https://github.com/rob-luke/"
                         "fnirs-apps-glm-pipeline/issues/1")
    stim_dur = raw_haemo.annotations.duration[0]
    print(f"Fitting HRF with duration {stim_dur} seconds.")

    # Create a design matrix
    design_matrix = make_first_level_design_matrix(raw_haemo,
                                                   stim_dur=stim_dur)

    # Append short channels mean to design matrix
    if short:
        for idx in range(len(sht_chans.ch_names)):
            design_matrix[f"short{idx}"] = sht_chans.copy().get_data()[idx]

    # Run GLM
    glm_est = run_GLM(raw_haemo, design_matrix)

    # Extract channel metrics
    cha = glm_to_tidy(raw_haemo, glm_est, design_matrix)
    cha["ID"] = ID  # Add the participant ID to the dataframe

    # Convert to uM for nicer plotting below.
    cha["theta"] = [t * 1.e6 for t in cha["theta"]]

    # Define channels in each region of interest
    # Then generate the correct indices for each pair
    groups = dict(AllChannels=range(len(raw_haemo.ch_names)))
    # Compute region of interest results from channel data
    roi = pd.DataFrame()
    for idx, col in enumerate(design_matrix.columns):
        roi = roi.append(glm_region_of_interest(glm_est, groups, idx, col))
    roi["ID"] = ID  # Add the participant ID to the dataframe
    roi["theta"] = [t * 1.e6 for t in roi["theta"]]

    return raw_haemo, cha, roi


########################################
# Main script
########################################

print(" ")
df_cha = pd.DataFrame()
df_roi = pd.DataFrame()
for id in ids:
    for task in tasks:
        b_path = BIDSPath(subject=id, task=task,
                          root=f"{args.input_datasets}",
                          datatype="nirs", suffix="nirs",
                          extension=".snirf")
        try:
            raw, cha, roi = individual_analysis(b_path, id,
                                           short=args.short_regression)
            p_out = b_path.update(root=f"{args.output_location}",
                                       extension='.csv',
                                       suffix='glm',
                                       check=False)
            p_out.fpath.parent.mkdir(exist_ok=True, parents=True)

            if args.export_drifts is False:
                cha = cha[~cha.Condition.str.contains("drift")]
                cha = cha[~cha.Condition.str.contains("constant")]
            if args.export_shorts is False:
                cha = cha[~cha.Condition.str.contains("short")]
            print(f"Writing subject results to: {p_out.fpath}")
            cha.to_csv(p_out.fpath, index=False)
            df_cha = df_cha.append(cha)
            df_roi = df_roi.append(roi)
        except FileNotFoundError:
            print(f"Unable to process {b_path.fpath}")
          
          
if len(ids) > 2:
  
  print("Computing group level channel results")
  df_cha = df_cha.query("Chroma in ['hbo']")
  ch_model = smf.mixedlm("theta ~ -1 + ch_name:Chroma:Condition",
                         df_cha, groups=df_cha["ID"]).fit(method='nm')
  ch_model_df = statsmodels_to_results(ch_model)
  group_path = f"{args.output_location}/group.csv"
  ch_model_df.to_csv(group_path)


  print("Computing group level result per conditions as single ROI.")
  df_roi = df_roi[~df_roi.Condition.str.contains("drift")]
  df_roi = df_roi[~df_roi.Condition.str.contains("constant")]
  df_roi = df_roi[~df_roi.Condition.str.contains("short")]
  roi_model = smf.mixedlm("theta ~ -1 + ROI:Condition:Chroma",
                          df_roi, groups=df_roi["ID"]).fit(method='nm')
  sys.stdout = open(f"{args.output_location}/stats.txt", "w")
  print(roi_model.summary())
  sys.stdout.close()
