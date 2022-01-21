#!/usr/bin/env python3
import numpy as np
import sys
import argparse
import pandas as pd
from mne_bids import BIDSPath, read_raw_bids, get_entity_vals
from glob import glob
import os.path as op
from mne.preprocessing.nirs import optical_density, beer_lambert_law
from mne_nirs.statistics import run_glm
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.channels import get_short_channels, get_long_channels
import statsmodels.formula.api as smf
import os
import subprocess
from mne.utils import logger
import mne
from pathlib import Path
from datetime import datetime
import json
import hashlib
from pprint import pprint

__version__ = "v0.3.4"

def fnirsapp_glm(command, env={}):
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


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser(description='GLM')
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
parser.add_argument('--session-label',
                    help='The label(s) of the session(s) that should be '
                    'analyzed. The label corresponds to '
                    'ses-<session-label> from the BIDS spec (so it does '
                    'not include "ses-"). If this parameter is not provided '
                    'all sessions should be analyzed. Multiple sessions '
                    'can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('--task-label',
                    help='The label(s) of the tasks(s) that should be '
                    'analyzed. If this parameter is not provided '
                    'all tasks should be analyzed. Multiple tasks '
                    'can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('--short-regression', type=str2bool, default=True,
                    help='Include short channels as regressors.')
parser.add_argument('--sample-rate', type=float, default=0.6,
                    help='Sample rate to resample data to (Hz).')
parser.add_argument('--export-drifts', type=str2bool, default=False,
                    help='Export the drift coefficients in csv.')
parser.add_argument('--export-shorts', type=str2bool, default=False,
                    help='Export the short channel coefficients in csv.')
parser.add_argument('-v', '--version', action='version',
                    version='BIDS-App Scalp Coupling Index version '
                    f'{__version__}')
args = parser.parse_args()


def create_report(app_name=None, pargs=None):

    exec_rep = dict()
    exec_rep["ExecutionStart"] = datetime.now().isoformat()
    exec_rep["ApplicationName"] = app_name
    exec_rep["ApplicationVersion"] = __version__
    exec_rep["Arguments"] = vars(pargs)

    return exec_rep

exec_files = dict()
exec_rep = create_report(app_name="fNIRS-Apps: GLM Pipeline", pargs=args)
pprint(exec_rep)

mne.set_log_level("INFO")
logger.info("\n")

########################################
# Extract parameters
########################################

logger.info("Extracting subject metadata.")
subs = []
if args.subject_label:
    logger.info("    Subject data provided as input argument.")
    subs = args.subject_label
else:
    logger.info("    Subject data will be extracted from data.")
    subs = get_entity_vals(args.input_datasets, 'subject')
logger.info(f"        Subjects: {subs}")


logger.info("Extracting session metadata.")
sess = []
if args.session_label:
    logger.info("    Session data provided as input argument.")
    sess = args.session_label
else:
    logger.info("    Session data will be extracted from data.")
    sess = get_entity_vals(args.input_datasets, 'session')
if len(sess) == 0:
    sess = [None]
logger.info(f"        Sessions: {sess}")


logger.info("Extracting tasks metadata.")
tasks = []
if args.task_label:
    logger.info("    Task data provided as input argument.")
    tasks = args.task_label
else:
    logger.info("    Session data will be extracted from data.")
    tasks = get_entity_vals(args.input_datasets, 'task')
logger.info(f"        Tasks: {tasks}")


########################################
# Report Sections
########################################

def individual_analysis(bids_path, ID, srate=0.6, short=True):

    raw_intensity = read_raw_bids(bids_path=bids_path, verbose=True)

    # Convert signal to haemoglobin and resample
    raw_od = optical_density(raw_intensity)
    raw_haemo = beer_lambert_law(raw_od)
    logger.info(f"    Resampling to {srate} Hz")
    raw_haemo.resample(srate, verbose=True)

    if short:
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
    logger.info(f"    Fitting HRF with duration {stim_dur} seconds.")

    # Create a design matrix
    design_matrix = make_first_level_design_matrix(raw_haemo,
                                                   stim_dur=stim_dur)

    # Append short channels mean to design matrix
    if short:
        for idx in range(len(sht_chans.ch_names)):
            design_matrix[f"short{idx}"] = sht_chans.copy().get_data()[idx]

    # Run GLM
    glm_est = run_glm(raw_haemo, design_matrix)

    # Extract channel metrics
    channel_df = glm_est.to_dataframe()
    channel_df["ID"] = ID  # Add the participant ID to the dataframe

    # Convert to uM for nicer plotting below.
    channel_df["theta"] = [t * 1.e6 for t in channel_df["theta"]]

    # Define channels in each region of interest
    # Then generate the correct indices for each pair
    groups = dict(AllChannels=range(len(raw_haemo.ch_names)))
    # Compute region of interest results from channel data
    roi = glm_est.to_dataframe_region_of_interest(groups, design_matrix.columns)

    roi["ID"] = ID  # Add the participant ID to the dataframe
    roi["theta"] = [t * 1.e6 for t in roi["theta"]]

    return raw_haemo, channel_df, roi


########################################
# Main script
########################################

print(" ")

df_cha = pd.DataFrame()
df_roi = pd.DataFrame()

for sub in subs:
    for task in tasks:
        for ses in sess:

            logger.info(f"Processing: sub-{sub}/ses-{ses}/task-{task}")
            exec_files[f"sub-{sub}_ses-{ses}_task-{task}"] = dict()

            b_path = BIDSPath(subject=sub, task=task, session=ses,
                              root=f"{args.input_datasets}",
                              datatype="nirs", suffix="nirs",
                              extension=".snirf")
            try:

                exec_files[f"sub-{sub}_ses-{ses}_task-{task}"]["FileName"] = str(b_path.fpath)
                exec_files[f"sub-{sub}_ses-{ses}_task-{task}"]["FileHash"] = hashlib.md5(open(b_path.fpath, 'rb').read()).hexdigest()

                raw, cha, roi = individual_analysis(b_path, sub,
                                                    short=args.short_regression,
                                                    srate=args.sample_rate)
                p_out = b_path.update(root=f"{args.output_location}",
                                      extension='.csv',
                                      suffix='glm', check=False)
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
          
          
if len(subs) > 2:

  print("Computing group level result per conditions as single ROI.")
  df_roi = df_roi[~df_roi.Condition.str.contains("drift")]
  df_roi = df_roi[~df_roi.Condition.str.contains("constant")]
  df_roi = df_roi[~df_roi.Condition.str.contains("short")]
  roi_model = smf.mixedlm("theta ~ -1 + ROI:Condition:Chroma",
                          df_roi, groups=df_roi["ID"]).fit(method='nm')
  print(roi_model.summary())
  sys.stdout = open(f"{args.output_location}/stats.txt", "w")
  print(roi_model.summary())
  sys.stdout.close()

exec_rep["Files"] = exec_files
exec_path = f"{args.input_datasets}/execution"
exec_rep["ExecutionEnd"] = datetime.now().isoformat()

Path(exec_path).mkdir(parents=True, exist_ok=True)
with open(f"{exec_path}/{exec_rep['ExecutionStart'].replace(':', '-')}-glm_pipline.json", "w") as fp:
    json.dump(exec_rep, fp)

