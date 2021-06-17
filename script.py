#!/usr/bin/env python3
import numpy as np
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
from mne_nirs.statistics import statsmodels_to_results

__version__ = "v0.0.1"


parser = argparse.ArgumentParser(description='Quality Reports')
parser.add_argument('--bids_dir', default="/bids_dataset", type=str,
                    help='The directory with the input dataset '
                    'formatted according to the BIDS standard.')
parser.add_argument('--short_regression', type=bool, default=True,
                    help='Include short channels as regressors.')
parser.add_argument('--export_drifts', type=bool, default=False,
                    help='Export the drift coefficents in csv.')
parser.add_argument('--export_shorts', type=bool, default=False,
                    help='Export the short channel coefficents in csv.')
parser.add_argument('--participant_label',
                    help='The label(s) of the participant(s) that should be '
                    'analyzed. The label corresponds to '
                    'sub-<participant_label> from the BIDS spec (so it does '
                    'not include "sub-"). If this parameter is not provided '
                    'all subjects should be analyzed. Multiple participants '
                    'can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('--task_label',
                    help='The label(s) of the tasks(s) that should be '
                    'analyzed. If this parameter is not provided '
                    'all tasks should be analyzed. Multiple tasks '
                    'can be specified with a space separated list.',
                    nargs="+")
parser.add_argument('-v', '--version', action='version',
                    version='BIDS-App Scalp Coupling Index version '
                    f'{__version__}')
args = parser.parse_args()


########################################
# Extract parameters
########################################


ids = []
# only for a subset of subjects
if args.participant_label:
    ids = args.participant_label
# for all subjects
else:
    subject_dirs = glob(op.join(args.bids_dir, "sub-*"))
    ids = [subject_dir.split("-")[-1] for
           subject_dir in subject_dirs]
    print(f"No participants specified, processing {ids}")


tasks = []
if args.task_label:
    tasks = args.task_label
else:
    all_snirfs = glob("/bids_dataset/**/*_nirs.snirf", recursive=True)
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

    raw_intensity = read_raw_bids(bids_path=bids_path, verbose=False)

    # Convert signal to haemoglobin and resample
    raw_od = optical_density(raw_intensity)
    raw_haemo = beer_lambert_law(raw_od)
    raw_haemo.resample(srate)

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

    return raw_haemo, cha


########################################
# Main script
########################################

print(" ")
df_cha = pd.DataFrame()
for id in ids:
    for task in tasks:
        b_path = BIDSPath(subject=id, task=task,
                          root="/bids_dataset",
                          datatype="nirs", suffix="nirs",
                          extension=".snirf")
        try:
            raw, cha = individual_analysis(b_path, id,
                                           short=args.short_regression)
            p_out = b_path.update(root='/bids_dataset/derivatives/'
                                       'fnirs-apps-glm-pipeline/',
                                       extension='.csv',
                                       suffix='glm',
                                       check=False)
            p_out.fpath.parent.mkdir(exist_ok=True, parents=True)

            if args.export_drifts is False:
                cha = cha[~cha.Condition.str.contains("drift")]
            if args.export_shorts is False:
                cha = cha[~cha.Condition.str.contains("short")]
            cha.to_csv(p_out.fpath, index=False)
            df_cha = df_cha.append(cha)
        except FileNotFoundError:
            print(f"Unable to process {b_path.fpath}")

df_cha = df_cha.query("Chroma in ['hbo']")
ch_model = smf.mixedlm("theta ~ -1 + ch_name:Chroma:Condition",
                       df_cha, groups=df_cha["ID"]).fit(method='nm')
ch_model_df = statsmodels_to_results(ch_model)
group_path = '/bids_dataset/derivatives/fnirs-apps-glm-pipeline/group.csv'
ch_model_df.to_csv(group_path)