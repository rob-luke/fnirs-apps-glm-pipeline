# fNIRS App: GLM Pipeline

[![build](https://github.com/rob-luke/fnirs-apps-glm-pipeline/actions/workflows/ghregistry.yml/badge.svg)](https://github.com/rob-luke/fnirs-apps-glm-pipeline/actions/workflows/ghregistry.yml)

http://fnirs-apps.org : Portable fNIRS neuroimaging pipelines that work with BIDS datasets.

This app runs a GLM pipeline on your data.
The pipeline converts the data to optical density and then applies the Beer Lambert Law conversion.
The data is resampled to 0.6 Hz.
A GLM is applied (TODO: expose more parameters here) using the duration of the stimulus convolved with a glover HRF.
If `--short_regression` is specified the short channels will be added as regressors to the design matrix for the GLM computation.
Drift components will be added to the design matrix using a cosine model including frequencies up to 0.01 Hz (TODO: make user specified parameter).

The results of the GLM will be exported per subject as a tidy csv file per run.
A summary is generated of the group level results by running a mixed effects model on the data and exported as a text file.

## Usage

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app
```

By default the app will process all subject and tasks.
You can modify the behaviour of the script using the options below.

## Arguments

|                   | Required | Default | Note                                                |
|-------------------|----------|---------|-----------------------------------------------------|
| short_regression  | optional | True    | Include short channels as regressor.                |
| export_drifts     | optional | False   | Export the drift coefficents in csv.                |
| export_shorts     | optional | False   | Export the short channel coefficents in csv.        |
| participant_label | optional | []      | Participants to process. Default is to process all. |
| task_label        | optional | []      | Tasks to process. Default is to process all.        |


For example

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app --short_regression=True --export_shorts=True
```

## Updating

To update to the latest version run.

```bash
docker pull ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app
```


Acknowledgements
----------------

This package uses MNE-Python, MNE-BIDS, and MNE-NIRS under the hood. Please cite those package accordingly.

MNE-Python: https://mne.tools/dev/overview/cite.html

MNE-BIDS: https://github.com/mne-tools/mne-bids#citing

MNE-NIRS: https://github.com/mne-tools/mne-nirs#acknowledgements
