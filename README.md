# fNIRS App: GLM Pipeline

[![build](https://github.com/rob-luke/fnirs-apps-glm-pipeline/actions/workflows/ghregistry.yml/badge.svg)](https://github.com/rob-luke/fnirs-apps-glm-pipeline/actions/workflows/ghregistry.yml)

Portable fNIRS neuroimaging pipelines that work with BIDS datasets. See http://fnirs-apps.org

This app runs a GLM pipeline on your data.
The pipeline converts the data to optical density and then applies the Beer Lambert Law conversion.
The data is resampled to 0.6 Hz.
A GLM is applied (TODO: expose more parameters here) using the duration of the stimulus convolved with a glover HRF.
If `--short_regression` is specified the short channels will be added as regressors to the design matrix for the GLM computation.
Drift components will be added to the design matrix using a cosine model including frequencies up to 0.01 Hz (TODO: make user specified parameter).
The results of the GLM will be export as a tidy csv file per run.

## Usage

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app
```

By default the app will process all subject and tasks.
You can modify the behaviour of the script using the options below.

## Arguments

|                  | Required | Default | Note                                         |
|------------------|----------|---------|----------------------------------------------|
| short_regression | optional | NA      | Include short channels as regressor.         |
| export_drifts    | optional | 0       | Export the drift coefficents in csv.         |
| export_shorts    | optional | end     | Export the short channel coefficents in csv. |


For example

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app --short_regression=True --export_shorts=True
```

Acknowledgements
----------------

This package uses MNE-Python, MNE-BIDS, and MNE-NIRS under the hood. Please cite those package accordingly.

MNE-Python: https://mne.tools/dev/overview/cite.html

MNE-BIDS: https://github.com/mne-tools/mne-bids#citing

MNE-NIRS: https://github.com/mne-tools/mne-nirs#acknowledgements
