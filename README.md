# fNIRS App: GLM Pipeline

[![build](https://github.com/rob-luke/fnirs-apps-glm-pipeline/actions/workflows/ghregistry.yml/badge.svg)](https://github.com/rob-luke/fnirs-apps-glm-pipeline/actions/workflows/ghregistry.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4998538.svg)](https://doi.org/10.5281/zenodo.4998538)


This [*fNIRS App*](http://fnirs-apps.org) will runs a GLM pipeline on your data and export subject and group level results.

**Feedback is welcome!!** Please let me know your experience with this app by raising an issue. 

## Usage

To run the app you must have [docker installed](https://docs.docker.com/get-docker/). See here for details about [installing fNIRS Apps](http://fnirs-apps.org/details/). You do NOT need to have MATLAB or python installed, and you do not need any scripts.

To run the app you must inform it where the `bids_dataset` to be formatted resides.
This is done by passing the app the location of the dataset using the `-v` command.
To run this app use the command:

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app
```

By default the app will process all subject and tasks.
You can also specify additional parameters by passing arguments to the app. A complete list of arguments is provided below.
A more complete example that only runs on participant 3 and also specifies that short channel regression should not be used, can be run as:

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app \
    --subject-label 03 \
    --short-regression False
```

## Pipeline details

1. The pipeline converts the data to optical density.
2. Then applies the Beer Lambert Law conversion.
4. The data is resampled
    - to 0.6 Hz. (TODO: make optional variable)
6. A GLM is applied
    - The model used is a glover HRF convolved with boxcar of stimulus duration
    - (TODO: expose more parameters here)
    - If `--short-regression True` is specified the short channels will be added as regressors to the design matrix for the GLM computation.
    - Drift components will be added to the design matrix using a cosine model including frequencies up to 0.01 Hz (TODO: make user specified parameter).
7. The results of the GLM will be exported per subject as a tidy csv file per run.
8. A mixed effects model is then run on the individual data to produce a summary result
    - First, for each subject the channels are combined to a single region of interest
    - The model is then applied that includes condition and chromaphore as a factor with id as a random variable
    - A summary is exported as a text file


## Arguments

|                   | Required | Default | Note                                                |
|-------------------|----------|---------|-----------------------------------------------------|
| short-regression  | optional | True    | Include short channels as regressor.                |
| export-drifts     | optional | False   | Export the drift coefficients in csv.               |
| export-shorts     | optional | False   | Export the short channel coefficients in csv.       |
| subject-label     | optional | []      | Sujects to process. Default is to process all.      |
| task-label        | optional | []      | Tasks to process. Default is to process all.        |


For example

```bash
docker run -v /path/to/data/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app --short-regression=True --export-shorts=True --subject-label 02 04
```

## Updating

To update to the latest version run.

```bash
docker pull ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app
```

Or to run a specific version:

```bash
docker run -v /path/:/bids_dataset ghcr.io/rob-luke/fnirs-apps-glm-pipeline/app:v1.4.2
```


Acknowledgements
----------------

This app is directly based on BIDS Apps and BIDS Execution. Please cite those projects when using this app.

BIDS Apps: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005209

BIDS Execution: https://github.com/bids-standard/bids-specification/issues/313

This app uses MNE-Python, MNE-BIDS, and MNE-NIRS under the hood. Please cite those package accordingly.

MNE-Python: https://mne.tools/dev/overview/cite.html

MNE-BIDS: https://github.com/mne-tools/mne-bids#citing

MNE-NIRS: https://github.com/mne-tools/mne-nirs#acknowledgements
