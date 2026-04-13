# Tutorials

This folder contains tutorials for the 2026 Colorado Neuropixels course. The tutorials are designed to guide you through the full workflow of working with Neuropixels electrophysiology data: from understanding the raw binary files produced by SpikeGLX, through spike sorting and quality control, to the analysis of neural responses.

The tutorials are organized into four sections, meant to be followed in the order listed below.

---

## Table of Contents

1. [Before You Begin](#before-you-begin)
2. [Spike Sorting Pipeline](#1-spike-sorting-pipeline)
3. [Clock Synchronization](#2-clock-synchronization)
4. [AL032 Chronic Dataset Tutorials](#3-al032-chronic-dataset-tutorials)
5. [Visual Dataset Course Tutorials](#4-visual-dataset-course-tutorials)
6. [Troubleshooting Common Issues](#troubleshooting-common-issues)

---

## Before You Begin

### Install the Python environment

> [!IMPORTANT]
> Complete installation instructions for software for the course are included in the `sheets` folder of this repository:
> [`../sheets/u24_student_software_installation.pdf`](../sheets/u24_student_software_installation.pdf)
>
> This document describes how to set up the `npix-test` conda environment that includes all required packages (NumPy, SciPy, pandas, Matplotlib, SpikeInterface, Kilosort and others).

### Launch Jupyter from the correct environment or use ``VSCode`` to run these notebooks

``Jupyter notebook`` is a web-based notebook environment for interactive computing. ``VSCode`` can also run notebooks.

After the environment is installed, activate it before starting Jupyter:

```bash
conda activate npix-test
jupyter notebook
```

### Download the tutorial data

Several tutorials require data to be downloaded before the notebook can be run.
The AL032 chronic dataset section begins with a dedicated download notebook (`00_download_data.ipynb`) that automates this step. Run that notebook before attempting the others in that section.

Alternatively, you can download the data and unzip from this link:
	 - [AL32 Output](https://drive.google.com/file/d/1cnxe4GcTI4recrZA3vY52SwfwOry5fU-/view?usp=drive_link)
	 - [Sync example](https://drive.google.com/file/d/1_4va0mQJ3ma31CHsGrB7PU2VaCJkRmJv/view?usp=drive_link)   

---

## 1. Spike Sorting Pipeline

**Folder:** `pipelines/`

These materials demonstrate how to go from raw SpikeGLX recordings to sorted spike data using a combination of established tools.

### `catgt_spikeinterface_pipeline.py` — Full Preprocessing and Sorting Pipeline

This Python script runs a complete spike sorting pipeline in the following stages:

1. **CatGT** — Bandpass filtering (300–10,000 Hz), common average referencing (CAR), and phase correction (tshift)
2. **Motion correction** (optional) — Uses the DREDge algorithm through SpikeInterface
3. **Kilosort4** — Template-matching spike sorting
4. **Quality metrics** — Waveform extraction, amplitude estimates, false positive rate, and a summary metrics table

**Before running, open the script and configure the following:**

| Variable | Location | Description |
|---|---|---|
| `catGTPath` | Line ~33, top of file | Full path to the folder containing the CatGT executable (`runit.bat` on Windows, `runit.sh` on Linux) |
| `raw_data_parent` | Inside `main()` | Path to the folder containing the raw SpikeGLX recording |
| `run_name` | Inside `main()` | Recording name, excluding the gate and trigger suffixes |
| `gate_str` | Inside `main()` | Gate index (e.g., `'0'` for `_g0`) |
| `prb_ind` | Inside `main()` | Probe index (e.g., `0` for `imec0`) |
| `output_parent` | Inside `main()` | Folder where all pipeline output will be written |
| `b_catgt` | Inside `main()` | Set to `False` to skip CatGT if it has already been run |
| `b_sort` | Inside `main()` | Set to `False` to skip sorting (useful to re-run metrics only) |
| `b_useDREDge` | Inside `main()` | Set to `True` to enable motion correction before sorting |

Run the script from the `npix-test` environment:

```bash
conda activate npix-test
python catgt_spikeinterface_pipeline.py
```

### `02_spike_interface_motion.ipynb` — Motion Correction with SpikeInterface

This notebook demonstrates how to apply motion correction to a SpikeGLX recording using `spikeinterface.correct_motion()` with the `"rigid_fast"` preset. It is useful for understanding what motion correction does before incorporating it into a full pipeline.

**Before running, configure:**
- The path to the SpikeGLX recording folder at the top of the notebook

---

## 2. Clock Synchronization

**Folder:** `synchronization/`

### `sglx_clock_calibration_and_sync.ipynb` — Clock Calibration and Multi-Stream Synchronization

When recording with Neuropixels, data are acquired simultaneously from multiple streams — for example, AP and LFP bands from each probe, and auxiliary analog or digital inputs from an NI DAQ or OneBox. Each stream is driven by an independent hardware clock. These clocks run at slightly different rates, and over a long recording (hours), this causes the streams to drift out of alignment by tens to hundreds of milliseconds.

This notebook explains:
- Why clock drift occurs and how large it can be in practice
- How to measure the true sample rate of each stream using a shared synchronization pulse
- How to use SpikeGLX's built-in clock calibration workflow to correct metadata
- How to use CatGT to extract synchronization edges and TPrime to translate event times from one stream's clock to another

**Tools used:** SpikeGLX (for clock calibration), CatGT (edge extraction), TPrime (event time translation)

**Before running, configure:**
- The path to the recording folder at the top of the notebook

---

## 3. AL032 Chronic Dataset Tutorials

**Folder:** `AL_chronic_dataset/`

These notebooks walk through a full analysis of a publicly available chronic Neuropixels dataset. The data come AL032 (recording date 2019-11-21), a mouse implanted with a 4-shank Neuropixels 2.0 probe in primary visual cortex (V1). The recordings were made by Anna Lebedeva and Michael Okun in the Carandini/Harris lab at UCL and are part of Steinmetz et al. 2021. Only shank 0 is used here. The dataset was prepared by Jennifer Colonell (HHMI Janelia).

Run the notebooks in the order listed below.

---

### `00_download_data.ipynb` — Download the Dataset

Downloads the AL032 dataset from Google Drive. This must be completed before running any other notebook in this section.

**Before running, configure:**
- `output_path` — the local folder where the data should be saved (default is the current directory)

Available datasets (passed as the first argument to `download_dataset()`):

| Name | Contents |
|---|---|
| `'chronic_stimulus'` | Stimulus event times and sorting output |
| `'chronic_sorting_output'` | Preprocessed binary and sorting output |
| `'chronic_raw'` | Raw binary files and metadata |

---

### `01_file_formats_plot_raw.ipynb` — Binary File Formats and Raw Data

Introduces the SpikeGLX binary file format. Explains how neural data are stored as a two-dimensional array of 16-bit integers (channels × samples), how the accompanying `.meta` file describes the recording parameters, and how to convert raw integer values to microvolts using the probe gain and voltage range.

**Before running, configure:**
- Path to a `.ap.bin` file from the AL032 dataset

---

### `01_preprocessing_basic_filter.ipynb` — Bandpass Filtering and Common Average Reference

Demonstrates how to apply a software bandpass filter (Butterworth, zero phase) and a common average reference (CAR) to the raw recording. The CAR subtracts the median signal across all channels at each time point, removing artifacts that are coherent across the probe.

**Before running, configure:**
- Path to a `.ap.bin` file from the AL032 dataset

---

### `01_preprocessing_tshift.ipynb` — Phase Correction

Neuropixels probes sample all channels through a time-division multiplexer. This means that, within each sample period, individual channels are measured at slightly different points in time. For high-frequency signals, this introduces a measurable phase offset between channels. This notebook demonstrates how to correct these offsets in the frequency domain using a Fourier-domain phase shift.

**Before running, configure:**
- Path to a `.ap.bin` file from the AL032 dataset

---

### `02_load_output.ipynb` — Loading Spike Sorting Output

Shows how to load the output produced by Kilosort4 and inspected with Phy. Covers the key files (spike times, cluster assignments, templates, waveforms, quality metrics), how to filter units by quality thresholds, and how to produce drift raster plots and depth maps.

**Before running, configure:**
- Path to the sorting output folder (the folder containing `spike_times.npy`, `spike_clusters.npy`, etc.)

---

### `03_visual_response.ipynb` — Visual Responses

Analyzes how individual neurons in V1 respond to natural image stimuli. Covers peri-stimulus time histograms (PSTHs), single-trial rasters, response fingerprints across the 112 stimulus images, and hierarchical clustering of units by response profile.

**Before running, configure:**
- Path to the sorting output folder
- Path to the stimulus events file (`events.csv`)

---

### `05_triggered_LFP.ipynb` — Triggered LFP and Current Source Density

Loads the low-frequency band recording (`.lf.bin`, 2.5 kHz), filters it to the LFP band (3–150 Hz), and computes the stimulus-triggered average LFP across depth. Then computes the current source density (CSD) using the second spatial derivative of the LFP across channels, following Sakata and Harris (2009).

**Before running, configure:**
- Path to the `.lf.bin` file from the AL032 dataset
- Path to the stimulus events file (`events.csv`)

---

## 4. Visual Dataset Course Tutorials

**Folder:** `visual_dataset_course/`

These two notebooks cover the same core concepts as the AL chronic section, but use a different dataset recorded at UCLA in 2026 (cortex, Karolina Socha and Aditya Singh). The data were preprocessed with CatGT and sorted with Kilosort4. Stimuli consisted of drifting gratings at 8 orientations and 3 spatial frequencies.

---

### `01_file_formats_plot_raw.ipynb` — Binary File Formats and Raw Data

The same binary file format concepts as in the AL chronic section, applied to the visual dataset paths. A good starting point if you are working exclusively with the visual dataset.

**Before running, configure:**
- Path to the `.ap.bin` file for this dataset

---

### `visual_responses.ipynb` — Visual Responses

Analyzes grating responses across the recorded population. Includes PSTHs, trial rasters, orientation and spatial frequency tuning, and an interactive widget for browsing individual units.

**Before running, configure:**
- Path to the sorting output folder
- Path to the stimulus events file

---
