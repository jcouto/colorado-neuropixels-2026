# -*- coding: utf-8 -*-
"""
Demonstration spike sorting pipeline using a mix of tools.
CatGT is used for preprocessing (could add event extraction for use with TPrime).
  This is also a useful demo showing how to run a command line utility.
Dredge is called through SpikeInterface for drift correction
KS4 and some metrics calculations are done through SpikeIntterface

To run this pipeline:
1) Follow instructions in 'NP_workshop_software_install.pdf' to create
   the npix-test environment. Select that environment as your Python interpreter.
2) Install CatGT and copy the path to catGTPath
3) Set User params (top of main())

To adapt this code into a "real" pipeline, the calls in main() would be in 
a loop over run names.

@author: colonellj@janelia.hhmi.org
"""


import spikeinterface.full as si
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import subprocess
import os
import time
import sys

# Environment constants for CatGT
catGTPath = r"C:\Users\labadmin\Documents\CatGT-win" # path to CatGT folder, which contains the runit.bat or runit.sh file

# Simple CatGT caller. Many more options availble (see CatGT ReadMe).
# tshift/phase correction happens along with bandpass filtering
# If the probe string is a list, all those probes will be processed
# Specific notes for this implementation: -no_catgt_fld is for compatibility with SpikeInterface.
# Assumes car_mode is either 'gblcar' or 'gbldmx'.
# Any other CatGT options (e.g. edge extractions, etc) can be added by including in the other_cmds_str parameter, 
#    which is added to the command line as-is.

def run_CatGT(parent_folder, dest_folder, run_name, 
              prb_string = '0', 
              gate_string = '0',
              trigger_string = '0,0', 
              stream_string = 'ap',  
              bp_min = 300,
              bp_max = 10000,    
              car_mode = 'gblcar', 
              other_cmds_str = None):

    print('CatGT helper function')

    
    if sys.platform.startswith('win'):
        os_str = 'win'
        # build windows command line
        # catGTexe_fullpath = catGTPath.replace('\\', '/') + "/runit.bat"
        # call catGT directly with params. CatGT.log file will be saved lcoally
        # in current working directory (with the calling script.
        # Output will be saved directly in dest_folder, without creating a catgt_ subfolder.
        # This accomodates si.neo_get_streams, which requires the original SGLX naming.
        catGTexe_fullpath = catGTPath.replace('\\', '/') + "/runit.bat"
    elif sys.platform.startswith('linux'):
        os_str = 'linux'
        catGTexe_fullpath = catGTPath.replace('\\', '/') + "/runit.sh"
    else:
        print('unknown system, cannot run CatGt')
   
    # build filter string
    ap_filter_str = f'-apfilter=butter,12,{bp_min},{bp_max}'
    cmd_parts = list()
    
    cmd_parts.append(catGTexe_fullpath)
    cmd_parts.append('-dir=' + parent_folder)
    cmd_parts.append('-run=' + run_name)
    cmd_parts.append('-g=' + gate_string)
    cmd_parts.append('-t=' + trigger_string)
    cmd_parts.append('-prb=' + prb_string)
    if car_mode is not None:
        cmd_parts.append('-' + car_mode)
    cmd_parts.append('-' + stream_string)
    cmd_parts.append(ap_filter_str)
    if other_cmds_str is not None:
        cmd_parts.append(other_cmds_str)    # other options added by caller
    cmd_parts.append('-prb_fld -out_prb_fld -gfix=0.40,0.10,0.02 -no_catgt_fld')  # useful defaults, check formatting
    cmd_parts.append('-dest=' + dest_folder)
    
    # print('cmd_parts')

    catGT_cmd = ' '        # use space as the separator for the command parts
    catGT_cmd = catGT_cmd.join(cmd_parts[1:len(cmd_parts)]) #
    if os_str=='linux':
        # enclose the params in single quotes, so curly braces will not be interpreted by Linux
        catGT_cmd = f"{cmd_parts[0]} '{catGT_cmd}'"
    else:
        catGT_cmd = f"{cmd_parts[0]} {catGT_cmd}"
    
     
    print('CatGT command line:' + catGT_cmd)
    
    start = time.time()
    subprocess.Popen(catGT_cmd,shell='False').wait()

    execution_time = time.time() - start 
    
    print('total time: ' + str(np.around(execution_time,2)) + ' seconds')

    # check for output
    output_bin = os.path.join(dest_folder, f'{run_name}_g{gate_string}_imec{prb_string}', f'{run_name}_g{gate_string}_tcat.imec{prb_string}.ap.bin')
    if not os.path.exists(output_bin):
        print('CatGT failed.')
        print('Check CatGT.log for details.')
        print(f'CatGT path = catGTPath')
        ok = False
    else:
        ok = True
     
    return ok

def simple_fp_est(phy_folder, ref_per_ms=1.5, sample_rate=30000, min_isi_ms=0):
    # example function for simple measurement derived from spike train
    # uses the model from Hill, et al. J Neurosci 31:8699
    # modefied from the ecephys pipeline
    # input: phy_folder, including 
    # Note that the calculation assumes the spikes are orderd in time.

    # fpRate : rate of contaminating spikes as a fraction of overall rate
    #     A perfect unit has a fpRate = 0
    #     A unit with some contamination has a fpRate < 0.25 
    #     A unit with lots of contamination has a fpRate = 1.0

    spike_times = np.load(os.path.join(phy_folder,'spike_times.npy'))    # in samples
    spike_labels = np.load(os.path.join(phy_folder,'spike_clusters.npy'))
    
    # get the number of units from KS labels rather than labels in spike_clusters
    # because there can be skipped labels.
    unit_label_df = pd.read_csv(os.path.join(phy_folder,'cluster_KSLabel.tsv'))
    [n_unit, n_col] = unit_label_df.shape

    isi_threshold = ref_per_ms/1000
    min_isi = min_isi_ms/1000
    spike_times_sec = spike_times.astype(float)/sample_rate
    rec_time = np.max(spike_times_sec) - np.min(spike_times_sec)

    fp_array = np.zeros((n_unit,))

    for n in range(n_unit):

        spike_train = spike_times_sec[spike_labels==n]

        # spikes with very small inter-spike intervals are assumed to be duplicate detections
        # for data that has already been through duplicate detection (e.g. KS4 output)
        # set min_isi_ms=0
        duplicate_spikes = np.where(np.diff(spike_train) <= min_isi)[0]
        # always remove the 2nd of the two 
        spike_train = np.delete(spike_train, duplicate_spikes + 1)

        isis = np.diff(spike_train)

        num_spikes = len(spike_train)
        if num_spikes > 0:
            num_violations = sum(isis < isi_threshold) 
            violation_time = 2*num_spikes*(isi_threshold - min_isi_ms)
            total_rate = num_spikes/rec_time
            c = num_violations/(violation_time*total_rate)
            if c < 0.25:        # valid solution to quadratic eq. for fpRate:
                fp_array[n] = (1 - np.sqrt(1-4*c))/2
            else:               # no valid solution to eq, call fpRate = 1
                fp_array[n] = 1.0
        
    return fp_array

def amplitude_from_waveforms(wf):
    # example function for simple measurement derived from waveform shape
    pp_all_chan = np.max(wf, axis=1) - np.min(wf, axis=1)  # peak-to-peak for all channels
    amplitude = np.max(pp_all_chan, axis=1)
    peak_chan = np.argmax(pp_all_chan, axis=1)
    return amplitude, peak_chan


def custom_metrics_table(analyzer_folder, phy_folder):
    # create a csv file of useful metrics
    # save to the phy folder, where it can be read in by phy for filtering
    # a number of the quality metrics have already been saved to
    # \analyzer\extensions\quality_metrics
    mdf = pd.read_csv(os.path.join(analyzer_folder,'extensions','quality_metrics','metrics.csv'))
    mdf.rename(columns={'Unnamed: 0': 'cluster_id'}, inplace=True)

    # add columns for amplitude and peak_chan
    # average waveforms are returned from SpikeInterface in uV
    average_wf = np.load(os.path.join(analyzer_folder,'extensions','templates','average.npy'))
    amplitude,peak_chan = amplitude_from_waveforms(average_wf)
    mdf['amplitude'] = amplitude
    mdf['peak_chan'] = peak_chan

    # add columns for unit locations calculated by SpikeInterface
    unit_locations = np.load(os.path.join(analyzer_folder,'extensions','unit_locations','unit_locations.npy'))
    mdf['unit_x'] = unit_locations[:,0]
    mdf['unit_y'] = unit_locations[:,1]
    mdf['unit_z'] = unit_locations[:,2]

    # add column for ecephys/simple Hill false positive estimate
    simple_fp = simple_fp_est(phy_folder)
    mdf['ecephys_fp'] = simple_fp
  
    # write new metrics.csv to phy_folder
    mdf.to_csv(os.path.join(phy_folder,'metrics.csv'),index=False)
    




def main():
    # ------User params---------------------
    # Data details
    raw_data_parent = r'D:\pipeline_test_data'
    run_name = 'SC035_010720_ex'     # excludes gate and trigger
    gate_str = '0'
    prb_ind = 0
    output_parent=r'D:\si_testing\catgt_si_out'

    # CatGT params
    car_mode = "gblcar"
    other_cmds = None  # add commands to extract edges, do other preprocessing here

    b_useDREDge = False 
    # motion 'Preset' selection in SpikeInterface
    # 'dredge' runs the official implementation from Windolf2023
    si_motion_preset = 'dredge'

    # KS4 params. Note that we still run preprocessing in KS4
    # This allows us to use KS4's whitening (which is different from the SI version)
    ks4_params = si.get_default_sorter_params('kilosort4')
    ks4_params['do_CAR'] = False # skip CAR in kilosort
   
    job_kwargs = dict(n_jobs=4, chunk_duration='1s', progress_bar=True) # how to chunk and process data

    # What to run -- set both to True to run from scratch
    b_catgt = True  # set to false to re-run wihout re-running CatGT, for example to test sorting and metrics calculation
    b_sort = True # set to False to skip sorting and just open the Analyzer -- useful for testing metrics calculation alone

    #------End of user params-------------

    rg_name = f'{run_name}_g{gate_str}'
    catgt_out_folder = os.path.join(output_parent, rg_name)
    if not os.path.exists(catgt_out_folder):
        os.mkdir(catgt_out_folder)

    # sort and analyzer folders must not exist (avoids overwriting)
    sort_folder = os.path.join(catgt_out_folder,'kilosort4_output')
    analyzer_folder = os.path.join(catgt_out_folder,'analyzer')

    # run catGT
    if b_catgt:
        ok = run_CatGT(raw_data_parent, catgt_out_folder, run_name, prb_string=str(prb_ind), gate_string=gate_str,
                                     car_mode=car_mode, other_cmds_str=other_cmds)
        if not ok:
            print('CatGT failed, skipping sort.')
            return
   
    # Load the preprocessed spikeglx data 
    stream_names, stream_ids = si.get_neo_streams('spikeglx', catgt_out_folder)
    print(f'Found streams {stream_names} in folder {catgt_out_folder}')

    # find rec that matches the specified prb_ind
    targ_stream = f'imec{prb_ind}.ap'
    try:
        idx_sort = stream_names.index(targ_stream)
    except ValueError:
        idx_sort = -1
    
    
    rec = si.read_spikeglx(catgt_out_folder, stream_name=stream_names[idx_sort], load_sync_channel=False)  
    if b_useDREDge:
        rec_sort = si.correct_motion(recording=rec, preset=si_motion_preset)
        ks4_params['nblocks'] = 0 # skip KS motion correction
    else:
        rec_sort = rec
        ks4_params['nblocks'] = 1 # rigid correction in KS (because this is a sort probe segment)


    if b_sort:    
        # run ks4
        sorting = si.run_sorter('kilosort4', rec_sort, folder=sort_folder,
                            docker_image=False, verbose=True, **ks4_params)
        sorting = si.read_sorter_folder(sort_folder)
        
        analyzer = si.create_sorting_analyzer(sorting, rec_sort, sparse=True, format="memory")
        # compute waveforms and other output
        job_kwargs = dict(n_jobs=4, chunk_duration='1s', progress_bar=True) # how to chunk and process data
        analyzer.compute("random_spikes", method="uniform", max_spikes_per_unit=500)
        start_time = time.perf_counter()
        analyzer.compute("waveforms",  ms_before=0.6,ms_after=1.4, **job_kwargs)
        stop_time = time.perf_counter()
        print(f'time for waveform calculation: {stop_time-start_time}')
        analyzer.compute("templates", operators=["average", "median", "std"])
        analyzer.compute("noise_levels")
        analyzer.compute("correlograms")
        analyzer.compute("unit_locations")
        analyzer.compute("spike_amplitudes", **job_kwargs)       
        
        metric_names=['firing_rate', 'presence_ratio', 'snr', 'isi_violation', 'amplitude_cutoff']
        metrics = si.compute_quality_metrics(analyzer, metric_names=metric_names)
        
        # save
        analyzer_saved = analyzer.save_as(folder=analyzer_folder, format="binary_folder")

    else:
        analyzer = si.load_sorting_analyzer(folder=analyzer_folder)

    # export_report makes good pictures, but it takes a long time!     
    # export_report(sorting_analyzer=analyzer, output_folder=os.path.join(out_dir,"report"))
    custom_metrics_table(analyzer_folder, os.path.join(sort_folder,'sorter_output'))

    
    return
   

if __name__ == "__main__":
        main()  