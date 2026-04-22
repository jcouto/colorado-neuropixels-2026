from pathlib import Path
import os
import numpy as np
from os.path import join as pjoin
from tqdm import tqdm
from glob import glob
import pandas as pd

###############################################################################
################# HANDLE RAW BINARY FILES #####################################
###############################################################################

def map_binary(fname,nchannels,dtype=np.int16,
               offset = 0,
               mode = 'r',nsamples = None,transpose = False):
    ''' 
    dat = map_binary(fname,nchannels,dtype=np.int16,mode = 'r',nsamples = None)
    
Memory maps a binary file to numpy array.
    Inputs: 
        fname           : path to the file
        nchannels       : number of channels
        dtype (int16)   : datatype
        mode ('r')      : mode to open file ('w' - overwrites/creates; 'a' - allows overwriting samples)
        nsamples (None) : number of samples (if None - gets nsamples from the filesize, nchannels and dtype)
    Outputs:
        data            : numpy.memmap object (nchannels x nsamples array)
See also: map_spikeglx, numpy.memmap

    Usage:
Plot a chunk of data:
    dat = map_binary(filename, nchannels = 385)
    chunk = dat[:-150,3000:6000]
    
    import pylab as plt
    offset = 40
    fig = plt.figure(figsize=(10,13)); fig.add_axes([0,0,1,1])
    plt.plot(chunk.T - np.nanmedian(chunk,axis = 1) + offset * np.arange(chunk.shape[0]), lw = 0.5 ,color = 'k');
    plt.axis('tight');plt.axis('off');

    Joao Couto 2019
    '''
    dt = np.dtype(dtype)
    if not os.path.exists(fname):
        if not mode == 'w':
            raise(ValueError('File '+ fname +' does not exist?'))
        else:
            print('Does not exist, will create [{0}].'.format(fname))
            if not os.path.isdir(os.path.dirname(fname)):
                os.makedirs(os.path.dirname(fname))
    if nsamples is None:
        if not os.path.exists(fname):
            raise(ValueError('Need nsamples to create new file.'))
        # Get the number of samples from the file size
        nsamples = os.path.getsize(fname)/(nchannels*dt.itemsize)
    ret = np.memmap(fname,
                    mode=mode,
                    dtype=dt,
                    shape = (int(nsamples),int(nchannels)))
    if transpose:
        ret = ret.transpose([1,0])
    return ret

###############################################################################
################# HANDLE SYNCHRONIZATION AND THE SYNC BYTES ###################
###############################################################################

def unpackbits(x,num_bits = 16):
    '''
    unpacks numbers in bits.

    Joao Couto - April 2019
    '''
    xshape = list(x.shape)
    x = x.reshape([-1,1])
    to_and = 2**np.arange(num_bits).reshape([1,num_bits])
    return (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])


def unpack_npix_sync(syncdat,srate=1,output_binary = False):
    '''Unpacks neuropixels phase external input data
events = unpack_npix3a_sync(trigger_data_channel)    
    Inputs:
        syncdat               : trigger data channel to unpack (pass the last channel of the memory mapped file)
        srate (1)             : sampling rate of the data; to convert to time - meta['imSampRate']
        output_binary (False) : outputs the unpacked signal
    Outputs
        events        : dictionary of events. the keys are the channel number, the items the sample times of the events.

    Joao Couto - April 2019

    Usage:
Load and get trigger times in seconds:
    dat,meta = load_spikeglx('test3a.imec.lf.bin')
    srate = meta['imSampRate']
    onsets,offsets = unpack_npix_sync(dat[:,-1],srate);
Plot events:
    plt.figure(figsize = [10,4])
    for ichan,times in onsets.items():
        plt.vlines(times,ichan,ichan+.8,linewidth = 0.5)
    plt.ylabel('Sync channel number'); plt.xlabel('time (s)')
    '''
    dd = unpackbits(syncdat.flatten(),16)
    mult = 1
    if output_binary:
        return dd
    sync_idx_onset = np.where(mult*np.diff(dd,axis = 0)>0)
    sync_idx_offset = np.where(mult*np.diff(dd,axis = 0)<0)
    onsets = {}
    offsets = {}
    for ichan in np.unique(sync_idx_onset[1]):
        onsets[ichan] = sync_idx_onset[0][
            sync_idx_onset[1] == ichan]/srate
    for ichan in np.unique(sync_idx_offset[1]):
        offsets[ichan] = sync_idx_offset[0][
            sync_idx_offset[1] == ichan]/srate
    return onsets,offsets

###############################################################################
############################ LOAD PHY DATA ####################################
###############################################################################

def waveforms_position(
    waveforms,
    channel_positions,
    active_electrode_threshold=3,
    max_waveform_extent=100,
    ):
    '''Calculates the position of a unit in a set of channels using the center of mass.
Considers only electrodes that have over active_electrode_threshold (3) mad and 
are within max_waveform_extent (100um) from the principal (max) electrode.

Using the max_waveform_extent is useful when there is noise in the recording. 

centerofmass,peak_channels = waveforms_position(waveforms,channel_positions)

Inputs
------
waveforms : array [ncluster x nsamples x nchannels]
    average waveform for a cluster 
channel_positions : array [nchannels x 2]
    x and y coordinates of each channel

Return
-------
centerofmass: array [nchannels x 2]
    center of mass of the waveforms 
peak_channels array [nchannels x 1]
    peak channel of the waveform (the argmax of the absolute amplitude)

Joao Couto - spks 2023
    '''
    nclusters, nsamples, nchannels = waveforms.shape
    N = int(nsamples/4)
    peak_to_peak = waveforms.max(axis=1) - waveforms.min(axis=1)
    # get the threshold from the median_abs_deviation
    channel_mad = np.median(peak_to_peak/0.6745, axis = 1)
    active_electrodes = []
    center_of_mass = []
    peak_channels = []
    for i,w in enumerate(peak_to_peak):
        peak_channels.append(np.argmax(w)) # the peak channel is the index of the channel that has the largest deflection
        idx = np.where(
            (w>(channel_mad[i]*active_electrode_threshold)) & 
            (np.linalg.norm(channel_positions - channel_positions[np.argmax(w)],axis = 1) < max_waveform_extent))[0]
        active_electrodes.append(idx)
        if not len(idx): # then there are no active channels..
            center_of_mass.append([np.nan]*2)
            continue
        # compute the center of mass (X,Y) of the waveforms using only significant peaks
        com = [w[idx]*pos for pos in channel_positions[idx].T]
        center_of_mass.append(np.sum(com,axis = 1)/np.sum(w[idx]))
    return np.array(center_of_mass), np.array(peak_channels), active_electrodes 
    return np.array(center_of_mass), np.array(peak_channels), active_electrodes

def compute_spike_amplitudes(templates,whitening_matrix,spike_templates,spike_template_amplitudes, channel_positions):
    '''
    Compute the amplitude of each spike from the template fitting
    '''

    templates_raw = np.dot(templates,whitening_matrix)
    # compute the peak to peak of each template
    templates_peak_to_peak = (templates_raw.max(axis = 1) - templates_raw.min(axis = 1))
    # the amplitude of each template is the max of the peak difference for all channels
    templates_amplitude = templates_peak_to_peak.max(axis=1)
    templates_amplitude = templates_amplitude.copy()
    # Fix for when kilosort returns NaN templates, make them the average of all templates
    templates_amplitude[~np.isfinite(templates_amplitude)] = np.nanmean(templates_amplitude)
    # compute the center of mass (X,Y) of the templates
    template_position,template_channel, electrode_channels = waveforms_position(templates_raw, channel_positions)
    # get the spike positions and amplitudes from the average templates
    spike_amplitudes = np.take(templates_amplitude,spike_templates)*spike_template_amplitudes
    return spike_amplitudes

def estimate_spike_positions_from_features(spike_templates,spike_pc_features,template_pc_features_ind,channel_positions,consider_feature=0):
    '''
    Estimates the spike 2d location based on a feature e.g the PCs.
    
    This is adapted from the cortexlab/spikes repository to estimate spikes based on the PC features.

    Parameters
    ----------
    spike_templates: nspikes x templates used for each spike
    spike_pc_features: nspikes x nfeatures x nchannels
    template_pc_features_ind: indice of the channels for the templates nchannels
    channel_positions: position of each channel
    consider_feature: feature to consider

    Returns
    -------
    spike_locations: nspikes

    Joao Couto - spks 2023
    '''
    # channel index for each feature
    feature_channel_idx = np.take(template_pc_features_ind,spike_templates.flatten().astype(int),axis=0)
    # 2d coordinates for each channel feature
    feature_coords = np.take(channel_positions,feature_channel_idx.flatten().astype(int),axis=0).reshape([*feature_channel_idx.shape,*channel_positions.shape[1:]])
    # ycoords of those channels?
    pc_features = spike_pc_features[:,consider_feature].squeeze()**2 # take the first pc for the features
    spike_locations = (np.sum(feature_coords.transpose((2,0,1))*pc_features,axis=2)/np.sum(pc_features,axis=1)).T
    return spike_locations

def load_phy_folder(folder, analyzer_waveforms = None):
    # we load the results in a dictionary so we don't accidentally confuse results from different sessions
    if analyzer_waveforms is None:
        analyzer_waveforms = folder
    res = dict(
        # spiketimes and other
        spike_times = np.load(folder.rglob('spike_times.npy').__next__()),
        spike_clusters = np.load(folder.rglob('spike_clusters.npy').__next__()),
        spike_templates = np.load(folder.rglob('spike_templates.npy').__next__()),
        pc_features = np.load(folder.rglob('pc_features.npy').__next__()),
        pc_feature_ind = np.load(folder.rglob('pc_feature_ind.npy').__next__()),
        spike_template_amplitudes = np.load(folder.rglob('amplitudes.npy').__next__()),
        # metrics for each cluster
        metrics = pd.read_csv(folder.rglob('metrics.csv').__next__()),
        # waveforms
        channel_indices = np.load(folder.rglob('channel_map.npy').__next__()),
        channel_positions = np.load(folder.rglob('channel_positions.npy').__next__()),
        mean_waveforms = np.load(analyzer_waveforms.rglob('average.npy').__next__()),
        templates = np.load(folder.rglob('templates.npy').__next__()),
        whitening_mat_inv = np.load(folder.rglob('whitening_mat_inv.npy').__next__())
    )

    # estimate the amplitudes from the template fitting
    res['spike_amplitudes'] = compute_spike_amplitudes(templates = res['templates'],
                                                    whitening_matrix= res['whitening_mat_inv'],
                                                    spike_templates = res['spike_templates'],
                                                    spike_template_amplitudes = res['spike_template_amplitudes'],
                                                    channel_positions = res['channel_positions'])
    # estimate the positions from the template fitting features
    res['spike_positions'] = estimate_spike_positions_from_features(spike_templates=res['spike_templates'],
                                        spike_pc_features = res['pc_features'],
                                        template_pc_features_ind = res['pc_feature_ind'],
                                        channel_positions = res['channel_positions'])
    return res



def plot_drift_raster(spiketimes,spikepositions,spikeamplitudes, n_spikes_to_plot = 200000,cmap = 'Spectral_r',clim=[0,10000],markersize = 0.3):
    import pylab as plt

    plt.figure(figsize = [12,4])

    # randomly subsample n_spikes_to_plot spikes
    subsample = np.random.choice(np.arange(len(spiketimes),dtype=int),n_spikes_to_plot,replace = False)

    # sort by the amplitude so the color is seen
    subsample = subsample[np.argsort(spikeamplitudes[subsample])[::-1]]

    spikes = spiketimes[subsample]
    
    plt.scatter(spikes,
                spikepositions[subsample][:,1],
                markersize,
                spikeamplitudes[subsample],
                clim=clim, 
                cmap = cmap)

    plt.ylim([spikepositions[subsample][:,1].min(),spikepositions[subsample][:,1].max()])
    plt.xlim([spikes.min(),spikes.max()])
