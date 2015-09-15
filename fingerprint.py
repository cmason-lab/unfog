#!/usr/bin/python

"""
Fingerprinting function and accessories, code adapted from https://github.com/worldveil/dejavu/blob/master/dejavu/fingerprint.py
"""

from __future__ import division
import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import (generate_binary_structure, iterate_structure, binary_erosion)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from operator import itemgetter
import hashlib


# Define fingerprinting function
def fingerprint_events(events, plot_spectro=False, plot_spectro_name="fp_events.png", plot=False, plot_name="fp_peaks.png", fs=2, nfft=128, noverlap=64, min_peak_amplitude=10, peak_fan=15, max_hash_time=200, outputdir=""):
    fingerprint = []
    spectrum, freqs, event_t = mlab.specgram(events, Fs=fs, NFFT=nfft, noverlap=noverlap, detrend=mlab.detrend_none)  #, NFFT=wsize, Fs=Fs, window=mlab.window_hanning, noverlap=int(wsize * wratio))

    # TODO: There are two different ways to plot the spectrogram here. The first doesn't have the peaks plotted, while the second does. The x-axis needs to be fixed from "Time" to events and the indexing for the second plotting method is not linear.

    # http://stackoverflow.com/questions/24775339/financial-time-series-python-matplotlib-specgram-y-axis-displaying-period-ins
    # calculate the bin limits in time (x dir)
    # note that there are n+1 fence posts
    try:
        d_event_t = event_t[1] - event_t[0]
    except IndexError:
        d_event_t = 0
        #print("t length 1")
    t_edge = np.empty(len(event_t) + 1)
    t_edge[:-1] = event_t - d_event_t / 2.
    # however, due to the way the spectrogram is calculates, the first and last bins 
    # a bit different:
    t_edge[0] = 0
    t_edge[-1] = t_edge[0] + len(events) / fs

    # calculate the frequency bin limits:
    df = freqs[1] - freqs[0]
    freq_edge = np.empty(len(freqs) + 1)
    freq_edge[:-1] = freqs - df / 2.
    freq_edge[-1] = freq_edge[-2] + df

    # calculate the period bin limits, omit the zero frequency bin
    p_edge = 1. / freq_edge[1:]

    #fft_fig = plt.figure()
    #spectrum = spectrum[(freqs > 2)]
    #freqs = freqs[(freqs > 2)]
    """extent = 0, max(event_t), freqs[0], freqs[-1]
    ax = plt.subplot()
    ax.imshow(spectrum, extent=extent)
    ax.set_xlabel('Event')
    ax.set_ylabel('Frequency')
    ax.set_title('Spectrogram')
    #plt.gca().invert_yaxis()
    
    plt.savefig(outputdir + "/freq_by_event_t_" + str(fs) + "_" + str(nfft) + "_" + str(noverlap))
    """
    #if plot_spectro:
        # we'll plot both
    #    fig = plt.figure()
    #    ax1 = fig.add_subplot(111)
    #    ax1.pcolormesh(t_edge, freq_edge, spectrum)
    #    fig.savefig(outputdir + "/" + plot_spectro_name)

    # Now that spectrum is calculated, need to extract peaks
    struct = generate_binary_structure(2, 1)
    neighborhood = iterate_structure(struct, 5) #TEMP 7
    arr2D = spectrum
    local_max = maximum_filter(arr2D, footprint=neighborhood) == arr2D
    background = (arr2D == 0)
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
    detected_peaks = local_max - eroded_background

    # extract peaks
    amps = arr2D[detected_peaks]
    j, i = np.where(detected_peaks)
    # filter peaks
    amps = amps.flatten()
    peaks = zip(i, j, amps)
    peaks_filtered = [x for x in peaks if x[2] > min_peak_amplitude] # freq, time, amp
    # get indices for frequency and time
    frequency_idx = [x[1] for x in peaks_filtered]
    time_idx = [x[0] for x in peaks_filtered]

    if plot:
        # scatter of the peaks
        fig, ax = plt.subplots()
        ax.imshow(arr2D)
        ax.scatter(time_idx, frequency_idx)
        ax.set_xlabel('Time')
        ax.set_ylabel('Frequency')
        ax.set_title("Spectrogram")
        plt.gca().invert_yaxis()
        fig.savefig(outputdir + "/" + plot_name)

    local_maxima = zip(frequency_idx, time_idx)

    return generate_hashes(local_maxima, peak_fan, max_hash_time) #15 originally, 40 in combination with max hash time 400 was working well for enterobacter


def comp(seq):
    nuc_table = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'U':'A'}
    comp = [nuc_table[s] for s in list(seq)]
    return "".join(comp)


# Hashing function
def generate_hashes(peaks, fan_value, max_hash_time): #degree of pairing between fingerprints and neighbours 
    """
    Hash list structure:
    sha1_hash[0:20] time_offset
    [(e05b341a9b77a51fd26, 32), ... ]
    """
    MIN_HASH_TIME_DELTA = 0
    MAX_HASH_TIME_DELTA = max_hash_time #200 thresholds on how close fingerprints can be in sequence of events in order to be paired
    FINGERPRINT_REDUCTION = 40 #number of bits to discard from front of SHA1 hash (20 in dejavu) - this broke everything
    peaks.sort(key=itemgetter(1))
    for i in range(len(peaks)):
        for j in range(1, fan_value):
            if (i + j) < len(peaks):
                freq1 = peaks[i][0]
                freq2 = peaks[i + j][0]
                t1 = peaks[i][1]
                t2 = peaks[i + j][1]
                t_delta = t2 - t1
                if t_delta >= MIN_HASH_TIME_DELTA and t_delta <= MAX_HASH_TIME_DELTA:
                    h = hashlib.sha1("%s|%s|%s" % (str(freq1), str(freq2), str(t_delta)))
                    yield (h.hexdigest()[0:FINGERPRINT_REDUCTION], t1)

