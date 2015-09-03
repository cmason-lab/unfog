#!/usr/bin/python

"""
Fingerprinting function and accessories
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
def fingerprint_events(events, plot_spectro=False, plot_spectro_name="fp_events.png", plot=False, plot_name="fp_peaks.png", fs=2, nfft=128, noverlap=64, outputdir=""):
    fingerprint = []
    spectrum, freqs, t = mlab.specgram(events, Fs=fs, NFFT=nfft, noverlap=noverlap)  #, NFFT=wsize, Fs=Fs, window=mlab.window_hanning, noverlap=int(wsize * wratio))

    # TODO: There are two different ways to plot the spectrogram here. The first doesn't have the peaks plotted, while the second does. The x-axis needs to be fixed from "Time" to events and the indexing for the second plotting method is not linear.

    # http://stackoverflow.com/questions/24775339/financial-time-series-python-matplotlib-specgram-y-axis-displaying-period-ins
    # calculate the bin limits in time (x dir)
    # note that there are n+1 fence posts
    try:
        dt = t[1] - t[0]
    except IndexError:
        dt = 0
        #print("t length 1")
    t_edge = np.empty(len(t) + 1)
    t_edge[:-1] = t - dt / 2.
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

    if plot_spectro:
        # we'll plot both
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.pcolormesh(t_edge, freq_edge, spectrum)
        fig.savefig(outputdir + "/" + plot_spectro_name)

    # Now that spectrum is calculated, need to extract peaks
    struct = generate_binary_structure(2, 1)
    neighborhood = iterate_structure(struct, 7)
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
    amp_min = 10
    peaks_filtered = [x for x in peaks if x[2] > amp_min] # freq, time, amp
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

    return generate_hashes(local_maxima, fan_value=15)


def comp(seq):
    nuc_table = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'U':'A'}
    comp = [nuc_table[s] for s in list(seq)]
    return "".join(comp)


# Hashing function
def generate_hashes(peaks, fan_value=15):
    """
    Hash list structure:
    sha1_hash[0:20] time_offset
    [(e05b341a9b77a51fd26, 32), ... ]
    """
    MIN_HASH_TIME_DELTA = 0
    MAX_HASH_TIME_DELTA = 200
    FINGERPRINT_REDUCTION = 40
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

