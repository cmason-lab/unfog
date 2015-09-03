#!/usr/bin/python

"""
Uncovering Nanopore's Fingerprints of Genomes (UNFOG)

AKA: Shazam for nanopore data.
Based off of https://github.com/worldveil/dejavu.

Angela M Yu, 6/15/2015 - 7/3/2015, Alexa BR McIntyre, 9/1/2015 - 
"""

from __future__ import division

from fingerprint import fingerprint_events
from make_db import (run_db_builder)

#import multiprocessing
#import math
import h5py
import numpy as np
import cPickle as pickle #pickle in Python3
import glob
from collections import defaultdict
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import (generate_binary_structure, iterate_structure, binary_erosion)
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import re
from operator import itemgetter
#import rpy2.robjects.numpy2ri
#rpy2.robjects.numpy2ri.activate()
#from rpy2.robjects.packages import importr
import sys
import os
import getopt

##from rpy2.robjects import r


def increment_dict_counts(in_dict, update_dict):
    ret = in_dict.copy()
    for k in update_dict:
        try:
            ret[k] += update_dict[k]
        except KeyError:
            ret[k] = update_dict[k]
    return ret


def find_matches(samples, geno_db, unique=False):
    """
    Finds matches from reads to genomes
    """
    # TODO: Also report back the genome location (not in 5-mer scale)
    mapper = {}
    matches = []
    for hash, offset in samples:
        mapper[hash] = offset
    for h in mapper.keys():
        for g in geno_db:
            if h in geno_db[g]:
                offset = geno_db[g][h]
                matches.append((g, offset - mapper[h]))

    diff_counter = {}
    largest = 0
    largest_count = 0
    geno_id = []
    for tup in matches:
        gid, diff = tup
        if diff not in diff_counter:
            diff_counter[diff] = {}
        if gid not in diff_counter[diff]:
            diff_counter[diff][gid] = 0
        diff_counter[diff][gid] += 1
        if diff_counter[diff][gid] > largest_count:
            largest = diff
            largest_count = diff_counter[diff][gid]
            geno_id = [gid]
        elif diff_counter[diff][gid] == largest_count:
            geno_id.append(gid)
    if unique and len(geno_id) >1:
        return ([], -1, {})
    return (geno_id, largest_count, diff_counter)


def find_overlap_fp(diff_counter, diff_counter_comp):
    offset_overlap = set(diff_counter.keys()).intersection(set(diff_counter_comp.keys()))
    diff_overlap = {}
    largest_count = -1
    largest_offset = -1
    geno_id = []
    for o in offset_overlap:
        diff_overlap[o] = increment_dict_counts(diff_counter[o], diff_counter_comp[o])
        for g in diff_overlap[o]:
            if diff_overlap[o][g] > largest_count:
                largest_count = diff_overlap[o][g]
                largest_offset = o
                geno_id = [g]
            elif diff_overlap[o][g] == largest_count:
                geno_id.append(g)
    return (geno_id, largest_count, largest_offset, diff_overlap)


def merge_window(indices, window=100):
    res = []
    start_edge = indices[0]
    before = indices[0]
    for i in indices[1:]:
        if i - before >= window:
            res.append((start_edge, before))
            start_edge = i
        if i == indices[-1]:
            res.append((start_edge, indices[-1]))
        before = i
    return res


def find_closest_interval(interval_list, i):
    # assumes non-overlapping intervals
    closest_i = 0
    diff = sys.maxint
    for ii in range(len(interval_list)):
        s, e = interval_list[ii]
        if i >= s and i <= e:
            closest_i = ii
            break
        d = min(abs(s-i), abs(e-i))
        if d < diff:
            diff = d
            closest_i = ii
    return interval_list[closest_i]


def run_UNFOG(input_dir, output_dir, create_db, unique, fs, nfft, noverlap, outlier_thresh, overlap, genomes, plot_spectro, plot, only_genomes, num_processes):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plot_dir = output_dir + "/plot/"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    database_open = False
    read_files = glob.glob(input_dir + "/*fast5")
    print "Number of read files: " + str(len(read_files))
    genome_counts = defaultdict(float)
    index_error_count = 0
    outlier_window_count = 0
    if genomes[-1] == "/":
        genome_files_list = glob.glob(genomes + "*fasta*")
    else:
        genome_files_list = genomes.split(",")
    print "Number of Genomes: " + str(len(genome_files_list))
    for f in read_files:
        fname = re.findall("([^/]+\.fast5)", f)[0] 
        with h5py.File(f, "r") as hf:
            if not database_open: 
                if create_db: 
                # runs only for first read file, assumes that the same model is applied to all the reads analyzed - this is more or less true according to nanoporetech wiki for different pores, not sure about different channels and it seems to depend on the version of the workflow used   
                    db_builder_input = (genome_files_list, num_processes, fs, nfft, noverlap, only_genomes)
                    fingerprint_db, revcomp_fingerprint_db = run_db_builder(hf, db_builder_input, output_dir)
                    # TODO: implement f["Analyses"]["Basecall_2D_000"]["BaseCalled_complement/Model"].attrs.items()?
                    if fingerprint_db and revcomp_fingerprint_db:
                        database_open = True
                    else: 
                        continue #skip to next in list because hmm model files don't exist for this read
                    pickle.dump(fingerprint_db, open(output_dir+"/fingerprint_db.p", "wb"))
                    pickle.dump(revcomp_fingerprint_db, open(output_dir+"/revcomp_fingerprint_db.p", "wb"))
                    if only_genomes:
                        return({},[])
                else:
                    try: 
                        fingerprint_db = pickle.load(open(output_dir+"/fingerprint_db.p", "rb"))
                        revcomp_fingerprint_db = pickle.load(open(output_dir+"/revcomp_fingerprint_db.p", "rb"))
                    except IOError:
                        print("Reference fingerprint database doesn't yet exist! Try --create_db True")

            # Compare the fingerprints from the read to what's in the database
            read_key = (hf["Analyses/EventDetection_000/Reads"].keys())[0]
            f_events = [a[0] for a in hf["Analyses/EventDetection_000/Reads"][read_key]["Events"].value]
            #print len(f_events)
            #print max(f_events)
            #print [i for i,e in enumerate(f_events) if e == max(f_events)] 

            # TODO: Redo to try more outlier window based analysis before going into max in middle case (no window). Cases: 3+ windows, 2 windows, 1 window, and no window.
            try:
                #template_i = max(enumerate(f_events),key=lambda x: x[1])[0]
                outliers = [i for i,e in enumerate(f_events) if e >= outlier_thresh] + [len(f_events)]
                #print "Outliers"
                #print outliers

                #outliers_i = np.argmax(np.diff(outliers))
                outliers_window = merge_window(outliers, 100)
                #print outliers_window

                mid_outliers = find_closest_interval(outliers_window, len(f_events)/2)
                #print mid_outliers

                fp = [s for s in fingerprint_events(f_events[(outliers_window[0][1]+1):(mid_outliers[0])], plot_spectro_name="fp_events_%s.png"%(fname), plot_name="fp_peaks_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir)]
                comp_fp = [s for s in fingerprint_events(f_events[(mid_outliers[1]+1):outliers_window[-1][0]], plot_spectro_name="fp_events_comp_%s.png"%(fname), plot_name="fp_peaks_comp_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir)]
                outlier_window_count += 1
            except IndexError:
                max_event_i = [a[0] for a in enumerate(f_events) if a[1] == max(f_events[int(len(f_events)/4):int(3*len(f_events)/4)]) and a[0] >= int(len(f_events)/4) and a[0] < int(3*len(f_events)/4)]   
                #print max_event_i
                max_event_i = max_event_i[0]  # use first potential max in case there are multiple in the range
                fp = [s for s in fingerprint_events(f_events[:max_event_i], plot_spectro_name="fp_events_%s.png"%(fname), plot_name="fp_peaks_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir)]
                comp_fp = [s for s in fingerprint_events(f_events[(max_event_i+1):], plot_spectro_name="fp_events_comp_%s.png"%(fname), plot_name="fp_peaks_comp_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir)]
                index_error_count += 1
                
            #r.dtw
            g_id_match, count, diff_counter = find_matches(fp, fingerprint_db, unique)
            g_id_match_comp, count_comp, diff_counter_comp = find_matches(comp_fp, revcomp_fingerprint_db, unique)

            if overlap:
                geno_id, largest_count, largest_offset, diff_overlap = find_overlap_fp(diff_counter, diff_counter_comp)
                for g_id in geno_id:
                    genome_counts[g_id] += 1. / len(geno_id)
            elif not unique:
                for g_id in g_id_match:
                    genome_counts[g_id] += 1. / len(g_id_match)
                for g_id in g_id_match_comp:
                    genome_counts[g_id] += 1. / len(g_id_match_comp)
            else:
                if len(g_id_match) == 1 and len(g_id_match_comp) == 1:
                    if g_id_match[0] == g_id_match_comp[0]:
                        genome_counts[g_id_match[0]] += 1

    results = sorted(genome_counts.items(), key=itemgetter(1), reverse=True)
    reads_name = input_dir.split("/")[-2:]
    if len(reads_name[-1]) > 0:
        reads_name = reads_name[-1]
    else:
        reads_name = reads_name[0]
    outname = output_dir+"/"+"_".join([reads_name, str(overlap), str(outlier_thresh), str(noverlap), str(nfft), str(fs), str(unique)])

    with open(outname + "_results.txt", "w") as f:
        lines = ["\t".join([str(genome), str(score)]) for genome, score in results]
        f.write("\n".join(lines))

    with open(outname + "_read_case.txt", "w") as f:
        f.write("index_error_count\t" + str(index_error_count) + "\n")
        f.write("outlier_window_count\t" + str(outlier_window_count) + "\n")

    pickle.dump(genome_counts, open(outname+"_genome_counts.p", "wb"))
    pickle.dump(diff_counter, open(outname+"_diff_counter.p", "wb"))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ind = np.arange(len(results))
    genome_names = [".".join(r[0].split(".")[:-1]) for r in results]
    bars = ax.bar(ind, [r[1] for r in results], align='center')
    ax.set_ylabel("Count of Reads Aligned")
    ax.set_xlabel("Genome")
    xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in bars]
    plt.xticks(ind, genome_names,  ha='right', rotation=45)
    plt.gcf().subplots_adjust(bottom=0.45)
    fig.savefig(outname + "_results.png")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ind = np.arange(len(results))
    genome_names = [".".join(r[0].split(".")[:-1]) for r in results]
    bars = ax.bar(ind, [100*r[1]/sum(genome_counts.values()) for r in results], align='center')
    ax.set_ylabel("Percent of Reads Aligned")
    ax.set_xlabel("Genome")
    xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in bars]
    plt.xticks(ind, genome_names,  ha='right', rotation=45)
    plt.gcf().subplots_adjust(bottom=0.45)
    fig.savefig(outname + "_results_percent.png")


    return (genome_counts, diff_counter)


def main(argv):
    # Parse arguments
    try:
        opts, args = getopt.getopt(argv, "", ["input_dir=", "output_dir=", "genomes=", "overlap=", "outlier_thresh=", "noverlap=", "nfft=", "fs=", "unique=", "create_db=", "plot_spectro=", "plot=", "only_genomes", "num_processes="])
    except getopt.GetoptError as err:
        print str(err)
        sys.exit(2)
    options = dict(opts)
    if "--create_db" in options:
        create_db = options['--create_db'] == "True"
    else:
        create_db = True
    if "--unique" in options:
        unique = options['--unique'] == "True"
    else:
        unique = False
    if "--fs" in options:
        fs = int(options['--fs'])
    else:
        fs = 2
    if "--nfft" in options:
        nfft = int(options['--nfft'])
    else:
        nfft = 128
    if "--noverlap" in options:
        noverlap = int(options['--noverlap'])
    else:
        noverlap = 64
    if "--outlier_thresh" in options:
        outlier_thresh = int(options['--outlier_thresh'])
    else:
        outlier_thresh = 100
    if "--overlap" in options:
        overlap = options['--overlap'] == "True"
    else:
        overlap = True
    if "--num_processes" in options:
        num_processes = int(options['--num_processes'])
    else:
        num_processes = 1
    #print('create_db %s, unique %s, sampling rate %s, fft block size %s, overlap %s, outlier threshold %s' %(str(x) for x in [create_db,unique,fs,nfft,noverlap,outlier_thresh]))
    #print (str(x) for x in [create_db,unique,fs,nfft,noverlap,outlier_thresh])
    print('sampling rate = ' + str(fs)  + ', fft block size = ' + str(nfft) + ', overlap = ' + str(noverlap) + ', outlier threshold = ' +str(outlier_thresh) + ', number of processes = ' + str(num_processes))
    if "--genomes" in options:
        genomes = options['--genomes']
    else:
        genomes = "/zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/escherichia.coli.k12.m1655/escherichia.coli.k12.m1655.fasta,\
                    /zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/entero/sequence.fasta,\
                    /zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/halobacillus.halophilus.dsm.2266/halobacillus.halophilus.dsm2266.fasta,\
                    /zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/micrococcus.luteus.NCTC2665/micrococcus.luteus.NCTC2665.fasta,\
                    /zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/pseudomonas.flourescens/pseudomonas.fluorescens.f113.fasta,\
                    /zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/staphylococcus.epidermidis.pm221/staphylococcus_epidermidis_pm221.fasta"
    if "--plot_spectro" in options:
        plot_spectro = options["--plot_spectro"] == "True"
    else:
        plot_spectro = False
    if "--plot" in options:
        plot = options["--plot"] == "True"
    else:
        plot = False
    input_dir = options['--input_dir']
    output_dir = options['--output_dir']
    if "--only_genomes" in options:
        only_genomes = True
    else:
        only_genomes = False

    genome_counts, diff_counter = run_UNFOG(input_dir, output_dir, create_db, unique, fs, nfft, noverlap, outlier_thresh, overlap, genomes, plot_spectro, plot, only_genomes, num_processes)

    results = sorted(genome_counts.items(), key=itemgetter(1), reverse=True)
    print sorted(genome_counts.items(), key=itemgetter(1), reverse=True)
    print [100*r[1]/sum(genome_counts.values()) for r in results]
    

if __name__ == "__main__":
    main(sys.argv[1:])

"""
sum_count = sum([r[1] for r in results])
fig = plt.figure()
ax = fig.add_subplot(111)
ind = np.arange(len(results))
genome_names = [".".join(r[0].split(".")[:-1]) for r in results]
bars = ax.bar(ind, [100*r[1]/sum_count for r in results], align='center')
ax.set_ylabel("Percent of Reads Aligned")
ax.set_xlabel("Genome")
xticks_pos = [0.65*patch.get_width() + patch.get_xy()[0] for patch in bars]
plt.xticks(ind, genome_names,  ha='right', rotation=45)
plt.gcf().subplots_adjust(bottom=0.45)
plt.gcf().subplots_adjust(left=0.2)
fig.savefig(outname + "_results_percent.png")
"""
