#!/usr/bin/python

"""
Uncovering Nanopore's Fingerprints of Genomes (UNFOG)

AKA: Shazam for nanopore data.
Based off of https://github.com/worldveil/dejavu.

Angela M Yu, 6/15/2015 - 7/3/2015, Alexa BR McIntyre, 9/1/2015 - 
"""

from __future__ import division

from fingerprint_150924 import fingerprint_events
from make_db import run_db_builder
from match_fp_150924 import (find_matches, find_overlap_fp)
from find_hairpin import define_strands

import h5py
import numpy as np
import cPickle as pickle #pickle in Python3
import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import pylab as P
import re
from operator import itemgetter
import sys
import os
import getopt

def run_UNFOG(input_dir, output_dir, create_db, unique, fs, nfft, noverlap, outlier_thresh, overlap, genomes, plot_spectro, plot, only_genomes, num_processes, min_peak_amplitude, peak_fan, max_hash_time):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plot_dir = output_dir + "/plot/"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    database_open = False
    read_files = glob.glob(input_dir + "/*fast5")
    num_reads = len(read_files)
    print "Number of read files: " + str(num_reads)
    genome_counts = defaultdict(float)
    index_error_count = 0
    outlier_window_count = 0
    if genomes[-1] == "/":
        genome_files_list = glob.glob(genomes + "*fasta*")
    else:
        genome_files_list = genomes.split(",")
    print "Number of Genomes: " + str(len(genome_files_list))
    reads_out_of_bounds = 0
    reverse_strand_mapped = 0
    parameter_suffix = str(fs)+"_"+str(nfft)+"_"+str(noverlap)+"_"+str(len(genome_files_list))+"_"+str(min_peak_amplitude)+"_"+str(peak_fan)+"_"+str(max_hash_time)
    for f in read_files:
        fname = re.findall("([^/]+\.fast5)", f)[0] 
        try:
            with h5py.File(f, "r") as hf:
                if not database_open: 
                    if create_db: 
                        db_builder_input = (genome_files_list, num_processes, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time)
                        fingerprint_db, fp_db_comphmm, revcomp_fingerprint_db, revcomp_fp_db_temphmm = run_db_builder(hf, db_builder_input, output_dir)
                        if fingerprint_db and fp_db_comphmm and revcomp_fingerprint_db and revcomp_fp_db_temphmm:
                            print(f)
                            database_open = True
                        else: 
                            continue #skip to next in read_files list if hmm model files don't exist for this read
                        pickle.dump(fingerprint_db, open(output_dir+"/fingerprint_db_"+parameter_suffix+".p", "wb"))
                        pickle.dump(fp_db_comphmm, open(output_dir+"/fingerprint_db_comphmm_"+parameter_suffix+".p", "wb"))
                        pickle.dump(revcomp_fingerprint_db, open(output_dir+"/revcomp_fingerprint_db_"+parameter_suffix+".p", "wb"))
                        pickle.dump(revcomp_fp_db_temphmm, open(output_dir+"/revcomp_fingerprint_db_temphmm_"+parameter_suffix+".p", "wb"))
                        if only_genomes:
                            return({},[])
                    else:
                        try: 
                            print("Opening database at " + output_dir+"/fingerprint_db_"+parameter_suffix+".p")
                            fingerprint_db = pickle.load(open(output_dir+"/fingerprint_db_"+parameter_suffix+".p", "rb"))
                            fp_db_comphmm = pickle.load(open(output_dir+"/fingerprint_db_comphmm_"+parameter_suffix+".p", "rb"))
                            revcomp_fingerprint_db = pickle.load(open(output_dir+"/revcomp_fingerprint_db_"+parameter_suffix+".p", "rb"))
                            revcomp_fp_db_temphmm = pickle.load(open(output_dir+"/revcomp_fingerprint_db_temphmm_"+parameter_suffix+".p", "rb"))
                            database_open = True
                            print('Database opened')
                        except IOError:
                            print("Reference fingerprint database doesn't yet exist! Try --create_db True")
                            raise SystemExit

                # Compare the fingerprints from the read to what is found in the database
                read_key = (hf["Analyses/EventDetection_000/Reads"].keys())[0]
                if hf.attrs.items()[0][1] == 0.5: #check Metrichor file version
                    f_events = [a[2] for a in hf["Analyses/EventDetection_000/Reads"][read_key]["Events"].value]
                else:
                    f_events = [a[0] for a in hf["Analyses/EventDetection_000/Reads"][read_key]["Events"].value]
                if len(f_events) < 100:
                    reads_out_of_bounds += 1
                    continue

                start_index_template = 20 #hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("start_index_temp")
                end_index_comp = len(f_events)-20 #hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("end_index_comp")
                mid_outliers = [int(len(f_events)/2)-50,int(len(f_events)/2)+50] #[hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("end_index_temp"), hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("start_index_comp")]
                template_events, comp_events, index_error_count = define_strands(f_events, outlier_thresh, index_error_count, start_index_template, mid_outliers, end_index_comp)
               
                fp = [s for s in fingerprint_events(template_events, plot_spectro_name="fp_events_%s.png"%(fname), plot_name="fp_peaks_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)]
                comp_fp = [s for s in fingerprint_events(comp_events, plot_spectro_name="fp_events_comp_%s.png"%(fname), plot_name="fp_peaks_comp_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)]
                
                g_id_match, count, diff_counter = find_matches(fp, fingerprint_db, unique)
                g_id_match_comp, count_comp, diff_counter_comp = find_matches(comp_fp, revcomp_fingerprint_db, unique)
                rev_g_id_match, rev_count, rev_diff_counter = find_matches(fp, revcomp_fp_db_temphmm, unique)
                rev_g_id_match_comp, rev_count_comp, rev_diff_counter_comp = find_matches(comp_fp, fp_db_comphmm, unique)
                if num_reads == 1: #print extra info
                    print(len(f_events))
                    #print(max(f_events))
                    #print [i for i,e in enumerate(f_events) if e == max(f_events)]
                    print(start_index_template,mid_outliers[0])
                    print(hf["Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"].value)
                    print(diff_counter)
                    print(diff_counter_comp)
                    print(rev_diff_counter)
                    print(rev_diff_counter_comp)

                if overlap:
                    geno_id, largest_count, largest_offset, diff_overlap = find_overlap_fp(diff_counter, diff_counter_comp)
                    rev_geno_id, rev_largest_count, rev_largest_offset, rev_diff_overlap = find_overlap_fp(rev_diff_counter, rev_diff_counter_comp)
                    if rev_largest_count > largest_count:
                        reverse_strand_mapped += 1
                        geno_id, largest_count, largest_offset, diff_overlap = rev_geno_id, rev_largest_count, rev_largest_offset, rev_diff_overlap
                    for g_id in geno_id:
                        genome_counts[g_id] += 1. / len(geno_id)
                elif unique:
                    if len(g_id_match) == 1 and len(g_id_match_comp) == 1:
                        if g_id_match[0] == g_id_match_comp[0]:
                            genome_counts[g_id_match[0]] += 1
                else: 
                    for g_id in g_id_match:
                        genome_counts[g_id] += 1. / len(g_id_match)
                    for g_id in g_id_match_comp:
                        genome_counts[g_id] += 1. / len(g_id_match_comp)
        except IOError:
            print("error loading read")
            continue
    print(str(reads_out_of_bounds) + " out of bounds")
    print(str(reverse_strand_mapped) + " mapped to reverse strand")
    results = sorted(genome_counts.items(), key=itemgetter(1), reverse=True)

    reads_name = input_dir.split("/")[-2:]
    if len(reads_name[-1]) > 0:
        reads_name = reads_name[-1]
    else:
        reads_name = reads_name[0]
    outname = output_dir+"/"+"_".join([reads_name, str(overlap), str(outlier_thresh), str(noverlap), str(nfft), str(fs), str(unique), str(min_peak_amplitude), str(peak_fan), str(max_hash_time)])

    with open(outname + "_results.txt", "w") as f:
        lines = ["\t".join([str(genome), str(score)]) for genome, score in results]
        f.write("\n".join(lines))

    with open(outname + "_read_case.txt", "w") as f:
        f.write("index_error_count\t" + str(index_error_count) + "\n")
        f.write("outlier_window_count\t" + str(outlier_window_count) + "\n")

    pickle.dump(genome_counts, open(outname+"_genome_counts.p", "wb"))

    diff_list = []

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ind = np.arange(len(results))
    genome_names = [" ".join(r[0].split("/")[-2:]) for r in results]
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

    print(index_error_count)

    return (genome_counts, diff_counter)


def main(argv):
    # Parse arguments
    try:
        opts, args = getopt.getopt(argv, "", ["input_dir=", "output_dir=", "genomes=", "overlap=", "outlier_thresh=", "noverlap=", "nfft=", "fs=", "unique=", "create_db=", "plot_spectro=", "plot=", "only_genomes", "num_processes=", "min_peak_amplitude=", "peak_fan=", "max_hash_time="])
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
    if "--min_peak_amplitude" in options:
        min_peak_amplitude = int(options['--min_peak_amplitude'])
    else:
        min_peak_amplitude = 10
    if "--peak_fan" in options:
        peak_fan = int(options['--peak_fan'])
    else:
        peak_fan = 15
    if "--max_hash_time" in options:
        max_hash_time = int(options['--max_hash_time'])
    else:
        max_hash_time = 200
    print('sampling rate = ' + str(fs)  + ', fft block size = ' + str(nfft) + ', overlap = ' + str(noverlap) + ', outlier threshold = ' +str(outlier_thresh) + ', number of processes = ' + str(num_processes) + ', minimum FFT peak amplitude = ' + str(min_peak_amplitude) + ', peak fan width = ' + str(peak_fan) + ', maximum hash time = ' + str(max_hash_time))
    if "--genomes" in options:
        genomes = options['--genomes']
    else:
        print 'must specify location of genomes following "--genomes" flag'
        sys.exit(2)
    if "--plot_spectro" in options:
        plot_spectro = options["--plot_spectro"] == "True"
    else:
        plot_spectro = False
    if "--plot" in options:
        plot = options["--plot"] == "True"
    else:
        plot = False
    input_dir = options['--input_dir']
    print('input: ' + input_dir)
    output_dir = options['--output_dir']
    if "--only_genomes" in options:
        only_genomes = True
    else:
        only_genomes = False

    genome_counts, diff_counter = run_UNFOG(input_dir, output_dir, create_db, unique, fs, nfft, noverlap, outlier_thresh, overlap, genomes, plot_spectro, plot, only_genomes, num_processes, min_peak_amplitude, peak_fan, max_hash_time)

    results = sorted(genome_counts.items(), key=itemgetter(1), reverse=True)
    print sorted(genome_counts.items(), key=itemgetter(1), reverse=True)
    print [100*r[1]/sum(genome_counts.values()) for r in results]
    

if __name__ == "__main__":
    main(sys.argv[1:])

