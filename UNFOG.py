#!/usr/bin/python

"""
Uncovering Nanopore's Fingerprints of Genomes (UNFOG)

AKA: Shazam for nanopore data.
Based off of https://github.com/worldveil/dejavu.

Angela M Yu, 6/15/2015 - 7/3/2015, Alexa BR McIntyre, 9/1/2015 - 
"""

from __future__ import division

from fingerprint import fingerprint_events
from make_db import run_db_builder
from match_fp import (find_matches, find_overlap_fp)
from find_hairpin import define_strands

#import multiprocessing
#import math
import h5py
import numpy as np
import cPickle as pickle #pickle in Python3
import glob
from collections import defaultdict
import matplotlib.pyplot as plt
import re
from operator import itemgetter
#import rpy2.robjects.numpy2ri
#rpy2.robjects.numpy2ri.activate()
#from rpy2.robjects.packages import importr
import sys
import os
import getopt

##from rpy2.robjects import r


def run_UNFOG(input_dir, output_dir, create_db, unique, fs, nfft, noverlap, outlier_thresh, overlap, genomes, plot_spectro, plot, only_genomes, num_processes, min_peak_amplitude, peak_fan, max_hash_time):
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
        try:
            with h5py.File(f, "r") as hf:
                #print(hf["Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"].value)
                #print(hf["Analyses/Basecall_2D_000/BaseCalled_2D/Alignment"])
                if not database_open: 
                    if create_db: 
                    # runs only for first read file, assumes that the same model is applied to all the reads analyzed - this is more or less true according to nanoporetech wiki for different pores, not sure about different channels and it seems to depend on the version of the workflow used   
                        db_builder_input = (genome_files_list, num_processes, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time)
                        fingerprint_db, revcomp_fingerprint_db = run_db_builder(hf, db_builder_input, output_dir)
                        # TODO: implement f["Analyses"]["Basecall_2D_000"]["BaseCalled_complement/Model"].attrs.items()?
                        if fingerprint_db and revcomp_fingerprint_db:
                            print(f)
                            database_open = True
                        else: 
                            continue #skip to next in read_files list if hmm model files don't exist for this read
                        pickle.dump(fingerprint_db, open(output_dir+"/fingerprint_db_"+str(fs)+"_"+str(nfft)+"_"+str(noverlap)+"_"+str(len(genome_files_list))+"_"+str(min_peak_amplitude)+"_"+str(peak_fan)+"_"+str(max_hash_time)+".p", "wb"))
                        pickle.dump(revcomp_fingerprint_db, open(output_dir+"/revcomp_fingerprint_db_"+str(fs)+"_"+str(nfft)+"_"+str(noverlap)+"_"+str(len(genome_files_list))+"_"+str(min_peak_amplitude)+"_"+str(peak_fan)+"_"+str(max_hash_time)+".p", "wb"))
                        if only_genomes:
                            return({},[])
                    else:
                        try: 
                            fingerprint_db = pickle.load(open(output_dir+"/fingerprint_db_"+str(fs)+"_"+str(nfft)+"_"+str(noverlap)+"_"+str(len(genome_files_list))+"_"+str(min_peak_amplitude)+"_"+str(peak_fan)+"_"+str(max_hash_time)+".p", "rb"))
                            revcomp_fingerprint_db = pickle.load(open(output_dir+"/revcomp_fingerprint_db_"+str(fs)+"_"+str(nfft)+"_"+str(noverlap)+"_"+str(len(genome_files_list))+"_"+str(min_peak_amplitude)+"_"+str(peak_fan)+"_"+str(max_hash_time)+".p", "rb"))
                            database_open = True
                            print('Database opened')
                        except IOError:
                            print("Reference fingerprint database doesn't yet exist! Try --create_db True")

                # Compare the fingerprints from the read to what's in the database
                read_key = (hf["Analyses/EventDetection_000/Reads"].keys())[0]
                f_events = [a[0] for a in hf["Analyses/EventDetection_000/Reads"][read_key]["Events"].value]
                event_times = [a[3] for a in hf["Analyses/EventDetection_000/Reads"][read_key]["Events"].value]
                """if len(f_events) > 1000 and len(f_events) < 10000:
                    event_file = open(output_dir+'/enterobacter_events_example','w')
                    event_file.write(','.join([str(x) for x in f_events]))
                    event_file.close()
                    time_file = open(output_dir+'/enterobacter_event_times_example','w')
                    time_file.write(','.join([str(x) for x in event_times]))
                    time_file.close()
                    break"""
                if np.mean(f_events) > 1000 or max(f_events) < outlier_thresh:
                    continue
                #print len(f_events)
                #print max(f_events)
                #print [i for i,e in enumerate(f_events) if e == max(f_events)] 

                # TODO: Redo to try more outlier window based analysis before going into max in middle case (no window). Cases: 3+ windows, 2 windows, 1 window, and no window.

                start_index_template = hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("start_index_temp")
                end_index_comp = hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("end_index_comp")
                mid_outliers = [hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("end_index_temp"), hf["Analyses/Basecall_2D_000/Summary/split_hairpin"].attrs.get("start_index_comp")]
                template_events, comp_events, index_error_count = define_strands(f_events, outlier_thresh, index_error_count, start_index_template, mid_outliers, end_index_comp)
                #if index_error_count > 0: #TEMP
                #    break
                fp = [s for s in fingerprint_events(template_events, plot_spectro_name="fp_events_%s.png"%(fname), plot_name="fp_peaks_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)]
                comp_fp = [s for s in fingerprint_events(comp_events, plot_spectro_name="fp_events_comp_%s.png"%(fname), plot_name="fp_peaks_comp_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, plot_spectro=plot_spectro, plot=plot, outputdir=plot_dir, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)]
                
                #r.dtw
                g_id_match, count, diff_counter = find_matches(fp, fingerprint_db, unique)
                g_id_match_comp, count_comp, diff_counter_comp = find_matches(comp_fp, revcomp_fingerprint_db, unique)
                #TODO: check if find better matches mapping template to reverse complements of genomes, complements to forward strands
                #if not g_id_match or not g_id_match_comp:
                #    print('rev')
                #g_id_match, count, diff_counter = find_matches(comp_fp, fingerprint_db, unique)
                #g_id_match_comp, count_comp, diff_counter_comp = find_matches(fp, revcomp_fingerprint_db, unique)

               # print(len(diff_counter))
               # print(len(diff_counter_comp))
                if overlap:
                    geno_id, largest_count, largest_offset, diff_overlap = find_overlap_fp(diff_counter, diff_counter_comp)
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
            continue
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
    #pickle.dump(diff_counter, open(outname+"_diff_counter.p", "wb"))

    #TODO: move plotting to its own script
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

