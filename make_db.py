#!/usr/bin/python

"""
Functions for building a database of genomes
"""
import os
import math
import multiprocessing
import cPickle as pickle
import glob
from collections import defaultdict
import gzip

from fingerprint_150924 import (comp, fingerprint_events)

def list_to_dict(li):  
    dct = {}  
    for item in li:  
        if item[0] in dct:  
            dct[item[0]].append(item[1])  
        else:  
            dct[item[0]] = [item[1]]  
    return(dct)  

def genomes_to_fingerprint_db(genome_files, model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, outputdir=""):
    """ converts genomes in a list of fasta genome files to fingerprints, returns a dict of reference genomes and associated fingerprint values """

    fingerprint_database = defaultdict(dict)
    fingerprint_database_comphmm = defaultdict(dict)
    revcomp_fingerprint_database = defaultdict(dict)
    revcomp_fingerprint_database_temphmm = defaultdict(dict)
    template_model_dict = dict([(a[0], a[1]) for a in model[0]])
    comp_model_dict = dict([(a[0], a[1]) for a in model[1]])
    nt_bases = set('ATCG')
    genome_ind = 0
    genome_files = list(set(genome_files)) #ensures unique entries in file list
    for g in genome_files:
        genome_ind +=1
        if len(genome_files) < 20:
            print g
        fname = g.split(".gz") 
        if len(fname) > 1: 
            gzf = gzip.open(g, "rb")
            lines = gzf.read().split("\n")
            gzf.close()
            fname = fname[0]
        else:
            fname = fname[0]
            with open(fname, "r") as gf:
                lines = gf.readlines()
        chr_ind = [i for i,s in enumerate(lines) if s.startswith('>')]
        for i in range(len(chr_ind)):
            if i+1 == len(chr_ind):
                seq = "".join(lines[(1+chr_ind[i]):])
            else:
                seq = "".join(lines[(1+chr_ind[i]):chr_ind[i+1]])
            seq = seq.replace('\n', '').upper()
            kmer_temphmm = []
            kmer_comphmm = []
            kmer_comp_comphmm = []
            kmer_comp_temphmm = []
            for i in range(0, len(seq)-4):
                if all((seq_i in nt_bases) for seq_i in seq[i:i+5]):
                    kmer_temphmm.append(template_model_dict[seq[i:i+5]])
                    kmer_comphmm.append(comp_model_dict[seq[i:i+5]])
                    kmer_comp_comphmm.append(comp_model_dict[comp(seq[i:i+5])])
                    kmer_comp_temphmm.append(template_model_dict[comp(seq[i:i+5])])
        kmer_comp_comphmm.reverse()  # reverse in order to keep same offset as the template
        kmer_comp_temphmm.reverse()
        fp = fingerprint_events(kmer_temphmm, plot_spectro_name="fp_events_%s.png"%(fname), plot_name="fp_peaks_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)
        g_fp = {key:value for (key, value) in fp} # if key not in g_fp else key.append[value]}
        #g_fp = list_to_dict(fp)

        fp_comphmm = fingerprint_events(kmer_comphmm, plot_spectro_name="fp_events_comphmm_%s.png"%(fname), plot_name="fp_peaks_comphmm_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)
        g_fp_comphmm = {key:value for (key,value) in fp_comphmm}
        #g_fp_comphmm = list_to_dict(fp_comphmm)

        fp_comp = fingerprint_events(kmer_comp_comphmm, plot_spectro_name="fp_events_comp_%s.png"%(fname), plot_name="fp_peaks_comp_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)
        g_fp_comp = {key:value for (key, value) in fp_comp}
        #g_fp_comp = list_to_dict(fp_comp)

        fp_comp_temphmm = fingerprint_events(kmer_comp_temphmm, plot_spectro_name="fp_events_comp_temphmm_%s.png"%(fname), plot_name="fp_peaks_comp_temphmm_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)
        g_fp_comp_temphmm = {key:value for (key, value) in fp_comp_temphmm}
        #g_fp_comp_temphmm = list_to_dict(fp_comp_temphmm)
        if only_genomes:
            pickle.dump(g_fp, open(outputdir+fname+"_genome_forward.p", "wb"))
            pickle.dump(g_fp_comphmm, open(outputdir+fname+"_genome_forward_comphmm.p", "wb"))
            pickle.dump(g_fp_comp, open(outputdir+fname+"_genome_complement.p", "wb"))
            pickle.dump(g_fp_temphmm, open(outputdir+fname+"_genome_complement_temphmm.p", "wb"))
        else: 
            fingerprint_database[fname] = g_fp
            fingerprint_database_comphmm[fname] = g_fp_comphmm
            revcomp_fingerprint_database[fname] = g_fp_comp
            revcomp_fingerprint_database_temphmm[fname] = g_fp_comp_temphmm
    return (fingerprint_database, fingerprint_database_comphmm, revcomp_fingerprint_database, revcomp_fingerprint_database_temphmm)


def mp_db_builder(genomes, nprocs, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, outputdir):
    """ distributes list of genomes across processes then adds resulting fingerprints to database """
    def worker(genomes, out_q, hmm_model=hmm_model, fs=fs, nfft=nfft, noverlap=noverlap, only_genomes=only_genomes, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time, outputdir=outputdir):
        outtup  = genomes_to_fingerprint_db(genomes, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, outputdir)
        out_q.put(outtup)
        
    nproc_min = min(nprocs, len(genomes))
    fp_dict = defaultdict(dict)
    fp_dict_comp = defaultdict(dict)
    revcomp_fp_dict = defaultdict(dict)
    revcomp_fp_dict_temp = defaultdict(dict)

    if nproc_min == 1:
        tmp_dict_f_temp, tmp_dict_f_comp,  tmp_dict_r_comp, tmp_dict_r_temp  = genomes_to_fingerprint_db(genomes, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, min_peak_amplitude, outputdir)
        fp_dict.update(tmp_dict_f_temp)
        fp_dict_comp.update(tmp_dict_f_comp)
        revcomp_fp_dict.update(tmp_dict_r_comp)
        revcomp_fp_dict_temp.update(tmp_dict_r_temp)

    else:
        out_q = multiprocessing.Queue()
        chunksize = int(math.ceil(len(genomes) / float(nprocs)))
        procs = []

        for i in range(nproc_min):
            p = multiprocessing.Process(
                    target=worker,
                    args=(genomes[chunksize * i:chunksize * (i + 1)],
                          out_q, hmm_model))
            procs.append(p)
            p.start()

        # Collect all results into dicts for forward and reverse strands. 
        fp_dict = defaultdict(dict) 
        fp_dict_comp = defaultdict(dict)
        revcomp_fp_dict = defaultdict(dict)
        revcomp_fp_dict_temp = defaultdict(dict)
        for proc in procs:
            tmp_dict_f, tmp_dict_f_comp, tmp_dict_r, tmp_dict_r_temp = out_q.get()
            fp_dict.update(tmp_dict_f)
            fp_dict_comp.update(tmp_dict_f_comp)
            revcomp_fp_dict.update(tmp_dict_r)
            revcomp_fp_dict_temp.update(tmp_dict_r_temp)

    # Wait for all worker processes to finish
        for p in procs:
            p.join()
   
    print('Database created')
    return(fp_dict,fp_dict_comp,revcomp_fp_dict,revcomp_fp_dict_temp)


def run_db_builder(hf, db_builder_input, output_dir):
    """ wrapper build function to call from main script """
    if os.path.isfile(output_dir+"/hmm_model.p"):
        hmm_model = pickle.load(open(output_dir+"/hmm_model.p","rb"))
    else:
        try: 
            hmm_model = [hf["Analyses/Basecall_2D_000/BaseCalled_template/Model"].value]
            hmm_model.append(hf["Analyses/Basecall_2D_000/BaseCalled_complement/Model"].value)
        except KeyError:
            return({},{},{},{})
        pickle.dump(hmm_model, open(output_dir+"/hmm_model.p", "wb"))
    genome_dir = output_dir + "genomes/"
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
    genome_files_list, num_processes, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time = db_builder_input
    fingerprint_db, fp_db_comphmm, revcomp_fingerprint_db, revcomp_fp_db_temphmm = mp_db_builder(genome_files_list, num_processes, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, genome_dir)
    return(fingerprint_db, fp_db_comphmm, revcomp_fingerprint_db, revcomp_fp_db_temphmm)

