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

from fingerprint import (comp, fingerprint_events)


def genomes_to_fingerprint_db(genome_files, model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, outputdir=""):
    """ converts genomes in a list of fasta genome files to fingerprints, returns a dict of reference genomes and associated fingerprint values """
    #TODO: Implement a true database, not a dict of dicts which was used for being able to quickly test the Shazam concept.
    fingerprint_database = defaultdict(dict)
    revcomp_fingerprint_database = defaultdict(dict)
    template_model_dict = dict([(a[0], a[1]) for a in model[0]])
    comp_model_dict = dict([(a[0], a[1]) for a in model[1]])
    nt_bases = set('ATCG')
    genome_ind = 0
    genome_files = list(set(genome_files)) #ensures unique entries in file list
    for g in genome_files:
        genome_ind +=1
        #if genome_ind%100 == 0:
        #    print genome_ind
        if len(genome_files) < 20:
            print g
        fname = g.split(".gz") #re.findall("([^/]+)\.f[^./]*a(.gz)?", g)[0]
        if len(fname) > 1: # and (len(fname) != 2 and '' not in fname): #FIX THIS
            gzf = gzip.open(g, "rb")
            lines = gzf.read().split("\n")
            gzf.close()
            fname = fname[0]
        else:
            fname = fname[0]
            with open(fname, "r") as gf: #MOD FROM open(g, "r")
                lines = gf.readlines()
        chr_ind = [i for i,s in enumerate(lines) if s.startswith('>')]
        for i in range(len(chr_ind)):
            if i+1 == len(chr_ind):
                seq = "".join(lines[(1+chr_ind[i]):])
            else:
                seq = "".join(lines[(1+chr_ind[i]):chr_ind[i+1]])
            seq = seq.replace('\n', '').upper()
            kmer = []
            kmer_comp = []
            for i in range(0, len(seq)-4):
                if all((seq_i in nt_bases) for seq_i in seq[i:i+5]):
                    kmer.append(template_model_dict[seq[i:i+5]])
                    kmer_comp.append(comp_model_dict[comp(seq[i:i+5])])
        kmer_comp.reverse()  # reverse in order to keep same offset as the template
        """if g == '/zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/entero/sequence.fasta':
            fname = 'enterobacter'
            genome_event_file = open(outputdir+fname+"_genome_reverse_events","w")
            genome_event_file.write(','.join([str(k) for k in kmer_comp]))
            genome_event_file.close()
            break""" 
        #print str(kmer[0]) + ': kmer example'
        fp = fingerprint_events(kmer, plot_spectro_name="fp_events_%s.png"%(fname), plot_name="fp_peaks_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)
        g_fp = {key:value for (key, value) in fp}
        fp_comp = fingerprint_events(kmer_comp, plot_spectro_name="fp_events_comp_%s.png"%(fname), plot_name="fp_peaks_comp_%s.png"%(fname), fs=fs, nfft=nfft, noverlap=noverlap, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time)
        g_fp_comp = {key:value for (key, value) in fp_comp}
        #if g == '/zenodotus/masonlab/noa2019_scratch/nanopore/alignment/genomes/entero/sequence.fasta':
        if only_genomes:
            #fname = 'enterobacter'
            pickle.dump(g_fp, open(outputdir+fname+"_genome_forward.p", "wb"))
            pickle.dump(g_fp_comp, open(outputdir+fname+"_genome_complement.p", "wb"))
        else: 
            fingerprint_database[fname] = g_fp
            revcomp_fingerprint_database[fname] = g_fp_comp
    return (fingerprint_database, revcomp_fingerprint_database)


def mp_db_builder(genomes, nprocs, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, outputdir):
    """ distributes list of genomes across processes then adds resulting fingerprints to database """
    def worker(genomes, out_q, hmm_model=hmm_model, fs=fs, nfft=nfft, noverlap=noverlap, only_genomes=only_genomes, min_peak_amplitude=min_peak_amplitude, peak_fan=peak_fan, max_hash_time=max_hash_time, outputdir=outputdir):
        outtup  = genomes_to_fingerprint_db(genomes, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, outputdir)
        out_q.put(outtup)
        
    nproc_min = min(nprocs, len(genomes))
    fp_dict = defaultdict(dict)
    revcomp_fp_dict = defaultdict(dict)
    if nproc_min == 1:
        tmp_dict_f, tmp_dict_r = genomes_to_fingerprint_db(genomes, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, min_peak_amplitude, outputdir)
        fp_dict.update(tmp_dict_f)
        revcomp_fp_dict.update(tmp_dict_r)
        
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
        # TODO: check single instances of genomes in genome list
        fp_dict = defaultdict(dict) 
        revcomp_fp_dict = defaultdict(dict)
        for proc in procs:
            tmp_dict_f, tmp_dict_r = out_q.get()
            fp_dict.update(tmp_dict_f)
            revcomp_fp_dict.update(tmp_dict_r)

    # Wait for all worker processes to finish
        for p in procs:
            p.join()
   
    print('Database created')
    #print(fp_dict.keys())
    #print(revcomp_fp_dict.keys())
    return(fp_dict,revcomp_fp_dict)


def run_db_builder(hf, db_builder_input, output_dir):
    """ wrapper build function to call from main script """
    if os.path.isfile(output_dir+"/hmm_model.p"):
        hmm_model = pickle.load(open(output_dir+"/hmm_model.p","rb"))
    else:
        try: 
            hmm_model = [hf["Analyses/Basecall_2D_000/BaseCalled_template/Model"].value]
            hmm_model.append(hf["Analyses/Basecall_2D_000/BaseCalled_complement/Model"].value)
        except KeyError:
            return({},{})
        pickle.dump(hmm_model, open(output_dir+"/hmm_model.p", "wb"))
    genome_dir = output_dir + "genomes/"
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
    genome_files_list, num_processes, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time = db_builder_input
    fingerprint_db, revcomp_fingerprint_db = mp_db_builder(genome_files_list, num_processes, hmm_model, fs, nfft, noverlap, only_genomes, min_peak_amplitude, peak_fan, max_hash_time, genome_dir)
    return(fingerprint_db, revcomp_fingerprint_db)

