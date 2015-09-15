#!/usr/bin/python

"""
functions for matching fingerprints from a sample to a reference set
"""

from __future__ import division

from collections import Counter
import matplotlib as mpl


def increment_dict_counts(in_dict, update_dict):
    """
    Adds counts stored in two dictionaries, returns the summed dictionary
    """
    ret = in_dict.copy()
    for k in update_dict:
        try:
            ret[k] += update_dict[k]
        except KeyError:
            ret[k] = update_dict[k]
    return ret


def find_matches(samples, geno_db, unique):
    """
    Finds genomes that match a sample read and returns the best one(s) 
    * measures offset for genome - offset for sample for each match
    * makes dict of offset differences for each genome by number of corresponding matches
    * returns the genome(s) with the highest number of matches for any particular offset difference  
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
                #print(offset)
                matches.append((g, offset - mapper[h])) #this doesn't take into account length of sample at all, and what are reasonable distances between samples in terms of absolute genomic positions of matches
    #print(matches)
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
    if len(geno_id) > 1: #added this Sep 4, actually made things worse for unique, no change for non unique for mgrg short test
        ordered_geno = Counter(geno_id).most_common()
        geno_id = [x[0] for x in ordered_geno if x[1] == ordered_geno[0][1]]
        #print (geno_id, largest_count)
    if unique and len(geno_id) >1: #TODO: figure out how to return best of possible matches
        return ([], -1, {})
    return (geno_id, largest_count, diff_counter)


def find_overlap_fp(diff_counter, diff_counter_comp):
    """
    Checks if same offsets present for the forward and reverse strands
    * if yes, increases the counts for those diffs in the diff_counter by the number in the complement 
    * checks if new count for the overlapping offset is greater than the previously recorded largest count
    * returns the new genomes(s) with the highest match count for a particular offset
    """
    offset_overlap = set(diff_counter.keys()).intersection(set(diff_counter_comp.keys()))
    diff_overlap = {}
    largest_count = -1
    largest_offset = -1
    geno_id = []
    for o in offset_overlap:
        diff_overlap[o] = increment_dict_counts(diff_counter[o], diff_counter_comp[o])
        for g in diff_overlap[o]: #TODO change so overlaps only count if same genome  
            if diff_overlap[o][g] > largest_count:
                largest_count = diff_overlap[o][g]
                largest_offset = o
                geno_id = [g]
            elif diff_overlap[o][g] == largest_count:
                geno_id.append(g)
    return (geno_id, largest_count, largest_offset, diff_overlap)

