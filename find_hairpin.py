import sys
import numpy as np
from itertools import izip

def merge_window(indices, window=100):
    """
    finds windows of closely located indices from a list
    """
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
    """
    finds the interval closest to the index i
    """
    # assumes non-overlapping intervals
    closest_i = 0
    diff = sys.maxint
    for ii in range(len(interval_list)):
        s, e = interval_list[ii]
        if i >= s and i <= e:
            closest_i = ii
            break
        d = min(abs(s-i), abs(e-i))  #huh? 
        if d < diff:
            diff = d
            closest_i = ii
    return interval_list[closest_i]

def find_max_delta_index(list_a,list_b,larger,threshold):
    max_delta = max([larger*(y - x) for (x, y) in izip(list_a, list_b)])
    if max_delta > threshold:
        max_delta_ind = [y[0] for y in enumerate(list_b) if larger*(y[1] - list_a[y[0]]) == max_delta]
        return(max_delta_ind[0])
    else:
        return -1


def define_strands(f_events, outlier_thresh, index_error_count, start_index_template, mid_outliers, end_index_comp):
    """
    Finds the hairpin from the central peak surpassing the provided outlier threshold if possible, else finds the maximum central event.
    Returns f_events sorted into template and complementary strands
    """
    try:
        #outlier_thresh = max(int(np.mean(f_events[int(len(f_events)/8):7*int(len(f_events)/8)])*1.15), outlier_thresh)
        #outliers = [i for i,e in enumerate(f_events) if e >= outlier_thresh] + [len(f_events)]
        #outliers_window = merge_window(outliers, 100)
        #print(outliers_window)
        #mid_outliers = find_closest_interval(outliers_window, len(f_events)/2)

        #outlier_thresh = 4
        #start_index_template = max(find_max_delta_index(f_events[:int(len(f_events)/4)], f_events[1:int(len(f_events)/4)+1], -1, outlier_thresh), 0)
        #end_index_comp = min(len(f_events)-find_max_delta_index(f_events[3*int(len(f_events)/4):-1], f_events[3*int(len(f_events)/4)+1:], 1, outlier_thresh), len(f_events))
        #mid_outliers = [int(len(f_events)/4)+find_max_delta_index(f_events[int(len(f_events)/4):int(len(f_events)/2)], f_events[int(len(f_events)/4)+1:int(len(f_events)/2)+1], 1, outlier_thresh)]
        #mid_outliers.append(int(len(f_events)/2)+find_max_delta_index(f_events[int(len(f_events)/2):3*int(len(f_events)/4)], f_events[int(len(f_events)/2)+1:3*int(len(f_events)/4)+1], -1, outlier_thresh))
        #print(start_index_template,mid_outliers,end_index_comp,len(f_events))

        template_events = f_events[start_index_template:(mid_outliers[0])]
        comp_events = f_events[(mid_outliers[1]+1):end_index_comp]
    except TypeError: 
        #max_event_i = [a[0] for a in enumerate(f_events) if a[1] == max(f_events[int(len(f_events)/4):int(3*len(f_events)/4)]) and a[0] >= int(len(f_events)/4) and a[0] < int(3*len(f_events)/4)]
        #max_event_i = max_event_i[0] 
        #try: #change threshold, reattempt window-based analysis
        #    max_event_val = f_events[max_event_i]
        #    outlier_thresh = max_event_val - 15
        #    outliers = [i for i,e in enumerate(f_events) if e >= outlier_thresh] + [len(f_events)]
        #    outliers_window = merge_window(outliers, 100)
        #    mid_outliers = find_closest_interval(outliers_window, len(f_events)/2)
        #    template_events = f_events[(outliers_window[0][1]+1):(mid_outliers[0])]
        #    comp_events = f_events[(mid_outliers[1]+1):outliers_window[-1][0]]
        #except IndexError: #use max in middle  
        index_error_count += 1  
        template_events = [] #f_events[:max_event_i]
        #print(max_event_val, outlier_thresh, f_events)
        comp_events = [] #f_events[(max_event_i+1):]
    return(template_events, comp_events, index_error_count)


