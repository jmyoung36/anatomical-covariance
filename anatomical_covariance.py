#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 15:08:00 2018

@author: jonyoung
"""

from nibabel.freesurfer import load
from nibabel.freesurfer.io import read_annot, read_morph_data
import numpy as np
#from make_sparse_inverse_covariance import get_covariance
from anatomical_covariance_utils import get_native_cortical_data, get_label_data, get_bg_label
from scipy.stats import entropy, gaussian_kde, wasserstein_distance

# metadata dir to find bins
metadata_dir = '/home/jonyoung/IoP_data/Data/EU-GEI/metadata/'

# calculate intersection of two histograms
def histogram_intersection(histogram_1, histogram_2) :
    
    paired_hists = zip(histogram_1, histogram_2)
    intersection = np.sum(list(map(lambda x: np.min(x), paired_hists)))       
    return intersection

def symmetric_KL(histogram_1, histogram_2) :
    
    # remove zeros with smoothing
#    nonzero = np.logical_not(np.logical_or(histogram_1 == 0, histogram_2 == 0))
#    histogram_1 = histogram_1[nonzero]
#    histogram_2= histogram_2[nonzero]
    histogram_1 = histogram_1 + 0.0001
    histogram_1 = histogram_1 / np.sum(histogram_1)
    histogram_2 = histogram_2 + 0.0001
    histogram_2 = histogram_2 / np.sum(histogram_2)
    
    return np.mean((entropy(histogram_1, histogram_2), entropy(histogram_2, histogram_1)))

# calculate an anatomical covariance matrix for a single subjects based on
# pairwise histogram intersection of cortical metric values between within
# all cortical ROIs in the selected atlas
def anatomical_covariance_matrix(subject_dir, annot_type, n_bins, data_type, hemispheres, lh_ROI_labels = None, rh_ROI_labels = None) :
    
    method = 'intersection'
    
    # read in the thickness data and labels for rh and lh
    lh_cortical_data, rh_cortical_data = get_native_cortical_data(subject_dir, data_type)
    lh_annot, rh_annot = get_label_data(subject_dir, annot_type)
    
    # remove label for unknown from annotations and thickness datacortical atlases
    #bg_label = bg_label_dict[annot_type]
    bg_label = get_bg_label(annot_type)
    lh_cortical_data = lh_cortical_data[np.logical_not(lh_annot == bg_label)]
    rh_cortical_data = rh_cortical_data[np.logical_not(rh_annot == bg_label)]
    lh_annot = lh_annot[np.logical_not(lh_annot == bg_label)]
    rh_annot = rh_annot[np.logical_not(rh_annot == bg_label)]

    # get set of labels for lh and rh
    # unless they are passed to the function
    if (lh_ROI_labels is None ) : 
        lh_ROI_labels = np.unique(lh_annot)
    if (rh_ROI_labels is None ) :     
        rh_ROI_labels = np.unique(rh_annot)
          
    if (hemispheres == 'left') | (hemispheres == 'both') :
        
        lh_n_ROIs = len(lh_ROI_labels)
        n_ROIs = lh_n_ROIs
        lh_histograms = np.zeros((lh_n_ROIs, n_bins))
    
        # fill in histograms
        for i in range(lh_n_ROIs) :
            
            ROI_label = lh_ROI_labels[i]
            ROI_thickness_data = lh_cortical_data[lh_annot == ROI_label]
            
            # different ranges for thickness and lgi
            # 5mm is max permitted value for thickness
            # 8 is highest observed value of LGI, rounded up
            if data_type == 'thickness' :
                
                histogram = np.histogram(ROI_thickness_data, bins=n_bins, range=(0, 5), density=False)
                
            elif data_type == 'pial_lgi' :
                
                # get bin edges for 1024
                bin_edges_1024 = np.load(metadata_dir + 'pial_lgi_1024_bins.npy')
                
                # convert to the edges for the desired number of bins
                skip = 1024/n_bins
                bin_edges = bin_edges_1024[0::skip]
                histogram = np.histogram(ROI_thickness_data, bins=bin_edges, range=(0, 3), density=False)
                #histogram = np.histogram(ROI_thickness_data, bins=n_bins, range=(0, 3), density=False)
                
            # in case histogram is empty
            if np.sum(histogram[0]) > 0 :
                histogram = histogram[0]/float(np.sum(histogram[0]))
            else :
                histogram = histogram[0]
            lh_histograms[i, :] = histogram
            
        histograms = lh_histograms
    
    if (hemispheres == 'right') | (hemispheres == 'both') :
        
        rh_n_ROIs = len(rh_ROI_labels)
        n_ROIs = rh_n_ROIs
        rh_histograms = np.zeros((rh_n_ROIs, n_bins))
    
        for i in range(rh_n_ROIs) :
            
            ROI_label = rh_ROI_labels[i]
            ROI_thickness_data = rh_cortical_data[rh_annot == ROI_label]
            # different ranges for thickness and lgi
            # 5mm is max permitted value for thickness
            # 8 is highest observed value of LGI, rounded up
            if data_type == 'thickness' :
                
                histogram = np.histogram(ROI_thickness_data, bins=n_bins, range=(0, 5), density=False)
                
            elif data_type == 'pial_lgi' :
                
                # get bin edges for 1024
                bin_edges_1024 = np.load(metadata_dir + 'pial_lgi_1024_bins.npy')
                
                # convert to the edges for the desired number of bins
                skip = 1024/n_bins
                bin_edges = bin_edges_1024[0::skip]
                histogram = np.histogram(ROI_thickness_data, bins=bin_edges, range=(0, 3), density=False)
                #histogram = np.histogram(ROI_thickness_data, bins=n_bins, range=(0, 3), density=False)
                
            # in case histogram is empty
            if np.sum(histogram[0]) > 0 :
                histogram = histogram[0]/float(np.sum(histogram[0]))
            else :
                histogram = histogram[0]
            rh_histograms[i, :] = histogram
            
        histograms = rh_histograms
    

    if hemispheres == 'both' : 
        
        # stack the sets of histograms together
        histograms = np.vstack((lh_histograms, rh_histograms))
        n_ROIs = lh_n_ROIs + rh_n_ROIs
    
    # fill covariance matrix with pairwise histogram intersections
    anatomical_covariance_matrix = np.zeros((n_ROIs, n_ROIs))
    for i in range(n_ROIs) :
    
        for j in range(i+1) :
        
            histogram_1 = histograms[i, :]
            histogram_2 = histograms[j, :]
            if method == 'intersection' :
                
                cov = histogram_intersection(histogram_1, histogram_2)
            
            elif method == 'EMD' :
                
                cov = wasserstein_distance(histogram_1, histogram_2)
                
            elif method == 'EMD_similarity' :
                
                cov = 1.0 / (1.0 + wasserstein_distance(histogram_1, histogram_2))
                
            #cov = symmetric_KL(histogram_1, histogram_2)
            anatomical_covariance_matrix[i, j] = cov
            anatomical_covariance_matrix[j, i] = cov
            
    return anatomical_covariance_matrix, lh_ROI_labels, rh_ROI_labels

# calculate regional mean of cortical data (thickness or LGI) with specified atlas
def cortical_regional_means(subject_dir, annot_type, data_type, hemisphere, lh_ROI_labels = None, rh_ROI_labels = None) :
    
    # read in the cortical data and labels for rh and lh
    lh_cortical_data, rh_cortical_data = get_native_cortical_data(subject_dir, data_type)
    lh_annot, rh_annot = get_label_data(subject_dir, annot_type)
    
    # remove label for unknown from annotations and data cortical atlases
    #bg_label = bg_label_dict[annot_type]
    bg_label = get_bg_label(annot_type)
    lh_cortical_data = lh_cortical_data[np.logical_not(lh_annot == bg_label)]
    rh_cortical_data = rh_cortical_data[np.logical_not(rh_annot == bg_label)]
    lh_annot = lh_annot[np.logical_not(lh_annot == bg_label)]
    rh_annot = rh_annot[np.logical_not(rh_annot == bg_label)]
    
    # get set of labels for lh and rh
    # unless they are passed to the function
    if (lh_ROI_labels is None ) : 
        lh_ROI_labels = np.unique(lh_annot)
    if (rh_ROI_labels is None ) :     
        rh_ROI_labels = np.unique(rh_annot)
    
    # create data structure for mean
    if (hemisphere == 'left') | (hemisphere == 'both') :
    
        lh_n_ROIs = len(lh_ROI_labels)
        n_ROIs = lh_n_ROIs
        lh_ROI_mean_data = np.zeros((lh_n_ROIs, 1))
        for i in range(lh_n_ROIs) :
        
            ROI_label = lh_ROI_labels[i]
            ROI_data = lh_cortical_data[lh_annot == ROI_label]
            lh_ROI_mean_data[i] = np.mean(ROI_data)
            ROI_mean_data = lh_ROI_mean_data
    
    if (hemisphere == 'right') | (hemisphere == 'both') :
    
        rh_n_ROIs = len(rh_ROI_labels)
        n_ROIs = rh_n_ROIs
        rh_ROI_mean_data = np.zeros((rh_n_ROIs, 1))
        for i in range(rh_n_ROIs) :
        
            ROI_label = rh_ROI_labels[i]
            ROI_data = rh_cortical_data[rh_annot == ROI_label]
            rh_ROI_mean_data[i] = np.mean(ROI_data)
            ROI_mean_data = rh_ROI_mean_data
    
    if hemisphere == 'both' :
        
        ROI_mean_data = np.vstack((lh_ROI_mean_data, rh_ROI_mean_data))
        
    return ROI_mean_data, lh_ROI_labels, rh_ROI_labels
        

    
    
        
        