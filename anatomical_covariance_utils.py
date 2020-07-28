#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:28:56 2019

@author: jonyoung
"""

from nibabel.freesurfer.io import read_annot, read_morph_data
from nibabel.freesurfer.mghformat import load

# do the subject dirs contain FS directory structure?
FS_dir_structure = True

# dictionary listing which label represents background / no ROI in each .annot file
bg_label_dict = {'aparc':-1, 'aparc2009':-1, 'DK':-1, 'HCP':0, 'BN':-1, 'Schaefer_100':0, 
                  'Schaefer_200':0, 'Schaefer_300':0, 'Schaefer_400':0, 'Schaefer_500':0, 
                  'Schaefer_600':0, 'Schaefer_800':0, 'Schaefer_1000':0, 'Power':0, 'hcp-mmp':0,
                  'gordon':0, 'arslan':0, 'shen':0, 'Cambridge':0}

def get_bg_label(atlas) :
    
    return bg_label_dict[atlas]

# import the lh and rh thickness/pial_lgi etc for a subject in dir and return them
def get_native_cortical_data(subject_dir, data_type) :
    
    if FS_dir_structure :
    
        lh_cortical_data = read_morph_data(subject_dir + '/surf/lh.' + data_type)
        rh_cortical_data = read_morph_data(subject_dir + '/surf/rh.' + data_type)
        
    else :
        
        lh_cortical_data = read_morph_data(subject_dir + '/lh.' + data_type)
        rh_cortical_data = read_morph_data(subject_dir + '/rh.' + data_type)
    
    return lh_cortical_data, rh_cortical_data

# import the lh and rh labels for a subject in dir and return them
def get_label_data(subject_dir, annot_type) :
    
    # dictionary turning an annotation type into a filename
    annot_dict = {'aparc':'aparc.annot', 'aparc2009':'aparc.a2009s.annot', 'DK':'aparc.DKTatlas.annot', 
                  'HCP':'HCP-MMP1.annot', 'BN':'BN_Atlas.annot', 'Schaefer_100':'Schaefer2018_100Parcels_7Networks_order.annot', 
                  'Schaefer_200':'Schaefer2018_200Parcels_7Networks_order.annot', 'Schaefer_300':'Schaefer2018_300Parcels_7Networks_order.annot', 
                  'Schaefer_400':'Schaefer2018_400Parcels_7Networks_order.annot', 'Schaefer_500':'Schaefer2018_500Parcels_7Networks_order.annot', 
                  'Schaefer_600':'Schaefer2018_600Parcels_7Networks_order.annot', 'Schaefer_800':'Schaefer2018_800Parcels_7Networks_order.annot', 
                  'Schaefer_1000':'Schaefer2018_1000Parcels_7Networks_order.annot', 'Power':'power.annot', 'hcp-mmp':'hcp-mmp.annot',
                  'gordon':'gordon333dil.annot', 'arslan':'arslan_res347.annot', 'shen':'shen.annot', 'Cambridge':'500.aparc.annot'}
    
    if FS_dir_structure :
        
        lh_annot = read_annot(subject_dir + '/label/lh.' + annot_dict[annot_type])
        rh_annot = read_annot(subject_dir + '/label/rh.' + annot_dict[annot_type])
        
    else :
    
        lh_annot = read_annot(subject_dir + '/lh.' + annot_dict[annot_type])
        rh_annot = read_annot(subject_dir + '/rh.' + annot_dict[annot_type])
    
    return lh_annot[0], rh_annot[0]