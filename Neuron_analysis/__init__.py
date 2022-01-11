'''Package for dealing with single neuron annotated points in corresponding atlases

'''

import os
import re

import numpy as np
import pandas as pd


import tkinter as tk
import tkinter.filedialog as fdialog
from tkinter import simpledialog

__all__ = ["atlas", "points", "analysis_tools"]

global atlas_labels
atlas_labels=atlas_labels=pd.read_csv('labels.csv')

global fullatlas_name
#fullatlas_name='D:\Allenbrainatlas\ARA_25_micron_mhd_ccf2017\\annotation_25.mhd'

bregma= [227.5, 0,  215]
# Bregma location of 25um atlas identified by me in x y z (points)
# corresponds with ML, DV, AP

def find_mousename(a_str):
    '''find the mouse name from a string.
    example D:/AL142 will return AL142'''
    m=re.search('\D{2}[0-9]{3}', a_str)
    return m[0]

def find_crop(name):
    ''' takes in atlas or tempalte name that has xxx-xxx in its name.
    finds the missing section in z coronal plane to transform back to original full span(1-528)
    For exmaple, for atlas_104-400.mhd, will output 103,128'''
    m= re.search("[0-9]{1,3}.[0-9]{1,3}",name)[0]
    to_add= m.split('-')
    lead= int(to_add[0])-1
    trail= 528 - int(to_add[1])+1-1
    return lead, trail

def get_pt_natlas(dspoint_name,outdir, full=False):
    '''Read the downsampled points and get the corresponding atlas name
    '''
    with open(dspoint_name,'r') as output:
        outputpoint= output.readlines()
    
    all_points=[]
    if 'outputpoints' in dspoint_name:
        files=os.listdir(outdir)
        atlas=[i for i in files if re.search("atlas.+\.mhd",i)][0]
        atlas_name=os.path.join(outdir, atlas)

        for lines in outputpoint:
            m=re.search("(?:OutputIndexFixed = \[ )([0-9]+ [0-9]+ [0-9]+)", lines).groups(0)
            this_line= str(m[0]).split(' ')
            mypoints= [int(stuff) for stuff in this_line]
            all_points.append(mypoints)
    else:
        atlas_name=outdir+'/ara2sample_atlas/result.mhd'
        for lines in outputpoint[2:]:
            this_line= lines.split (' ')
            mypoints= [int(stuff) for stuff in this_line]
            all_points.append(mypoints)
    
    if full== True:
        atlas_name=fullatlas_name
    else: 
        pass
    return all_points, atlas_name

def give_me_name(x):
    ''' return the name of brain region give the provided id
    will return an empty string if the id is not found in atlas label'''
    found=atlas_labels[atlas_labels.id==x]
    if found.size ==0:
        out=' '
    else:
        out=found.acronym.item()
    return out

def stereotaxis(data_point, ml=1):
    '''Converts data value(1 value) to stereotaxic coordinates for 25um/pixel scale
    outputs ML and AP positions in mm to the nearest 2 places
    ml=1 when computing medial-lateral points(ie, x coordinate in a coronal view)
    ml=0 when computing anterior posterior points(ie. z coordinate in a coronal view)
    '''
    
    if ml==1:
        out=  abs(bregma[0] - data_point) *25 /1000
    else:
        out= (bregma[2]- data_point) *25 /1000
    return round(out,2)


    
    
    