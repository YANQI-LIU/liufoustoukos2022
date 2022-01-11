'''Contains tools for analysis
    convert_anno: Converts annotation from pixels to ums. This is usually the first step for points analysis. Resolution hardcoded to be 0.8 x 0.8 x 5 um/pixel
    downsample_anno: downsamples annotation that has been already converted(by convert_anno) and resampled(in vaa3d). Resample step =1 and goal =25, hardcoded
    get_pt_natlas: Read the downsampled annotation file and identify the associated atlas. This is the starting point of many subsequent analysis
    make_pd: outputs an excel file with brain region and amount of points & endings within
    check_points: check if all points have an index within the atlas labels
    findID_original: assign atlas ID to origional non-downsampled, non-transformed points
    make_tif: make a tif file of all the downsampled, registered point. Nice to have when merge with brain volume.

'''

import os
import re

import numpy as np
import pandas as pd
import SimpleITK as sitk
import warnings

import skimage
from skimage import io


import tkinter as tk
import tkinter.filedialog as fdialog
from tkinter import simpledialog

global atlas_labels
atlas_labels=atlas_labels=pd.read_csv('D:\Allenbrainatlas\ARA_25_micron_mhd_ccf2017\labels.csv')

global fullatlas_name
fullatlas_name='D:\Allenbrainatlas\ARA_25_micron_mhd_ccf2017\annotation_25.mhd'


def find_mousename(dir):
    '''find the mouse name from out dir.
    example D:/AL142 will return AL142'''
    m=re.search('\D{2}[0-9]{3}', dir)
    return m[0]

def read_eswc():
    # read the final version of eswc file generated from vaa3d, after convert to real physical distances and resampled
    # return headings and annotations in list forms
    resample_file=fdialog.askopenfile(title='Select the final converted and resampled eswc file').name 
        
    with open(resample_file,'r') as resampled_anno:
        resampled_anno_data=resampled_anno.readlines()
    # heading is stored in anno_data[2], 1st line basically useless

    headings=resampled_anno_data[2].rstrip('\n').replace(' ', '').split(',')
    annotations=[lines[0:-5].split(' ') for lines in resampled_anno_data[3:]]
    #slight modification on replacing and stripping due to the format of the resampled swc

    return headings, annotations

def downsample_anno():
    ''' Downsample the annotation that has been already converted to um and resampled
        Creates two files: main downsample file and list of indexes for endings
    '''
    outdir = fdialog.askdirectory(initialdir='D:/', title='Please select the output directory')
    resample_file=fdialog.askopenfile(initialdir='M:/analysis/Yanqi_Liu/Annotations/', title='Select the eswc file containing the converted and resampled annotations').name 
    
    resampled_xyz=1
    goal_xyz=25
    #resampled_xyz = simpledialog.askfloat("Input", "What is the x, y and z resolution in um after converting and resampling?",
    #                              minvalue=0.0, maxvalue=100)
    #goal_xyz = simpledialog.askfloat("Input", "What do you want to downsample the resolution to '(in um)' ?",
    #                              minvalue=10, maxvalue=100)
    name=find_mousename(outdir)

    ratioxyz=goal_xyz/resampled_xyz

    message = (f"Resampled annotation step size is {resampled_xyz} um in x y z. "
    f"downsampling to {goal_xyz} um. "
    f"dowmsample ratio is {ratioxyz}.")

    print(message)

    input("Press Enter to continue...")

    with open(resample_file,'r') as resampled_anno:
        resampled_anno_data=resampled_anno.readlines()
    # heading is stored in anno_data[2], 1st line basically useless

    headings=resampled_anno_data[2].rstrip('\n').replace(' ', '').split(',')
    resampled_annotations=[lines[0:-5].split(' ') for lines in resampled_anno_data[3:]]
    #slight modification on replacing and stripping due to the format of the resampled swc

    resampled_annotation_df=pd.DataFrame(resampled_annotations, columns=headings)

    ds_x= pd.to_numeric(resampled_annotation_df['x'])
    ds_x=ds_x/ratioxyz
    ds_xround=ds_x.astype(int).astype(str)
    ds_y= pd.to_numeric(resampled_annotation_df['y'])
    ds_y=ds_y/ratioxyz
    ds_yround=ds_y.astype(int).astype(str)

    ds_z= pd.to_numeric(resampled_annotation_df['z'])
    ds_z=ds_z/ratioxyz
    ds_zround=ds_z.astype(int).astype(str)

    ds_coordinates= pd.DataFrame(columns=['x','y','z'])
    ds_coordinates['x']=ds_xround
    ds_coordinates['y']=ds_yround
    ds_coordinates['z']=ds_zround

    q = [' '.join(x) for x in zip(ds_xround,ds_yround,ds_zround)]

    if name+'D' in resample_file:
        out_name= name+ f'D_{goal_xyz}voxel_trace_1umStepsize.txt'
    else:
        out_name= name+ f'_{goal_xyz}voxel_trace_1umStepsize.txt'
    print(out_name)

    num_row=len(resampled_annotation_df.index)
    f=open(outdir+'/'+out_name,'w+')
    f.write('point'+'\n')
    f.write(str(num_row)+'\n')
    for lines in q:
        f.write(lines+'\n')
    f.close()

    print('Downsampling annotation done.')

    print('Finding endings of the annotations...')
    list_parent= resampled_annotation_df['pid'].to_numpy()

    parent_index=np.argwhere(list_parent=='-1')
    ending_index=parent_index[:-1]+1
    ending_index=np.insert(ending_index, 0, 0)
    print(f'There are {len(ending_index)} endings')

    endings_df=resampled_annotation_df.iloc[ending_index]

    m=re.search('\D{2}[0-9]{3}[D]', resample_file)
    if m :
        out_name_endings= outdir[3:]+ f'D_endings.csv'
    else:
        out_name_endings= outdir[3:]+ f'_endings.csv'
    print(outdir+'/'+out_name_endings)

    np.savetxt(outdir+'/'+out_name_endings, ending_index, delimiter=",", fmt='%i')
    print('File saved in the output directory.')
    return

def get_pt_natlas(dspoint_name,outdir):
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
    return all_points, atlas_name

def make_pd(allpoints, ending_indices,dir):
    ''' 
    Reads in all points as well as ending indices and formulates a pd structure.
    Input: all downsampled points (in transformix compatible format), corresponding indicies of endings 
    ouputs: a pandas dataframe with anatomical regions and their corresponding total points count and ending points count, list of atlas ID for each point
    
    '''
        
    ending_indices = np.genfromtxt(ending_indices, delimiter=',', dtype='int')
            
    with open(allpoints,'r') as output:
        outputpoint= output.readlines()
    
    all_points=[]
    
    if 'outputpoints' in allpoints:
        files=os.listdir(dir)
        atlas=[i for i in files if re.search("atlas.+\.mhd",i)][0]
        atlas_name=os.path.join(dir, atlas)

        for lines in outputpoint:
            m=re.search("(?:OutputIndexFixed = \[ )([0-9]+ [0-9]+ [0-9]+)", lines).groups(0)
            this_line= str(m[0]).split(' ')
            mypoints= [int(stuff) for stuff in this_line]
            all_points.append(mypoints)
    else:
        atlas_name=dir+'/ara2sample_atlas/result.mhd'
        for lines in outputpoint[2:]:
            this_line= lines.split (' ')
            mypoints= [int(stuff) for stuff in this_line]
            all_points.append(mypoints)

    atlas= sitk.ReadImage(atlas_name)
        
    points_in_atlas=[int(atlas[i]) for i in all_points ]
    #find an ID for all points
    endings_in_atlas=[points_in_atlas[i] for i in ending_indices]
    # find the IDs associated with endings

    unique_id=set(points_in_atlas)

    our_regions=atlas_labels.loc[atlas_labels['id'].isin (unique_id)]

    id_withcounts=[]
    for i in unique_id:
        id_withcounts.append([i, points_in_atlas.count(i), endings_in_atlas.count(i)])

    new_df= pd.DataFrame(id_withcounts, columns=['id', 'Total_counts','Endings_counts'])
    our_regionWcounts=pd.merge(atlas_labels, new_df)
    return our_regionWcounts.sort_values(by=['Total_counts']), points_in_atlas

def check_points(points_in_atlas):
    '''Checks whether all your points' ID is within the atlas labels
    Input: matching ID of the points (this is the second output from analysis_tools.make_pd)
    '''
    id_inatlas=[]
    for x in atlas_labels['id']:
        intID = int(x)
        id_inatlas.append(intID)

    # need to format this first ourselves,otherwise problematic for 0 and very large numbers (idk why)    

    num_of_zeros = [i for i, x in enumerate(points_in_atlas) if x == 0]
    # find the indices for which carries an id =0
    
    unique_id=set(points_in_atlas)
    
    for id_inbrain in unique_id:
        if id_inbrain not in id_inatlas:
            if id_inbrain==0:
                print(f'There are {len(num_of_zeros)} points with ID= {id_inbrain}, this index is outside of the brain, consider possible suboptimal image registration')
            else: 
                print(id_inbrain,'this index does not exist in allen reference atlas, see https://github.com/ChristophKirst/ClearMap/issues/37')
            warnings.warn('Some points do not have corresponding labels')
    return 


def findID_origional(origional_points, points_in_atlas,dir,axon=True):
    '''
    Associate origional points with atlas ID and name, default is axon
    Inputs: origional onverted resampled points and matching atlas ID (this is the second output from analysis_tools.make_pd)
    Ouputs a pd dataframe with x, y, z, atlas ID, region name, colour
    '''
    #Load the origional points
    with open(origional_points,'r') as anno:
        anno_data=anno.readlines()
    # heading is stored in anno_data[2], 1st line basically useless

    headings=anno_data[2].rstrip('\n').replace(' ', '').split(',')
    annotations=[lines.rstrip('0 1\n').split(' ') for lines in anno_data[3:]]
    #slight modification on replacing and stripping due to the format of the resampled swc
    annotation_df=pd.DataFrame(annotations, columns=headings)

    points_with_id= pd.DataFrame (zip(annotation_df['x'],annotation_df['y'], annotation_df['z'],points_in_atlas ), columns=['x', 'y','z', 'atlasID'])
    if axon==True:
        out_name=dir+'/resamp_oripoints_withID.csv'
    else:
        out_name=dir+'/D_resamp_oripoints_withID.csv'
    points_with_id.to_csv (out_name, index = None, header=True) #Don't forget to add '.csv' at the end of the path
    
    #Creates colour for each region
    uniqueID=np.unique(points_with_id['atlasID'])
    colour= np.linspace(1,np.size(uniqueID)+1, num=np.size(uniqueID),dtype='int')
    colourdict=dict(zip(uniqueID,colour))
    
    # find region name based on ID
    namedict=dict(zip(atlas_labels['id'],atlas_labels['name']))
    points_with_id['name'] = points_with_id['atlasID'].map(namedict)
    points_with_id['colour'] = points_with_id['atlasID'].map(colourdict)
    return points_with_id

