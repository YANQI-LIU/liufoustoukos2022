'''Contains functions working with points preprocessing'''

import Neuron_analysis as na
import pandas as pd
import numpy as np

import warnings
import re

import tkinter as tk
import tkinter.filedialog as fdialog
from tkinter import simpledialog

import SimpleITK as sitk
import skimage
from skimage import io

def convert_anno():
    xy = 0.6
    # 0.8 for 2pt samples, 0.6 for COLM images
    z = 5
    outdir = fdialog.askdirectory(title='Please select the output directory')
    anno_file=fdialog.askopenfile(initialdir=outdir, title='Select the eswc file containing the annotations').name 
    print('Converting pixels to ums....')

    anno=open(anno_file,'r')
    anno_data=anno.readlines()
    # heading is stored in anno_data[2], 1st line basically useless

    headings=anno_data[2].rstrip('\n').split(' ')
    annotations=[lines.rstrip('\n').split(' ') for lines in anno_data[3:]]

    annotation_df=pd.DataFrame(annotations, columns=headings)

    annotation_df['x']=pd.to_numeric(annotation_df['x'])*xy
    annotation_df['y']=pd.to_numeric(annotation_df['y'])*xy
    annotation_df['z']=pd.to_numeric(annotation_df['z'])*z

    tfile = open(anno_file+'converted.eswc', 'a')
    tfile.write(anno_data[0])
    tfile.write(anno_data[1])
    tfile.write(anno_data[2])
    tfile.write(annotation_df.to_string(header=False, index=False))
    tfile.close()

    print('Done!')
    return

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
    name=na.find_mousename(outdir)

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

#     print('Finding endings of the annotations...')
#     list_parent= resampled_annotation_df['pid'].to_numpy()

#     parent_index=np.argwhere(list_parent=='-1')
#     ending_index=parent_index[:-1]+1
#     ending_index=np.insert(ending_index, 0, 0)
#     print(f'There are {len(ending_index)} endings')

#     endings_df=resampled_annotation_df.iloc[ending_index]

#     m=re.search('\D{2}[0-9]{3}[D]', resample_file)
#     if m :
#         out_name_endings= outdir[3:]+ f'D_endings.csv'
#     else:
#         out_name_endings= outdir[3:]+ f'_endings.csv'
#     print(outdir+'/'+out_name_endings)

#     np.savetxt(outdir+'/'+out_name_endings, ending_index, delimiter=",", fmt='%i')
#     print('File saved in the output directory.')
    return

def refill_section(all_points,name):
    ''' Adds the cropped z points back to the begining
    input: all_points, name of cropped atlas
    (usually these are the output from common function get_pt_natlas)
    '''
    lead, trail= na.find_crop(name)
    all_points_full=[]
    for i in all_points:
        new_z=i[2]+lead
        all_points_full.append([i[0],i[1],new_z])
    return all_points_full


def make_pd(all_points, ending_indices,out_name, axon=1):
    ''' 
    For sample2ara only
    Takes in all points (as a list, usually the output of na.refill_section ) as well as ending indices (document name) and formulates a pd structure.
    Input: downsampled points (in transformix compatible format), name of corresponding indicies of endings 
    ouputs: a pandas dataframe with anatomical regions and their corresponding total points count and ending points count, list of atlas ID for each point
    
    '''
        
    ending_indices = np.genfromtxt(ending_indices, delimiter=',', dtype='int')
                
    atlas= sitk.ReadImage(na.fullatlas_name)
        
    points_in_atlas=[int(atlas[i]) for i in all_points ]
    #find an ID for all points
    endings_in_atlas=[points_in_atlas[i] for i in ending_indices]
    # find the IDs associated with endings

    unique_id=set(points_in_atlas)

    our_regions=na.atlas_labels.loc[na.atlas_labels['id'].isin (unique_id)]

    id_withcounts=[]
    for i in unique_id:
        id_withcounts.append([i, points_in_atlas.count(i), endings_in_atlas.count(i)])

    new_df= pd.DataFrame(id_withcounts, columns=['id', 'Total_counts','Endings_counts'])
    our_regionWcounts=pd.merge(our_regions, new_df)
    
    if axon==1:
        out_name=out_name + 'axons_region_with_counts.xls'
    else:
        out_name=out_name + 'dendrites_region_with_counts.xls'
    
    our_regionWcounts.to_excel(out_name,index=None,header=True)

    return our_regionWcounts.sort_values(by=['Total_counts']), points_in_atlas


def find_point_id(points, atlas_name):
    
    ''' 
    August 18 2021 update
    For ara2sample only
    Takes in all points (as a list, usually the output of na.refill_section or na.get_pt_natlas ) and formulates a pd structure.
    Input: downsampled points (in transformix compatible format), name of corresponding atlas
    ouputs:list of atlas ID for each point
    
    '''
    
    image= sitk.ReadImage(atlas_name)
    atlas =sitk.GetArrayFromImage(image)
        
    points_in_atlas=[int(atlas[i[2], i[1],i[0]]) for i in points]
    #find an ID for all points
    return points_in_atlas
    

def make_pd_ara2sample(points_in_atlas, atlas_labels, out_name, axon=1):
    ''' 
    August 18 2021 update
    For ara2sample only
    Takes in all points (as a list, usually the output of na.refill_section or na.get_pt_natlas ) and formulates a pd structure.
    Input: downsampled points (in transformix compatible format) and name of corresponding atlas
    ouputs: a pandas dataframe with anatomical regions and their corresponding total points count and ending points count
    
    '''
    unique_id=set(points_in_atlas)

    our_regions=atlas_labels.loc[atlas_labels['region_id'].isin (unique_id)]

    id_withcounts=[]
    for i in unique_id:
        id_withcounts.append([i, points_in_atlas.count(i)])

    new_df= pd.DataFrame(id_withcounts, columns=['region_id', 'Total_counts'])
    our_regionWcounts=pd.merge(our_regions, new_df)
    
    if axon==1:
        out_name=out_name + 'axons_region_with_counts.xls'
    else:
        out_name=out_name + 'dendrites_region_with_counts.xls'
    
    our_regionWcounts.to_excel(out_name,index=None,header=True)

    return our_regionWcounts.sort_values(by=['Total_counts'])


def check_points(points_in_atlas):
    '''Checks whether all your points' ID is within the atlas labels
    Input: matching ID of the points (this is the second output from na.make_pd)
    '''
    id_inatlas=[]
    for x in na.atlas_labels['id']:
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


def make_point_csv(all_points, points_in_atlas,out_name, axon=1):
    '''Makes a csv file with all the downsampled points and its associated id with a header
    example: [[x,y,z,atlasID],
              [83,75,223,981],...]'''
    dspoints_unpack=pd.DataFrame((all_points), columns=['x','y','z'])
    dspoints_with_id= pd.DataFrame(zip(dspoints_unpack['x'],dspoints_unpack['y'],dspoints_unpack['z'],points_in_atlas), columns=['x','y','z','atlas_ID'])
    
    
    if axon==1:
        out_name=out_name + '_axons.csv'
    else:
        out_name=out_name + '_dendrites.csv'
    
    dspoints_with_id.drop_duplicates().to_csv (out_name,index=None,header=True)
    return


def findID_origional(origional_points, points_in_atlas,out_name,axon=1):
    '''
    Associate origional points with atlas ID and name, default is axon
    Inputs: name of origional onverted resampled points and matching atlas ID (this is the second output from na.make_pd)
    saves a csv file with x, y, z, atlasID
    Ouputs a pd dataframe with x, y, z, atlas ID, region name, (colour?)
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
    if axon==1:
        out_name=out_name+'_oripoints_withID.csv'
    else:
        out_name=out_name+'D_oripoints_withID.csv'
    points_with_id.to_csv (out_name, index = None, header=True) #Don't forget to add '.csv' at the end of the path
    
    #Creates colour for each region
    #uniqueID=np.unique(points_with_id['atlasID'])
    #colour= np.linspace(1,np.size(uniqueID)+1, num=np.size(uniqueID),dtype='int')
    #colourdict=dict(zip(uniqueID,colour))
    
    # find region name based on ID
    #namedict=dict(zip(atlas_labels['id'],atlas_labels['name']))
    #points_with_id['name'] = points_with_id['atlasID'].map(namedict)
    #points_with_id['colour'] = points_with_id['atlasID'].map(colourdict)
    return points_with_id
    
    
    