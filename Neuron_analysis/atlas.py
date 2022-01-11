'''Contains functions that deals with both atlas and points'''


import Neuron_analysis as na
import SimpleITK as sitk

import skimage
from skimage import io

import numpy as np
import pandas as pd


def make_pd_from_cropped(allpoints, ending_indices,dir):
    ''' 
    Ara2sample use only
    Reads in all points as well as ending indices and formulates a pd structure.
    Input: name of downsampled points (in transformix compatible format), name of corresponding indicies of endings 
    ouputs: a pandas dataframe with anatomical regions and their corresponding total points count and ending points count, list of atlas ID for each point
    
    '''
    ending_indices = np.genfromtxt(ending_indices, delimiter=',', dtype='int')
            
    all_points, atlas_name= get_pt_natlas(all_points,dir)
    
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
    return print('Done checking!')


def make_tif(all_points, atlas_name, outname,axon=1):
    ''' Project downsampled points on to a tiff stack, useful for overlaping with brain or template (ie, in imageJ)
    input: downsampled points in a list containing x y z ordinates as int, directory containing it (this is also the output directory) and whether annotation is axon or not (default True)
    example: [[12, 13, 25],
             [13, 14, 25],...]
    
    output: a tiff stack with the same dimensions of the brain/template/atlas mhd files with downsampled points only
    each point has a value of the number of occurences (since downsampling combines multiple points as one)
    '''
        
    print('Starting to saving tif files..')
    
    atlas= sitk.ReadImage(atlas_name)
    svolume=np.zeros(atlas.GetSize())
    #columns, rows, planes
    
    zplanes=[]
    for i in all_points:
        zplanes.append( i[2])
    zplanes=np.unique(zplanes)
    temp=np.zeros(atlas.GetSize()[0:2])
    thepoints=np.asarray(all_points)

    for i in zplanes:
        index= thepoints[:,2]==i
        uindex,counts=np.unique(thepoints[index],return_counts=True, axis=0)
        for j, lines in enumerate(uindex):
            coord1,coord2=lines[0:2]
            temp[coord1][coord2]= counts[j]
        svolume[:,:,i]=temp #write this in 
        temp=np.zeros(atlas.GetSize()[0:2]) #reset the empty plane after each z
        
    
    coronal_planetmp= np.swapaxes(np.int16(svolume),0,2)
    #for some reason, if just save stuff as tiff, it will save x planes of yz view
    #here we shift the 3rd dimension with the first dimension to obtain xy view
    if axon==1:
        out_name=outname + '_axons.tif'
    else:
        out_name=outname + '_dendrites.tif'

    io.imsave(out_name,coronal_planetmp)
    return 


def make_tif_1(all_points, atlas_name, outname,axon=1):
    ''' Project downsampled points on to a tiff stack, useful for overlaping with brain or template (ie, in imageJ)
    
    August 18 2021
    Slightly modified version where each pixel will have a value of 1, instead of actual number of points.
    Useful for visualization (ie. max projection) so that all pixel axon or dendrite is equally strong in intensity
    
    
    input: downsampled points in a list containing x y z ordinates as int, directory containing it (this is also the output directory) and whether annotation is axon or not (default True)
    example: [[12, 13, 25],
             [13, 14, 25],...]
    
    output: a tiff stack with the same dimensions of the brain/template/atlas mhd files with downsampled points only
    each point has a value of the number of occurences (since downsampling combines multiple points as one)
    '''
        
    print('Starting to saving tif files..')
    
    atlas= sitk.ReadImage(atlas_name)
    svolume=np.zeros(atlas.GetSize())
    #columns, rows, planes
    
    zplanes=[]
    for i in all_points:
        zplanes.append( i[2])
    zplanes=np.unique(zplanes)
    temp=np.zeros(atlas.GetSize()[0:2])
    thepoints=np.asarray(all_points)

    for i in zplanes:
        index= thepoints[:,2]==i
        uindex,counts=np.unique(thepoints[index],return_counts=True, axis=0)
        for lines in uindex:
            coord1,coord2=lines[0:2]
            temp[coord1][coord2]= 1
        svolume[:,:,i]=temp #write this in 
        temp=np.zeros(atlas.GetSize()[0:2]) #reset the empty plane after each z
        
    
    coronal_planetmp= np.swapaxes(np.int16(svolume),0,2)
    #for some reason, if just save stuff as tiff, it will save x planes of yz view
    #here we shift the 3rd dimension with the first dimension to obtain xy view
    if axon==1:
        out_name=outname + '_axons.tif'
    else:
        out_name=outname + '_dendrites.tif'

    io.imsave(out_name,coronal_planetmp)
    return 