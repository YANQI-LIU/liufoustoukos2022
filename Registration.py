'''Registration to Allen mouse ccf
Jan 10, 2022 modified just for anatomy MS in frontiers of neuroanatomy such that the directory organization is more automatic

Requirements: mhd formats of your image stack and allen brain template in the correct z orientation (can be done in imageJ to flip the Y order. Image->transform->flipz)

User inputs: define output directory, moving image, fixed image, elastix parameter file (in order)

NOTE: This py file needs to be under the same directory as elastix to work
'''

import os

import subprocess

import tkinter.filedialog as fdialog

# this is the gui for finding directory and files

import shutil

askdir = fdialog.askdirectory(title='Please select the folder with mouse name')
tempdir= os.path.join(askdir, 'ara2sample')
os.mkdir(tempdir)
#makes folder for transformed atlas

atlasdir= os.path.join(askdir, 'ara2sample_atlas')
bfdir= os.path.join(askdir, 'ara2sample_bf')
os.mkdir(atlasdir)
os.mkdir(bfdir)
#makes folder now for transforming atlas and bf mask in the future

fixed_img=fdialog.askopenfile(initialdir=askdir, title='select the fixed image').name
# need the .name part because askopenfile returns an io.textiowrapper not the full str name
moving_img=fdialog.askopenfile(initialdir=askdir, title='select the moving image').name

#param1=fdialog.askopenfile(title='select the first parameter file', initialdir =os.path.join(askdir, 'elastix_params') ).name
#param2=fdialog.askopenfile(title='select the second parameter file',initialdir = os.path.join(askdir, 'elastix_params')).name

paramsdir= os.path.join(askdir, 'elastix_params')
params= os.listdir(paramsdir)
param1= os.path.join(paramsdir, params[0])
param2= os.path.join(paramsdir, params[1])

shutil.copy(param1, tempdir)
shutil.copy(param2, tempdir)

command_line= ['elastix -f '+ fixed_img+ ' -m ' + moving_img+ 
               ' -out '+ tempdir+ ' -p '+ param1+
              ' -p ' + param2]
print('formulated commandline: ' , command_line)

input("Press Enter to continue...")

command=['elastix', 
         '-f', fixed_img, 
         '-m', moving_img,
         '-out', tempdir,
         '-p', param1, 
         '-p', param2
        ]

subprocess.run(command, cwd= 'C:/Users/liu')

print('The elastix command was: ', command_line)

answer=input('Rerun the previous elastix command?[y/n]')
while answer not in ('y', 'n'):
    print('please enter y or n, case sensitive')
    answer=input('123')
if answer == 'y':
    subprocess.run(command, cwd= 'C:/Users/liu')
else: 
    print('Finished')