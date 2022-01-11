''' For ARA2sample: Transform the atlas using Transformix

User input: output directory, transformParameters, and atlas to be transformed

NOTE: This py file needs to be under the same directory as elastix to work
'''

import os

import subprocess

import tkinter.filedialog as fdialog


tempdir = fdialog.askdirectory(title='Please select a output directory (ALXXX/ara2sample_atlas or ara2sample_bf)')

transparam0=fdialog.askopenfile(initialdir=tempdir, title='select the transformparameter0 file').name
transparam1=fdialog.askopenfile(initialdir=tempdir, title='select the transformparameter1 file').name
atlasname=fdialog.askopenfile(initialdir=tempdir, title='select the corresponding atlas file').name

tparam1=open(transparam1,'r')
tparam1_data=tparam1.readlines()
tparam1.close()

new_transparam1_name=tempdir+'/tx_transparam1.txt'
new_tparam1=open(new_transparam1_name,'w+')

for lines in tparam1_data:
    if 'FinalBSplineInterpolationOrder' in lines:
        print(lines)
        print('Changing FinalBSplineInterpolationOrder 3 to 0')
        print(' ')
        new_line=lines.replace('3','0')
        new_tparam1.write(new_line)
    elif 'ResultImagePixelType ' in lines:
        print(lines)
        print('Changing ResultImagePixelType to int')
        new_line1=lines.replace('short', 'int')
        new_tparam1.write(new_line1)
    else:
        new_tparam1.write(lines)
new_tparam1.close()

command=['transformix', 
         '-out', tempdir,
         '-tp',new_transparam1_name, 
         '-in', atlasname
        ]

subprocess.run(command, cwd= 'C:/Users/liu')