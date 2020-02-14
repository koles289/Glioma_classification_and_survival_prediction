# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:40:35 2019

@author: Kristína Olešová
"""
import matplotlib
import numpy as np
import pandas as pd
import os
from PIL import Image  ### conda install pillow  
import matplotlib.pyplot as plt
from itertools import chain

import pymrmr
import IPython
from IPython.core.debugger import set_trace
from sklearn import preprocessing
import SimpleITK as sitk
from medpy.io import load
from radiomics import featureextractor
import re

import Dclr_of_class_variables as DCV





def Add_folder(path,folder):
    a=path+ '\\'+folder
    return a

def Create_class_structure(OWN_class,n):
    own_struct=[ OWN_class(n[i]) for i in range(len(n))]
    return own_struct;
    
def Create_path(BRATS,val,i,k,j):  
    return BRATS[val].foldername+'\\'+BRATS[val].type[i].foldername+'\\'+BRATS[val].type[i].grade[k].foldername+'\\'+BRATS[val].type[i].grade[k].foldername+DCV.end_part[j]

def Create_array(result,feature):
    a=[]
    for i in range(len(result)): 
        a.append(result[i].weightening[feature])
    return a

def Create_Image_repository():
    folder=os.listdir(DCV.main_path)
    BRATS =Create_class_structure(DCV.folder_content,folder)
    path_raw=Add_folder(DCV.main_path,folder[0])
    subfolder=os.listdir(path_raw)   #D:\data\Images\reálne data
    BRATS[0].type=Create_class_structure(DCV.folder_content,subfolder)
    BRATS[1].type=Create_class_structure(DCV.folder_content,subfolder)

    for k in range(len(folder)):
        for kk in range(len(subfolder)):
            path=Add_folder(DCV.main_path,folder[k])
            path_raw=Add_folder(path,subfolder[kk])
            subfolder4=os.listdir(path_raw) 
            BRATS[k].type[kk].grade=Create_class_structure(DCV.folder_content,subfolder4)
            subfolder3=os.listdir(path_raw)
            for i in range(len(subfolder3)):
                path=Add_folder(path_raw,subfolder3[i]) #D:\data\Images\reálne data\BRATS_HG0001............
                subfolder2=os.listdir(path)
                BRATS[k].type[kk].grade[i].image=Create_class_structure(DCV.folder_content,subfolder2)
                for j in range(len(subfolder2)):
                    BRATS[k].type[kk].grade[i].image[j]=subfolder2[j]
                    if k==1:
                        break
    return BRATS


def Create_ROI(Own_Structure) :
    for j in range(2):
        for i in range(len(Own_Structure[1].type[j].grade)):
            path_truth=Add_folder(DCV.main_path,Create_path(Own_Structure,1,j,i,4))
            image_data= sitk.ReadImage(path_truth)
            old=sitk.GetArrayFromImage(image_data)
            ROI_path=Add_folder(os.path.dirname(os.path.abspath(path_truth)),'ROI')
                
            ROI_I=old.copy();
            ROI_I[np.logical_or(ROI_I==1, ROI_I==3)] = 6;
            im=sitk.GetImageFromArray(ROI_I)
            im.SetDirection(image_data.GetDirection())
            im.SetOrigin(image_data.GetOrigin())
            im.SetSpacing(image_data.GetSpacing())
            
            pomm=os.path.basename(path_truth)
            pomm=pomm[:-10]
            
            sitk.WriteImage(im,Add_folder(ROI_path,pomm+'ROI_I.nii.gz'),True)
                
            ROI_II=ROI_I.copy(); 
            ROI_II[ROI_II == 4] = 6;
            im=sitk.GetImageFromArray(ROI_II)
            im.SetDirection(image_data.GetDirection())
            im.SetOrigin(image_data.GetOrigin())
            im.SetSpacing(image_data.GetSpacing())
            sitk.WriteImage(im,Add_folder(ROI_path,pomm+'ROI_II.nii.gz'),True)

            ROI_III=ROI_II.copy();
            ROI_III[ROI_III == 2] = 6;
            im=sitk.GetImageFromArray(ROI_III)
            im.SetDirection(image_data.GetDirection())
            im.SetOrigin(image_data.GetOrigin())
            im.SetSpacing(image_data.GetSpacing())
            sitk.WriteImage(im,Add_folder(ROI_path,pomm+'ROI_III.nii.gz'),True)
        
    return 0

def Create_ticks (df,indices):
	ticks=[]
	for i in range(5):
		ftr_name=re.sub('[^A-Z]', '', df.columns.values[indices[i]+1][1])
		if len(ftr_name)<2:
			pom_ftr=re.search(r'(?<=[A-Z])\w+',df.columns.values[indices[i]+1][1])
			ftr_name=ftr_name+pom_ftr.group(0)
		ticks.append(ftr_name+'\n'+ROI_I.columns.values[indices[i]+1][0])
	return ticks

def Saving_CSV(header,target,multi_dataset,ROI):
    multi_header = pd.MultiIndex.from_product([DCV.weight_seq,header],names=['weight_seq','features'])
    df = pd.DataFrame(multi_dataset,
                      index=target,
                      columns=multi_header)
    pom_path=DCV.newpath + '\\' + ROI +'.csv'
    df.to_csv(pom_path,sep=',')
    return multi_header

def Saving_normalized_CSV(HGG,LGG,ROI):
    header_pom=0
    for W_seq in (DCV.weight_seq):
        if header_pom==0:
            normalized_features,header=normalized_Clases(HGG,LGG,W_seq)
            target=[1]* len(HGG)
            target.extend([0]* len(LGG))
            header_pom=1
        else:
            normalized_features_step,_=normalized_Clases(HGG,LGG,W_seq)
            normalized_features=np.concatenate((normalized_features, normalized_features_step), axis=1)
    multi_header=Saving_CSV(header,target,normalized_features,ROI)
    return multi_header,target

def Radionomics_structure(extr,Own_Structure,file,glioma):
    grd=0;
    if glioma=='LGG':
        grd=1
    result = [ [[] for col in range(DCV.pocet_vahovani)] for col in range(len(Own_Structure[1].type[grd].grade))] 
    
    for i in range(len(Own_Structure[1].type[grd].grade)):
        path_truth=Add_folder(DCV.main_path,Create_path(Own_Structure,1,grd,i,4));
        ROI_pom=Add_folder(os.path.dirname(os.path.abspath(path_truth)),'ROI');
        pomm=os.path.basename(path_truth)
        ROI_path=Add_folder(ROI_pom,pomm[:-10]+file);
        for j in range(len(Own_Structure[0].type[grd].grade[i].image)):
            path_raw=Add_folder(DCV.main_path,Create_path(Own_Structure,0,grd,i,j))
            try:
                pom=extr.execute(path_raw,ROI_path)
            except ValueError:
                print (ROI_path)

            pom2=(list(pom.items()))
            del pom2[0:10]
            result[i][j]=list(pom.items())
            del pom,pom2
    return result


def normalized_Clases(features_HGG,features_LGG,wstr):
    poc_ind_radiopar=22;
    velikost=len(features_HGG)+len(features_LGG);
    features_matrix=np.zeros((velikost,len(features_HGG[0][0])-poc_ind_radiopar))
    #LGG=np.zeros((len(features_LGG[0][0])-poc_ind_radiopar,len(features_LGG)))
    wint=DCV.weight_seq.index(wstr)
    
    
    pocet_objektov=len(features_HGG[0][0]);
    if len(features_HGG[0][0])<len(features_LGG[0][0]):
        pocet_objektov=len(features_LGG[0][0]);
        
    
    header=[];
    
    i=poc_ind_radiopar;
    while i<pocet_objektov:
        header.append(features_HGG[0][1][i][0]);
        hg_pom=[];
        lg_pom=[];
        if i<=len(features_HGG[0][0]):
            hg_pom=np.array([li[wint][i][1] for li in features_HGG]);
        if i<=len(features_LGG[0][0]):    
            lg_pom=np.array([li[wint][i][1] for li in features_LGG]);
        
        x_pom=np.concatenate([hg_pom,lg_pom])
        features_matrix[:,i-poc_ind_radiopar] = preprocessing.normalize([x_pom])
        del x_pom
        i=i+1;
    return features_matrix,header

    
def Create_datasets(dt,Target_class,headers):
    
    dataset=pd.DataFrame(data=dt,    # values
            index=Target_class,  # 1st column as index
            columns=headers)
    
    return dataset
