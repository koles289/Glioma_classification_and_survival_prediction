# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 14:25:41 2019

@author: Kristína Olešová
"""

class folder_content():
    def __init__(self,foldername=''):
        self._item=1
        self.foldername=foldername
        #self._item=[None].itemnumber
    def __len__(self):
        return self._data_len
    def __setitem__(self, key, data):
          self._item[key] = data
    def __getitem__(self, key):
          return self._item[key]

end_part=['_flair.nii.gz','_t1.nii.gz','_t1ce.nii.gz','_t2.nii.gz','_seg.nii.gz'] 
List_of_ROI_names=['result_HGG_I','result_LGG_I','result_HGG_II','result_LGG_II','result_HGG_III','result_LGG_III']
weight_seq=['FLAIR','T1','T1C','T2']
newpath='C:\\Users\\Kristína Olešová\\Diplomka_Jupyter_Files\\excel' 
main_path='D:\Data\Images'
pocet_vahovani=4;
ROI_list=['ROI_I','ROI_II','ROI_III']