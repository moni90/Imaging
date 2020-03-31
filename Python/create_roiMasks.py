# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 15:55:09 2019

@author: Hirofumi
"""

#Code to convert manually segmented masks into glomerulus masks that can be used in matlab
#RoiSet.zip created by manual segmentation of 
#603_010419_100hz_4MValericA_2e10_2e12_144r4_dff_mean_adjInd.tiff

'''
/*  ImageJ/NIH Image 64 byte ROI outline header
    2 byte numbers are big-endian signed shorts
    
    0-3     "Iout"
    4-5     version (>=217)
    6-7     roi type (encoded as one byte)
    8-9     top
    10-11   left
    12-13   bottom
    14-15   right
    16-17   NCoordinates
    18-33   x1,y1,x2,y2 (straight line) | x,y,width,height (double rect) | size (npoints)
    34-35   stroke width (v1.43i or later)
    36-39   ShapeRoi size (type must be 1 if this value>0)
    40-43   stroke color (v1.43i or later)
    44-47   fill color (v1.43i or later)
    48-49   subtype (v1.43k or later)
    50-51   options (v1.43k or later)
    52-52   arrow style or aspect ratio (v1.43p or later)
    53-53   arrow head size (v1.43p or later)
    54-55   rounded rect arc size (v1.43p or later)
    56-59   position
    60-63   header2 offset
    64-       x-coordinates (short), followed by y-coordinates
*/
'''
import sys
from pathlib import Path
from pylab import rcParams
rcParams['figure.figsize'] = 10, 10

#Load RoiSet.zip created in ImageJ

import pickle as p
import cv2
import os
import numpy as np
import matplotlib.pyplot as plt
import zipfile

def read_roi(fileobj):
    #See https://imagej.nih.gov/ij/developer/source/ij/io/RoiDecoder.java.html for .roi file
    #decode byte literals to unicode (Necessary in python3 but not in python2)
    if fileobj[:4].decode() != 'Iout':
        raise IOError('Magic number not found {0}'.format(fileobj[:4]))

    y1 = ord(fileobj[9:10])
    x1 = ord(fileobj[11:12])
    y2 = ord(fileobj[13:14])
    x2 = ord(fileobj[15:16])

    frame = np.zeros(256**2)
    frame = frame.reshape(256, 256)
    frame[y1:y2, x1:x2] = 255
    return frame

def create_roiMasks(fname=None, path=r'C:\Users\hnaka\Dropbox\MATLAB\OMP_Gcamp\Segmentation\ManualSegmentation'):
    if '.zip' not in fname:
        fname=fname+'.zip'

    zf = zipfile.ZipFile(os.path.join(path, fname))
    roiMasks = []

    for roi in zf.namelist():
        rfile = zf.open(roi, 'r')
        fileobj = rfile.read()
        frame = read_roi(fileobj)
        roiMasks.append(frame)
    roiMasks = np.transpose(np.array(roiMasks), (1, 2, 0))
    return roiMasks
