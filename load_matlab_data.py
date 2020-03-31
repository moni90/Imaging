'''
Script to load .mat file in python
Preventing nested structures in the conversion of matlab structure to python dictionary
data = load_mat(filename=fname, folder_path=path)
Hirofumi Nakayama
6/3/2019
'''

import os
import pandas as pd
import numpy as np
import math
import random
import scipy.io
import h5py
import numbers
import h5py
import copy


def check_data_structure(d):
    # dt, dtnames, len_dt, types, shapes, pos_O = check_data_structure(d)
    tmp = d
    dt = []
    dtnames = []
    types = []
    shapes = []
    for i in range(10):
        try:
            dt.append(tmp.dtype)
            dtnames.append(tmp.dtype.names)
            types.append(type(tmp))
            shapes.append(np.shape(tmp))
            tmp = tmp[0]
        except:
            return dt, dtnames, [0 if x is None else len(x) for x in dtnames], types, shapes, [x.char is 'O' for x in dt]
    return dt, dtnames, [0 if x is None else len(x) for x in dtnames], types, shapes, [x.char is 'O' for x in dt]


def matlab_struct(d):
    # function to reformat matlab struct
    d2 = d[0][0]

    D = {}
    for k in d2.dtype.names:
        dt, dtnames, len_dt, types, shapes, pos_O = check_data_structure(d2[k])
        if type(d2[k][0]) is str or isinstance(d2[k][0], np.str_):
            D[k] = d2[k]
        else:
            # if d2[k][0][0].ndim == 1:
            #     #Not sure what is this for
            #     #Maybe for pure str variables in data
            #     D[k] = d2[k][0][0]
            # else:
            tmp = d2[k]
            if tmp.ndim == 1 or (tmp.ndim == 2 and np.shape(tmp)[0] == 1):
                D[k] = np.squeeze(d2[k][0])
            else:
                D[k] = np.squeeze(d2[k])

    return D


def matlab_cell_vals(d):
    # function to reformat matlab cell array
    # i-th element of matlab cell array D can be accessed by D[i][0]
    D = []
    for i in range(len(d)):
        d2 = d[i][0]

        D.append(d2)
    return D


# Todo: deal with struct.struct
def matlab_cell_struct(d):
    # function to reformat matlab cell array
    # i-th element of matlab cell array D can be accessed by D[i][0]
    D = []
    for i in range(len(d)):
        d2 = d[i][0]
        d3 = {}
        if type(d2[0]) is str or isinstance(d2[0], np.str_):
            D.append(d2[0])
        else:
            for k in d2.dtype.names:
                if d2[k][0][0].ndim == 1:
                    d3[k] = d2[k][0][0][0]
                else:
                    d3[k] = np.squeeze(d2[k][0][0])
            D.append(d3)
    return D


def load_mat(filename=None, folder_path=None):
    matdict = scipy.io.loadmat(os.path.join(folder_path, filename))

    # Todo
    # Try to read .mat file with -v7.3 format
    # matdict = {}
    # with h5py.File(os.path.join(folder_path,filename), 'r') as f:
    #     for k, data in zip(f.keys(),f.items()):
    #         matdict[k] = data

    data = {}

    for k in (kk for kk in matdict.keys() if '__' not in kk):
        print(k)
        if type(matdict[k]) is str or isinstance(matdict[k], np.str_):
            # exec('{KEY} = matdict["{KEY}"]'.format(KEY=k))
            data[k] = matdict[k]
        else:
            if len(matdict[k].dtype) > 0:
                # Matlab struct with multiple fields
                # exec('{KEY} = matlab_struct(matdict["{KEY}"])'.format(KEY=k))
                data[k] = matlab_struct(matdict[k])
            elif matdict[k].dtype == 'O':
                # Matlab cell array of struct
                # exec('{KEY} = matlab_cell_struct(matdict["{KEY}"])'.format(KEY=k))

                # Todo
                # choose one of them based on features in variables
                try:
                    data[k] = matlab_cell_struct(matdict[k])
                except:
                    data[k] = matlab_cell_vals(matdict[k])
            else:
                # exec('{KEY} = matdict["{KEY}"]'.format(KEY=k))
                data[k] = matdict[k]

    return data


def clean_TrialInfo(TrialInfo):
    #Clean up nested structure in TrialInfo and return dictionaries
    #TrialInfo is array containing trial_info of each sessions
    TrialInfo_dict = []
    for trial_info in TrialInfo:
        td={}
        for k in trial_info.dtype.names:
            td[k] = trial_info[0][0][k]
        TrialInfo_dict.append(td)
    return TrialInfo_dict