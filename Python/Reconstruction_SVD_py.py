'''
Python implementation of Reconstruction_SVD.m
Hirofumi Nakayama
8/20/2019
'''


import sys
path_parent = r'C:\Users\hnaka\Dropbox\Python_codes\BehaviorRig'
sys.path.append(path_parent)

import os
import numpy as np

from load_matlab_data import load_mat

import matplotlib.pyplot as plt

class SVD_data():
    def __init__(self,matdict):
        self.U = matdict['U']
        self.SV = matdict['SV']
        self.TrialInfo = matdict['TrialInfo']
        self.SM = matdict['SM']
        self.numSess = len(self.U)

    def calcDff(self,sess=0, roiMasks=None, pre=30, post=990):
        # np.shape(roiMask) = (256^2, #_of_ROIs) or (256,256, #_of_ROIs)
        # Calculate dff of session (sess)
        # dff (roi,t,odors)

        u = self.U[sess]
        sv = self.SV[sess]
        trial_info = self.TrialInfo[sess]

        if np.ndim(roiMasks)==3:
            roiMasks = np.reshape(roiMasks, (256**2,-1))

        u_roi = np.matmul(np.transpose(roiMasks),np.reshape(u,(256**2,-1)))

        stim_num = trial_info['stim_num']
        trials_unread = trial_info['trials_unread']
        num_roi = np.shape(roiMasks)[-1]

        fpt = np.shape(sv)[0] / len(stim_num)
        dff_tr = np.zeros((num_roi, pre+post, len(trials_unread)))
        tr_index = np.where(trials_unread == 0)[0]
        for tr in tr_index:
            # Time of the first frame in -100 to +200 frames image stacks
            # added 5ms so that frame intensity becomes the intensity of middle time point of each frame
            time_record = np.arange(1000 - pre, 1000 + post)  # different from matlab version due to 0-indexing

            frame_mask = np.arange(fpt * tr, fpt * (tr + 1)).astype('int')  # different from matlab version due to 0-indexing

            sv_tr = sv[frame_mask, :]
            sv_up = []  # upsampled sv, 3000 time bins (from -1000ms to +2000ms)
            for s in range(np.shape(sv)[1]):
                sv_up.append(np.interp(np.arange(0, len(frame_mask), 0.1), np.arange(len(frame_mask)), sv_tr[:, s]))
            sv_up = np.transpose(np.vstack(sv_up))
            sv_tr = sv_up[time_record, :]  # upsampled sv from -20ms to 1000ms

            #Assuming that pre=20
            dff_base = np.transpose(np.tile(np.mean(np.matmul(u_roi, np.transpose(sv_tr[1:41, :])),axis=1),(np.shape(dff_tr)[1],1)))

            # a(256 x 256 x 520) - b(256 x 256) caused broadcasting error
            # Need to permute axis of a to (520,256,256) to avoid that
            dff_tr[:,:,tr] = (np.matmul(u_roi, np.transpose(sv_tr)) - dff_base) / dff_base

        num_cond = len(set(stim_num))
        dff = np.zeros((np.shape(dff_tr)[0], pre+post, num_cond))
        for c in range(num_cond):
            ind = np.logical_and(stim_num == (c + 1), trials_unread == 0)
            dff[:, :, c] = np.mean(dff_tr[:, :, ind], axis=2)
        return dff

    def calcDff_summaryImages(self, sess=0, time_points=[[50,150],[200,500],[600,1000]]):

        #Calculate dff of session (sess)
        #dff (x,y,t,odors)
        pre = 30
        post = 990
        u = self.U[sess]
        sv = self. SV[sess]
        trial_info = self.TrialInfo[sess]

        stim_num = trial_info['stim_num']
        trials_unread = trial_info['trials_unread']

        fpt = np.shape(sv)[0] / len(stim_num)
        dff_tr=np.zeros((256,256,len(time_points),len(trials_unread)))
        tr_index = np.where(trials_unread==0)[0]
        for tr in tr_index:
            #Time of the first frame in -100 to +200 frames image stacks
            #added 5ms so that frame intensity becomes the intensity of middle time point of each frame
            time_record = np.arange(1000 - pre, 1000 + post)  # different from matlab version due to 0-indexing

            frame_mask = np.arange(fpt*tr,fpt*(tr+1)).astype('int') #different from matlab version due to 0-indexing

            sv_tr=sv[frame_mask,:]
            sv_up=[] #upsampled sv, 3000 time bins (from -1000ms to +2000ms)
            for s in range(np.shape(sv)[1]):
                sv_up.append(np.interp(np.arange(0,len(frame_mask),0.1),np.arange(len(frame_mask)),sv_tr[:,s]))
            sv_up = np.transpose(np.vstack(sv_up))
            sv_tr = sv_up[time_record,:] #upsampled sv from -20ms to 1000ms

            im_base = self.usv2img_mean(u,sv_tr[1:41,:])

            #a(256 x 256 x 520) - b(256 x 256) caused broadcasting error
            #Need to permute axis of a to (520,256,256) to avoid that
            offset=20 #sv_tr contains 20ms prior to inh onset, need to add this offset for slicing
            for i,t in enumerate(time_points):
                print('trial{0}, time_step{1}'.format(tr,i))
                sv_average = [np.mean(sv_tr, axis=0)]
                dff_tr[:, :,i,tr]=(self.usv2img_mean(u,sv_average)-im_base)/im_base

        num_cond = len(set(stim_num))

        dff=np.zeros((256,256,len(time_points),num_cond))
        for c in range(num_cond):
            ind = np.logical_and(stim_num == (c+1), trials_unread==0)
            dff[:,:,:,c]=np.mean(dff_tr[:,:, :, ind], axis=3)

        return dff

    def usv2img_mean(self, u,sv):
        #sv = (time x svals)
        #u = (x-pix x y-pix x svals)
        #img = (x-pix x y-pix), averaged image across time
        img = np.reshape(np.mean(np.matmul(sv,np.transpose(np.reshape(u,(256**2,-1)))), axis=0),(256,256))
        return img

    def usv2img_stack(self, u,sv):
        #sv = (time x svals)
        #u = (x-pix x y-pix x svals)
        #stack = (x-pix x y-pix x time)
        stack = np.reshape(np.transpose(np.matmul(sv,np.transpose(np.reshape(u,(256**2,-1))))),(256,256,-1))
        return stack

if __name__ == '__main__':
    filename = '1953_11sess_081719__SVD_registration_not_v73.mat'
    folder_path=r'C:\Users\hnaka\Dropbox\MATLAB\OMP_Gcamp\LocalData'
    matdict=load_mat(filename=filename, folder_path=folder_path)


    svd_rec = SVD_data(matdict)

    time_points = [[50, 150], [200, 500], [600, 1000]]
    #dff (x-pixel, y-pixel, time, odor)
    dff = svd_rec.calcDff_summaryImages(sess=6,time_points=time_points)

    OdorNames = matdict['TrialInfo'][0]['OdorNames']
    fig, ax = plt.subplots(3,4,figsize=(13,6))
    for i in range(3):
        for j in range(4):
            im = ax[i][j].imshow(dff[:,:,0,j+4*i])
            plt.colorbar(im,ax=ax[i][j])
    a=1
