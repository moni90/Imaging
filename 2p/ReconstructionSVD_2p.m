function [dff,trial_info,trials_read, meanAll]=ReconstructionSVD_2p(U,SV,roiMasks,trial_info,pre,post)
% This function reconstructs df/f using SVD compressed data. Mostly used
% for analyzing 2p imaging data
%
%Input
%U: spatial component obtained by SVD of raw imaging data
%SV: product of singular value and temporal component obtained by SVD of raw imaging data
%trial_info: information about imaging session
%roiMasks: binary mask of ROIs 2d(pixel x roi) or 3d matrix(pixel x pixel x roi)
%pre:time before inh_onset of each trial (ms)
%post:time after inh_onset of each trial (ms)
%
%Output
%dff: df/f (roi x time x trial)
%trial_info:structure containing trial information
%trials_read:0 or 1 to tell if the trials are correctly read
%meanAll: mean dff across all ROIs in each trial. Used for debugging
%
%Hirofumi Nakayama 2020

if ~exist('pre','var')
    %time before inh_onset of each trial (ms)
    pre=500;
end

if ~exist('post','var')
    %time after inh_onset of each trial (ms)
    post=1000;
end

inh_onset=trial_info.inh_onset;
frametrigger=trial_info.frametrigger;
num_trial=length(inh_onset);
%Temporal upsampling around inhalation onset
usfac=(1000/trial_info.fps);bin=1/usfac;

%roiMasks supposed to have shape of (pixel x roi)
%Reshape it if stack of binary images are given
if ndims(roiMasks)==3
   roiMasks = reshape(roiMasks,[],size(roiMasks,3)); 
end
Ucell=roiMasks'*U;

%number of frames pre- post- inhalation of each trial
pre_frame=40;
post_frame=150;
dff=[];
meanAll=[];
base = mean(SV)*Ucell';

for tr=1:num_trial
    inh=inh_onset(tr);
    inh_frame=find(frametrigger<inh, 1, 'last' );%frame in which inh_onset is included
    if ~isempty(inh_frame)&&inh_frame+post_frame-1<=size(SV,1)&&inh_frame-pre_frame>=1
        sv_range=inh_frame-pre_frame:inh_frame+post_frame-1;
        sv_upsample=interp1(SV(sv_range,:),1:bin:(pre_frame+post_frame));
        
        fcell=sv_upsample(single([inh-pre:inh+post-1])-frametrigger(inh_frame-pre_frame),:)*Ucell';
        
        %baseline subtraction using nonROI pixels
        %{
        f_nonRoi_centered_upsample = interp1(f_nonROI_centered(sv_range,:),1:bin:(pre_frame+post_frame))';
        fcell = fcell - f_nonRoi_centered_upsample(single([inh-pre:inh+post-1])-frametrigger(inh_frame-pre_frame),:)*weight;
        %}
        
        baseline=[-19:20]+pre; %for -20 to +20 from inhalation
        %dff(cellID, time from inh_onset(ms), trial)
        dff(:,:,tr)=((fcell-mean(fcell(baseline,:)))./mean(fcell(baseline,:)))';
        
        
        trials_read(tr)=true;
        meanAll(tr) = mean(mean(fcell(baseline,:)));
    else
        dff(:,:,tr)=zeros(size(Ucell,1),pre+post);
        ffprintf('trial %d inh_frame was not included in tiff stack',tr)
        trials_read(tr)=false;
        meanAll(tr) = 0;
    end
    
end

