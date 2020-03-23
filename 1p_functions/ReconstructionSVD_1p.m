function [dffmean,dffmean_es,dff,trial_info]=ReconstructionSVD_1p(U,SV,roiMasks,trial_info,opt)
% This function reconstruct df/f using SVD compressed data. Mostly used for
% analyzing 1p imaging data.
% Modified from ReconstructionSVD.m
%
% [dff_mean,dff_mean_es,dff,dff_img]=ReconstructionSVD_v2(U,SV,roi_vec,trial_info)
% [dff_mean,dff_mean_es,dff,dff_img]=ReconstructionSVD_v2(U,SV,roi_vec,trial_info,opt)
%
% Input
%     U: spatial component obtained by SVD of raw imaging data
%     SV: product of singular value and temporal component obtained by SVD of
%     raw imaging data (100 pre-inhalation frames (1sec) + 200 post-inhalation
%     frames (2sec))
%     trial_info: information about imaging session
%     roi_vec: binary mask of ROIs 2d matrix of (pixel x roi)
%     opt:parameters for reconstruction
%     opt.pre     :time in ms included before inhalation onset
%     opt.post    :time in ms included after inhalation onset
%     parmas.trig    :vector of times in ms from inhalation onset. Used if dff is triggered by timings other than inhalation onset
%     opt.base    :cell arrays containing of vectors representing timebins in
%     ms from inhalation onset. Used if baseline for calculating dff is
%     different from -20 to +20 ms of inhalation onset
%     opt.sniffFilterUsed: true or false. true if sniff filter is used to
%     record sniff. This will compensate for 10ms delay imposed by the filter ( default:true)
%
% Output:
%     dff_mean (roi x time x stimulus): trial averaged dff for each of odor
%     dff_mean_es(roi x time x stimulus): dff_mean subtracted by dff of empty
%     trials. dff_mean_es = dff_mean - repmat(dff_mean(:,:,1),1,1,size(dff_mean,3))
%     dff (roi x time x trial): single trial dff, use trial_info.stim_num to see
%     odors presented in each trial
%     dff_img: images of activated glomeruli at earlier timing from inahaltion
%     onset.
%     trial_info: updated trial_info, added trials_unread and dff_img
%
% Hirofumi Nakayama 2019

if ndims(U)==2
%     U(pixel x svals);
    img_size = sqrt(size(U,1));
    u=U;
elseif ndims(U)==3
%     U(pixel x pixel x svals)
    img_size = size(U,1);
    u=reshape(U,img_size^2,[]);
end

stim_num=trial_info.stim_num;
fpt=size(SV,1)/length(stim_num);



%roiMasks (pixels x roi) or (pixels x pixels x roi)
if ndims(roiMasks)==3
    roiMasks = reshape(roiMasks,[],size(roiMasks,3));
end
roiMasks=roiMasks./sum(roiMasks);
Ug=roiMasks'*u;

usfac=10;bin=1/usfac;

%Final dfg or dfgm have -20ms to +300ms in this processing (Modified
%021920, do not consider filter delay). It will be handled below.
pre=20;%time before inhalation onset included in dfg 20ms pre odor + 10ms offset due to amplifier
post=300;%time after inhalation onset included in dfg
trigOpt=false;
baseOpt=false;
sniffFilterUsed=true;
if exist('opt','var')
    if isfield(opt, 'pre')
        pre = opt.pre;
    end
    if isfield(opt, 'post')
        post = opt.post;
    end
    if isfield(opt, 'trig')
        %In case trig is given,
        %signals are aligned to trig(i) of each trial not inh onset
        trigOpt = true;
        trig = opt.trig;
    end
    if isfield(opt, 'base')
        baseOpt = true;
        base = opt.base;% assume cell arrays containing baseline period of each trial
    end
    if isfield(opt,'sniffFilterUsed')
        sniffFilterUsed = opt.sniffFilterUsed;
    end
end

if sniffFilterUsed
    %In case sniff filter is used during aquisition. Filter will impose 10ms
    %delay in recorded signal.
    sniffOffset = 10;
else
    sniffOffset = 0;
end

%Number of pre inhalation frames contained in data
pre_inh_frames=101;%original Trialframes were cropped to -100 to +200 frames from inhalation onset

%Detect trials in which LED was turned off in the middle
%Added 021320
meanF = reshape(mean(SV*Ug',2),fpt,[]);
meanF = meanF(100-pre/10:100+post/10,:);%Only care about time bins analuzed
trial_info.trials_unread(sum(meanF<200|meanF>5000)>0)=1; %Trials in which LED turned off or blink

ind_read=find(trial_info.trials_unread==0);
dff=[];%(glom x frame x tr)
for trial=1:length(ind_read)
    tr=ind_read(trial);
    time_first_frame=double((trial_info.inh_frames(tr)-pre_inh_frames-1)*10+5);%time of first frame of -100 to +200 frames
    %added 5ms so that frame intensity is that of middle time point of each frame
    
    %frme including inh_onset at trial t(r) =
    %fpt*(tr-1)+trial_info.inh_frames(tr)
    if trigOpt
        %In case trigger event other than first inhalation onset
        time_record=double(trial_info.lag_cam_inh(tr))+trig(tr)-sniffOffset+[-pre+1:post];
    else
        time_record=double(trial_info.lag_cam_inh(tr))-sniffOffset+[-pre+1:post];
    end
    range=(fpt*(tr-1)+1):fpt*tr;%cropped -100 to +200 frames from inh_onset
    
    %frame 10+1 (=100+1 in upsampled frame)in range corresponds to time_inh_onset_frame
    sv_up=interp1(SV(range,:),1:bin:length(range));
    
    time_upsample=[1:size(sv_up,1)]+time_first_frame-1;
    
    sv=sv_up(ismember(time_upsample,time_record),:);%time_record from -pre to +post ms (-20 - +200ms)
    
    %use -20 to +20ms from inhalation onset as baseline to calculate df/f
    if baseOpt
        %Subtract 10ms to compensate for delay imposed by sniff filter
        %base{tr} is ms time from inhalation onset
        base_period = base{tr} + pre;
    else
        %Assymetry is due to 10ms delay imposed by sniff filter
        base_period=(pre-20+1):(pre+20);
    end
    
    fg=sv*Ug';
    dff(:,:,tr)=((fg-mean(fg(base_period,:)))./mean(fg(base_period,:)))';%Use -20 to +20ms as baseline
    
    %Use 50 to 250ms to get good image of spatial patterns
    dfi(:,:,tr)=(reshape((mean(sv(71:271,:)*u'))',img_size,img_size)-reshape((mean(sv(1:41,:)*u'))',img_size,img_size))./(reshape((mean(sv(1:41,:)*u'))',256,256));
end
% trial_info.stim_num=trial_info.stim_num(ind_read);

num_cond=length(unique(stim_num));
dffmean=[];

for c=1:num_cond
    ind=find(stim_num==c&trial_info.trials_unread==0);
    dffmean(:,:,c)=mean(dff(:,:,ind),3);
    dff_img(:,:,c)=mean(dfi(:,:,ind),3);
    
end
trial_info.df_img=dff_img;
dffmean_es=dffmean-repmat(dffmean(:,:,1),1,1,num_cond);
