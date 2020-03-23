%demo_analysis_1p.m
%
%Hirofumi Nakayama 2020

%% Step1 - Read raw imaging data
%
%This part is for 1p imaging data taken at old 1p rig (RedShirt NeuroCCD)
%
SessNames{1}='1953_101419_100hz_8odors5acids_allMix_f10_149r4_';%Low conc (N2 flow = 10)
SessNames{2}='1953_101119_100hz_8odors5acids_allMix_149r4_';%High conc (N2 flow = 50)
SessNames{3}='1953_101719_100hz_8odors5acids_allMix_f50_149r4_';%High conc (N2 flow = 50)

H5Names{1} = '1953_101419-1.h5';
H5Names{2} = '1953_101119-3.h5';
H5Names{3} = '1953_101719-2.h5';

%parameters for reading raw data
path_raw_data='R:\Rinberglab\rinberglabspace\Users\Hiro\SMDATA';
path_h5='R:\Rinberglab\rinberglabspace\Users\Hiro\VoyeurData\HN10163';
fps=100;

for sess=1:numel(SessNames)
    % Read raw imaging data (binary file .tsm format)
    h5_name=H5Names{sess};
    path1 = fullfile(path_raw_data,SessNames{sess});
    path2 = fullfile(path_h5,H5Names{sess});
    [Trialframes,trial_info]=readDataFiles_tsm_h5(path1,path2,fps);
end

%% Step2 - Compress and denoise imaging data using SVD

path_trialframes = 'R:\Rinberglab\rinberglabspace\Users\Hiro\ImagingData_MATLAB';
%parameters for svd
[rr,cc] = meshgrid(1:256);
fov_mask = sqrt((rr-128).^2+(cc-128).^2)<=125;
num_svals = 200;
bwf = 0;
for sess=1:numel(SessNames)
    %load Trialframes and trial_info that are saved in step1
    load(fullfile(path_trialframes,strcat(SessNames{sess},'Trialframes_trial_info.mat')))
    
    % Compress and denoise imaging data using SVD
    [u,sv]=splitSVD(Trialframes,trial_info,fov_mask,num_svals,bwf,1); 
end

%% Step3 - Align frames from multiple sessions.
% Spatial components from SVD are aligned
% Before start this step, all sessions need to go through step1
% If single session data is analyzed, skip this step and go to Step4

fname = ''; %header that is added to data file when outputs are saved
[U, SV, TrialInfo, reg_params] = Registration_SVD(SessNames,fname,opt);

% Plot top spatial components of 1st session after registration
sess = 1;
figure('Position',get(0,'ScreenSize')*0.9)
for i=1:12
    subplot(3,4,i)
    imshow(mat2gray(U{sess}(:,:,i)));title(num2str(i))
end
%% Step4 - Calculate df/f using the result of SVD
% Create roi masks using manually segmented glomeruli in ImageJ
RoiSet_name = '1953_100hz_8odors5acids_allMix_149r4_RoiSet';
opt.path = 'R:\Rinberglab\rinberglabspace\Users\Hiro\1P\Segmentation\ManualSegmentation';
opt.img_size = 256;
[roiMask_id,roiMask_stack]=create_ROImask_manual(RoiSet_name,opt);

% In case running this code without step 1 - 2
path_data='R:\Rinberglab\rinberglabspace\Users\Hiro\ImagingData_MATLAB';
load(fullfile(path_data,'1953_8pure_16mixtures_SVD_registration.mat'))

% Calculalte df/f
clear opt
opt.pre = 20;
opt.post = 300;
sess = 1;
[dffmean,dffmean_es,dff,trial_info]=ReconstructionSVD_1p(U{sess},SV{sess},roiMask_stack,TrialInfo{sess},opt);

%Plot trial averaged df/f
%Each panel shows df/f of all roi in a given odor
figure('Position',get(0,'ScreenSize')*0.9)
for i=1:29
    subplot(5,6,i)
    plot(dffmean(:,:,i)');title(trial_info.OdorNames_cat{i})
end
