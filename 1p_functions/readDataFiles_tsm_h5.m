function [Trialframes,trial_info]=readDataFiles_tsm_h5(tsm_name,h5_name,fps)
%Modified from Load_DataFiles_tsm_h5.m

%In case sniff warping is done, provide information about the first sniff
%inh_onset, lag_cam_inh is later modified with the timing information
%to calculate real inh_frames
%
%Hirofumi Nakayama 2017

%% set path
path_original=pwd;
try
    [direct]=configPath();
    try
        path_data = fullfile(direct.smdata,tsm_name);
    catch
        path_data = fullfile(direct.rinberg_smdata,tsm_name);
    end
catch
    path_data = pwd;
end

%% Obtain file names of image stacks

cd(path_data)
movFileList = dir(['*.tsm']); %List of files with .tsm extension
Names={movFileList.name};
cd(path_original)

%% Collect trial information (odor, inhalation timings etc)
%Cam-trigger start when totalms-enddtrial>iti
%final valve(FV) is turned on during exhalation after baseline_dur passed after camtrigger
%FV is kept on for 2 seconds after first inhalation onset
[Sniff,data]=read_sniff_trialinfo(h5_name);%This function gives 'data' with trial-0 being removed

%Number of trials analyzed
if isfield(data,'tr_imaged')
    t_len=min(numel(Names),nnz(data.tr_imaged));
else
    t_len=min(numel(Names),length(data.camtime));
end
%stim:information about odor identity of each trial
%Maybe necessary to modify this for taregt/Non-target discrimination
stim_desc=data.stim_desc;
stim=stim_desc(:,5)';
StimType={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c'};
stim_num=arrayfun(@(x) find(ismember(StimType,x)),stim);

camtime=data.camtime; %timing of the initiation of frame acqisition in each frame
inh_onset=data.inh_onset; %timing of inhalation onset in each trial
lag_cam_inh=inh_onset-data.camtime;

%inh_frames: frame number at the time of inhalation (fps + few frames)
inh_frames=round(lag_cam_inh*fps/1000);

num_cond=length(unique(stim));

first_frame_numbers=inh_frames-1*fps+1;

if isfield(data,'tr_imaged')
    %For Target Non-target discrimination
    %crop data for imaged trials only
    %Not sure, it is necessary to keep non-imaged trials
    tr_imaged=data.tr_imaged;tr_imaged=tr_imaged;
    tr_used=find(tr_imaged)';
else
    %For imaging only
    tr_used=1:t_len;
end
stim=stim(tr_used);
stim_num=stim_num(tr_used);
inh_onset=inh_onset(tr_used);
lag_cam_inh=lag_cam_inh(tr_used);
inh_frames=inh_frames(tr_used);
camtime=[camtime(tr_used);100000+camtime(tr_used(end))];%Shape of tr_used different between imaging and 2afc

od=data.odor_names;

field={'olfa0_air',...
    'olfa0_nit',...
    'olfa0_odor',...
    'olfa0_vialconc',...
    'olfa1_air',...
    'olfa1_nit',...
    'olfa1_odor',...
    'olfa1_vialconc',...
    'f_dilutor'};
if all(isfield(data,field))
    f_dilutor=data.f_dilutor;
    f_dilutor(f_dilutor<=10)=0;
    for i=1:length(inh_onset)
        OD{i}={deblank(data.olfa0_odor(i,:)),deblank(data.olfa1_odor(i,:))};
        c0=1e-6*data.olfa0_vialconc(i)*(single(data.olfa0_nit(i))/1000)*(1000-single(f_dilutor(i)))/1000;
        c1=1e-6*data.olfa1_vialconc(i)*(single(data.olfa1_nit(i))/1000)*(1000-single(f_dilutor(i)))/1000;
        
        if contains(OD{i}{1},'None')||contains(OD{i}{1},'empty')
            c0=0;
        end
        
        if contains(OD{i}{2},'None')||contains(OD{i}{2},'empty')
            c1=0;
        end
        conc{i,1}=[c0,c1];
    end
else
    conc=data.odorconc(2:end);
    for i=1:length(conc)
        s=deblank(od(i,:));
        if strfind(s,'None')
            s=erase(s,'None');
            s=erase(s,'_');
        else
            %For morphing stimuli
            s=strsplit(s,'_');
        end
        OD{i}=s;
    end
end

for i=1:max(stim_num)
    ind=find(stim_num==i,1);
    OdorNames{i}=OD{ind};
    if iscell(OD{ind})
        OdorNames_cat{i} = sprintf('%s_%s',OD{ind}{1},OD{ind}{2});
    else
        OdorNames_cat{i} = OD{ind};
    end
    try
        OdorConcs{i}=conc{ind};
    catch
        OdorConcs{i}=conc(ind);
    end
end
% other timing information for debugging
% starttrial = inh_onset
% no_sniff=data.no_sniff;no_sniff=no_sniff(2:end);%Remove trial 0
% lag_cam_fv=data.fvOnTime-data.paramsgottime;lag_cam_fv=lag_cam_fv(2:end);
% fvOn_frames=round(lag_cam_fv*fps/1000);fvOn_frames=fvOn_frames(2:end);

%% Read frames from individual trials
%In some trials, camera was not triggered due to poor sniff
%Identify trial without frame acquistion using camtime
cam_interval=camtime(2:end)-camtime(1:end-1);
framesCollected=find(cam_interval>0);
Ind_Frame2Trial=zeros(length(cam_interval),1);
Ind_Trial2Frame=zeros(length(cam_interval),1);
for tr=1:length(cam_interval)
    if tr<=length(framesCollected)
        Ind_Frame2Trial(tr)=framesCollected(tr);
    else
        Ind_Frame2Trial(tr)=0;
    end
    Ind_Trial2Frame(tr)=nnz(framesCollected<=tr)*(cam_interval(tr)>0);
end

%inf_frames is the array specifying the lag between camtime and inhalatin
%onset usually (fps + few) frames
%Trialframes=cell(t_len,1);
Trialframes=cell(t_len,1);

frame_dur=1000/fps;
odor_frames=double(data.fvdur(1)/frame_dur);
baseline_frames=double(data.baseline_dur(1)/frame_dur);
total_frames=double(2*baseline_frames+odor_frames);

ErrorData=(inh_frames-baseline_frames+1)<=0;
Truncated=(inh_frames+odor_frames)>total_frames;
TooLongSniff=inh_frames>total_frames-baseline_frames; %Trial at which odor period was truncated to below 1sec

frames_read=1:t_len;

frames_read=frames_read(~ErrorData&~Truncated);
trials_unread=zeros(1,t_len,'int16');
for i=intersect(frames_read,find(Ind_Trial2Frame~=0))'
    sprintf('i=%d,name=%s',i,Names{Ind_Trial2Frame(i)})
    Trialframes{i}=read_tsm_int16(Names{Ind_Trial2Frame(i)},inh_frames(i)-baseline_frames+1:inh_frames(i)+odor_frames);
end
trials_unread(intersect(frames_read,find(Ind_Trial2Frame==0)))=1;
trials_unread(inh_frames==0)=1;%Trials in which frames are taken but communication between arduino and voyeur failed

if nnz(Truncated)>0
    for i=find(Truncated)'
        if Ind_Trial2Frame(i)~=0
            Trialframes{i}=read_tsm_int16(Names{Ind_Trial2Frame(i)},[inh_frames(i)-baseline_frames+1:total_frames,total_frames*ones(1,inh_frames(i)+odor_frames-total_frames)]);
        else
            trials_unread(i)=2;
        end
    end
end
trials_unread(intersect(find(~Truncated)',find(Ind_Trial2Frame==0)))=2;

%This shouldn't happen but it did due to unknown bug in a behavior rig
if nnz(ErrorData)>0
    for i=intersect(find(ErrorData),find(Ind_Trial2Frame~=0))'
        tmp=read_tsm_int16(Names{Ind_Trial2Frame(i)},1:inh_frames(i)+odor_frames);
        Trialframes{i}(:,:,baseline_frames+odor_frames-size(tmp,3)+1:baseline_frames+odor_frames)=tmp;
        Trialframes{i}(:,:,1:baseline_frames+odor_frames-size(tmp,3))=repmat(mean(tmp(:,:,1:inh_frames(i)),3,'native'),[1,1,baseline_frames+odor_frames-size(tmp,3)]);
    end
end
trials_unread(intersect(find(~ErrorData)',find(Ind_Trial2Frame==0)))=3;

trials_unread(TooLongSniff)=4;

if nnz(trials_unread)>0
    frame_size=size(Trialframes{find(trials_unread==0,1)});
end

for i=find(trials_unread)
    Trialframes{i}=zeros(frame_size,'int16');
end

trial_info.lag_cam_inh=lag_cam_inh;
trial_info.inh_frames=inh_frames;
trial_info.stim=stim;
trial_info.stim_num=stim_num;
trial_info.num_cond=num_cond;
trial_info.first_frame_numbers=first_frame_numbers;
trial_info.Ind_Frame2Trial=Ind_Frame2Trial;
trial_info.Ind_Trial2Frame=Ind_Trial2Frame;
trial_info.ErrorData=ErrorData;
trial_info.Truncated=Truncated;
trial_info.trials_unread=trials_unread;
trial_info.data=data;
trial_info.tsm_name=tsm_name;
trial_info.h5_name=h5_name;
trial_info.fps=fps;
trial_info.inh_onset=inh_onset;
trial_info.OdorNames=OdorNames;
trial_info.OdorNames_cat=OdorNames_cat;
trial_info.OdorConcs=OdorConcs;
trial_info.Sniff=Sniff;