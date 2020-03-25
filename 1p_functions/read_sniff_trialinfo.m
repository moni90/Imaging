function [Sniff,data]=read_sniff_trialinfo(h5_name,pre,post)
% Be careful about lost packets.
% Usually, duration of lost packet is 400ms. In some situation, duration
% calcualted from h5 file is smaller than 400ms but the 400ms packets need
% to be inserted to avoid overall shift in sniff signal
%
%This function do followings
% 1. read sniff
% 2. find real inh_onset
% 3. reshape data(remove trial0 and replace inh_onset)
%
% Hirofumi Nakayama 2017

if ~exist('pre','var')
    pre=1000;
end

if ~exist('post','var')
    post=2000;
end

try
%     h5 file exist in different location
    fullpath=pathFinderH5(h5_name);
catch
    fullpath = pwd;
end


data=h5read(fullpath,'/Trials');
data=removeTrial0(data);
%Need to check if this manipulation is ok for 2p data
inh_onset=double(data.inh_onset);
num_trial=length(inh_onset);

h_info=h5info(fullpath);
h_info=h_info.Groups;
Keys={h_info.Name};
Keys=Keys(2:end);%Remove trial-0
%%
record_onset=zeros(size(inh_onset));
% for i=1:numel(Keys)
for i=1:numel(Keys)
    %i
    data2=h5read(fullpath,strcat(Keys{i},'/Events'));
    
    packet_sent_time=data2.packet_sent_time;
    sniff_samples=data2.sniff_samples;    %array for the duration of each packet
    
    %record onset(i) represent the time of the first packet of sniff data
    %in trial (i). Those sniff data is stored in Trial(i)/Events/packets in
    %h5 file.
    record_onset(i)=packet_sent_time(1)-sniff_samples(1); %onset time for Sniff{trial(i)}
    
end

%%
%Inhalation timing within sniff traces of each trial
inh_onset_local=inh_onset-record_onset;

Sniff=zeros(num_trial,pre+post);
Sniff_time=zeros(num_trial,pre+post);

%check missing sniff packets
packet_sent_time=[];
sniff_samples=[];
trial_index = [];
trial_subind = [];
for i=1:num_trial
    %check missing sniff packet
    events = h5read(fullpath,strcat(Keys{i},'/Events'));
    packet_sent_time = [packet_sent_time;events.packet_sent_time];
    sniff_samples = [sniff_samples;events.sniff_samples];
    trial_index = [trial_index;i*ones(length(events.sniff_samples),1)];
    trial_subind = [trial_subind;[1:length(events.sniff_samples)]'];
end

lost_packet = find(diff(packet_sent_time)~=sniff_samples(2:end))+1;
sniff_all=[];
packet_length=[];
for i=1:num_trial
    st=zeros(1,pre+post);
    try
        sniffcell=h5read(fullpath,strcat(Keys{i},'/sniff'));
        if ~isempty(lost_packet)
            pos_lost=lost_packet(trial_index(lost_packet)==i);%position of lost packet in time from start
            
            %             subind=trial_subind(lost_packet);
            %             subind(trial_index(lost_packet)~=i)=[];
            if ~isempty(pos_lost)
                for j=pos_lost
                    lost_dur = (packet_sent_time(j)-packet_sent_time(j-1))-length(sniffcell{trial_subind(j)});
                    fprintf('trial %d, %d-th packet lost for %d ms\n',i,trial_subind(j),lost_dur)
                    if lost_dur==400
                        sniffcell{trial_subind(j)}=[single(zeros(lost_dur,1));sniffcell{trial_subind(j)}];
                    else
                        fprintf('Unusual lost_dur(%dms), 400ms packets were inserted instead\n', lost_dur)
                        sniffcell{trial_subind(j)}=[single(zeros(400,1));sniffcell{trial_subind(j)}];
                        data.unusual_lost_packet = true;
                    end
                end
                
            end
        end
        sniff0=cell2mat(sniffcell);
        packet_length(i) = length(sniff0); %Length of packets in each trial after resolving missing packets
        
        sniff_all=[sniff_all;sniff0];
        %pad with 10000 zeros before and after sniff
        sniff=[zeros(10000,1);sniff0;zeros(10000,1)];
        sniff_pos=[-9999:0,1:length(sniff),length(sniff)+1:length(sniff)+10000];
        
        sniff_range=inh_onset_local(i)-pre:inh_onset_local(i)+post-1;
        
        st=sniff(ismember(sniff_pos,sniff_range));
        
        Sniff(i,:)=st;
        
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    catch
        sprintf('Error in trial %d',i)
    end
    
end

%%
%Find real inhalation onset (=point of threshold crossing)
%Check sniff alignment to fvOnTime or inh_onset
inh_onset=data.inh_onset;
fvOnTime=data.fvOnTime;

inh_pos=findfirst(Sniff_time==double(inh_onset),2);
fv_pos=findfirst(Sniff_time==double(fvOnTime),2);

pn=Sniff>mean(Sniff(:));
% pn=Sniff>mean(Sniff,2);
% pn=Sniff>mode(Sniff(:));
filt=[ones(1,30),-1*ones(1,30)];
trig=zeros(size(Sniff));
for t=1:size(Sniff,2)-length(filt)
    for tr=1:size(Sniff,1)
        trig(tr,t+length(filt)/2)=pn(tr,t:t+length(filt)-1)*filt';
        trig(tr,1:fv_pos(tr)+10)=0;
    end
end
inh_pos2=findfirst(trig==30,2);
ind_no_inh = find(inh_pos2==0);
if ~isempty(ind_no_inh)
    %Perfect match to filt may not happen in noisy sniff data
    inh_pos2(ind_no_inh) = findfirst(trig(ind_no_inh,:)>20,2);
end
data.inh_onset_voyeur=data.inh_onset;
inh_onset=data.inh_onset+int32(inh_pos2-inh_pos);
data.inh_onset=inh_onset;

sniff_all = [sniff_all;zeros(10000,1)]; %To prevent index out of bounds at the last trial

inh_inSniff=inh_onset-record_onset(1);
for i=1:num_trial
    if inh_onset(i)>0
        Sniff2(i,:)=sniff_all(inh_inSniff(i)-pre:inh_inSniff(i)+post-1);
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    else
        Sniff2(i,:)   = zeros(1,pre+post);
        Sniff_time(i,:)=zeros(1,pre+post);
    end
end

%Check corrected sniff
figure('Name',sprintf('%s_readDataFiles_tsm_h5__read_sniff_trialinfo.m',h5_name));
subplot(211);plot(Sniff(:,901:1100)');line([100,100], ylim,'LineStyle','--');
subplot(212);plot(Sniff2(:,901:1100)');line([100,100], ylim,'LineStyle','--');
% subplot(211);plot(Sniff(22,1:3000)');line([1000,1000], ylim,'LineStyle','--');
% line([1,1]*inh_pos(22),ylim,'Color','r')
% line([1,1]*inh_pos2(22),ylim,'Color','g')
% line([1,1]*fv_pos(22),ylim,'Color','k')
%
% subplot(212);plot(Sniff2(22,1:3000)');line([1000,1000], ylim,'LineStyle','--');
% line([1,1]*inh_pos(22),ylim,'Color','r')
% line([1,1]*inh_pos2(22),ylim,'Color','g')
% line([1,1]*fv_pos(22),ylim,'Color','k')

Sniff=Sniff2;
a=1;


