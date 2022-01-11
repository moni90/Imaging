function data = read_behavioral_data(filename)
%this function reads the h5 file with behavioral data and save relevant
%quantities to compute the behavioral distance matrix

% filename = 'H:\Hiro data\Behavior\Odor_set1\3902\mouse3902_sess1_D2020_8_5T17_56_42.h5';
odors = {'Cinnamaldehyde_empty','Ethylbutyrate_empty','2MethylbutyricAcid_empty',...
    '22-DimethylbutyricAcid_empty','empty_CyclopentanecarboxylicAcid','empty_2-Heptanone',...
    'empty_IsobutyricAcid','empty_IsovalericAcid','Cinnamaldehyde_CyclopentanecarboxylicAcid',...
    'Cinnamaldehyde_2-Heptanone','Cinnamaldehyde_IsobutyricAcid',...
    'Cinnamaldehyde_IsovalericAcid','Ethylbutyrate_CyclopentanecarboxylicAcid',...
    'Ethylbutyrate_2-Heptanone','Ethylbutyrate_IsobutyricAcid',...
    'Ethylbutyrate_IsovalericAcid','2MethylbutyricAcid_CyclopentanecarboxylicAcid',...
    '2MethylbutyricAcid_2-Heptanone','2MethylbutyricAcid_IsobutyricAcid',...
    '2MethylbutyricAcid_IsovalericAcid','22-DimethylbutyricAcid_CyclopentanecarboxylicAcid',...
    '22-DimethylbutyricAcid_2-Heptanone','22-DimethylbutyricAcid_IsobutyricAcid',...
    '22-DimethylbutyricAcid_IsovalericAcid'};

info = h5info(filename);
data_list = {'/Events','/lick1','/lick2','/sniff'};
trials = h5read(filename,'/Trials');
for ii = 6:length(info.Groups)
    i = ii-5;
    info_temp = info.Groups(ii);
    for j = 1:length(data_list)
        dataset_name = strcat(info_temp.Name, data_list{j});
        if j==1
%             data(i).Events = h5read(filename, dataset_name);
        elseif j==2
            data(i).LickCenter = h5read(filename, dataset_name);
        elseif j==3
%             data(i).LickSide = h5read(filename, dataset_name);
        elseif j==4
%             data(i).Sniff = h5read(filename, dataset_name);
        end
    end
    data(i).Odor1_id = find(strcmp(odors,strcat(trials.olfa_0_2nd0x3Aodor(:,ii)','_',trials.olfa_1_2nd0x3Aodor(:,ii)')));
    data(i).Odor2_id = find(strcmp(odors,strcat(trials.olfas0x3Aolfa_00x3Aodor(:,ii)','_',trials.olfas0x3Aolfa_10x3Aodor(:,ii)')));
%     data(i).Odor1 = strcat(trials.olfa_0_2nd0x3Aodor(:,i)','_',trials.olfa_1_2nd0x3Aodor(:,i)');
%     data(i).Odor2 = strcat(trials.olfas0x3Aolfa_00x3Aodor(:,i)','_',trials.olfas0x3Aolfa_10x3Aodor(:,i)');
    
    response_onset = trials.starttrial(ii)+1500;
    go_trial = 0;
    if isempty(data(i).LickCenter)
        data(i).go_trial = go_trial;
    else
        lick_temp = data(i).LickCenter;
        for j = 1:length(lick_temp)
            this_lick = cell2mat(lick_temp(j));
            if (this_lick(1) >= response_onset) && (this_lick(1) <= (response_onset+1000))
                go_trial = 1;
                break;
            end
        end
        data(i).go_trial = go_trial;
    end
    
    
    
end
