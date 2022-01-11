%% Path to code and data folder
code_path = 'D:\mmoroni\github_repos\Imaging';
working_dir = pwd;
addpath(genpath(code_path));
path_to_dff = 'D:\mmoroni\DATA\olfaction_metric\1P\data\Behavior_odors\Set1';
%Latency
file_name = 'Latency_data_8pures_16mixtures_022220.mat';

%% Load data
load(fullfile(path_to_dff, file_name));

%% Some plots of responses and latencies
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_glom_plot = 6;
n_odor_plot = 3;%size(this_dff,3)
for id_mouse = 1%:length(Dff)
%     dff_temp = Dff{1}; %#glom * time * #odors * concentration
    for conc_level = 1:size(Dff{id_mouse},4)
        this_dff = Dff{id_mouse}(:,:,:,conc_level);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0;
        figure; imagesc((this_lat-200)*frame_period); colorbar;
        xlabel('odor'); ylabel('glomeruli'); caxis([0, nanmax((this_lat(:)-200)*frame_period)]);
        if conc_level == 1
            title('low concentration');            
        else
            title('high concentration');
        end
        for id_odor = 1:n_odor_plot
            dff_this_odor = this_dff(:,:,id_odor);
            is_active = 1-isnan(this_lat(:,id_odor));
            figure; imagesc(time, [], dff_this_odor); colorbar;
            hold on; line([0, 0], [1,size(this_lat,1)], 'Color','k','lineWidth',2);
            hold on; scatter(time(this_lat(find(is_active),id_odor)),find(is_active),'r');
            xlabel('time'); ylabel('glomeruli'); title(odors{id_odor});
            figure;
            for id_glom = 1:n_glom_plot
                subplot(n_glom_plot,1,id_glom);
                plot(time,dff_this_odor(id_glom,:),'k');
                if is_active(id_glom)
                    hold on; stem(time(this_lat(id_glom, id_odor)), dff_this_odor(id_glom,this_lat(id_glom, id_odor)),'r');
                end
                xlim([time(1), time(end)]);
                title(['Glom. ' num2str(id_glom)]);
            end
       
        end
    end
    
end

%% Fraction of active glomeruli
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 25;
fract_active = zeros(length(Dff), 2, n_odors);
for id_mouse = 1:length(Dff) %iterate over mice
    for conc_level = 1:size(Dff{id_mouse},4) %iterate over conc level
        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
        n_glom = size(this_dff,1);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
        for id_odor = 1:n_odors
            is_active = 1-isnan(this_lat(:,id_odor)); %find active glomeruli
            fract_active(id_mouse, conc_level, id_odor) = mean(is_active);          
        end
    end
end

figure;
id_subplot = 1;
for conc_level = 1:size(Dff{id_mouse},4)
    subplot(1,2, id_subplot);
    id_subplot = id_subplot+1;
    hold on; bar(1:1:length(Dff), nanmean(squeeze(fract_active(:, conc_level, :)),2),...
        'FaceColor',[.7 .7 .7],'EdgeColor','k');
    hold on; errorbar(1:1:length(Dff), nanmean(squeeze(fract_active(:, conc_level, :)),2), ...
        nanstd(squeeze(fract_active(:, conc_level, :)),[],2)/sqrt(size(fract_active,3)),...
        'k','lineStyle','none','lineWidth',1);
    ylim([.6,1.02]); xlabel('mouse ID'); ylabel('fraction active glom')
    if conc_level == 1
        title('low concentration');
    else
        title('high concentration');
    end
end

delta_fract = squeeze(fract_active(:,2,:)-fract_active(:,1,:));
figure;
scatter(ones(size(delta_fract,1),1),nanmean(delta_fract,2),30,[.7 .7 .7],'filled');
hold on; errorbar(1,nanmean(nanmean(delta_fract,2)),...
    nanstd(nanmean(delta_fract,2))/sqrt(length(nanmean(delta_fract,2))),...
    'k','lineWidth',2);
hold on; scatter(1,nanmean(nanmean(delta_fract,2)),100,'k','s','filled')
set(gca,'XTick',[]); ylabel('Fract active high-low conc')

%% Average response strength glomeruli
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 25;
dff_max = zeros(length(Dff), 2, n_odors);
for id_mouse = 1:length(Dff) %iterate over mice
    for conc_level = 1:size(Dff{id_mouse},4) %iterate over conc level
        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
        n_glom = size(this_dff,1);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
        for id_odor = 1:n_odors
            is_active = 1-isnan(this_lat(:,id_odor)); %find active glomeruli          
            dff_max(id_mouse, conc_level, id_odor) = nanmean(nanmax(this_dff(find(is_active),:,id_odor),[],2),1);          
        end
    end
end

plot_var = dff_max;
figure;
id_subplot = 1;
for conc_level = 1:size(Dff{id_mouse},4)
    subplot(1,2, id_subplot);
    id_subplot = id_subplot+1;
    hold on; bar(1:1:length(Dff), nanmean(squeeze(plot_var(:, conc_level, :)),2),...
        'FaceColor',[.7 .7 .7],'EdgeColor','k');
    hold on; errorbar(1:1:length(Dff), nanmean(squeeze(plot_var(:, conc_level, :)),2), ...
        nanstd(squeeze(plot_var(:, conc_level, :)),[],2)/sqrt(size(plot_var,3)),...
        'k','lineStyle','none','lineWidth',1);
    ylim([0,0.16]);xlabel('mouse ID'); ylabel('dff max')
    if conc_level == 1
        title('low concentration');
    else
        title('high concentration');
    end
end

delta_plot = squeeze(plot_var(:,2,:)-plot_var(:,1,:));
figure;
scatter(ones(size(delta_plot,1),1),nanmean(delta_plot,2),30,[.7 .7 .7],'filled');
hold on; errorbar(1,nanmean(nanmean(delta_plot,2)),...
     nanstd(nanmean(delta_plot,2))/sqrt(length(nanmean(delta_plot,2))),...
     'k','lineWidth',2);
hold on; scatter(1,nanmean(nanmean(delta_plot,2)),100,'k','s','filled')
set(gca,'XTick',[]); ylabel('Max dff high-low conc')

%% Neural distance based on timing
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
similarity_index=zeros(id_mouse, conc_level, n_odors, n_odors);
for id_mouse = 1:length(Dff) %iterate over mice
    for conc_level = 1:size(Dff{id_mouse},4) %iterate over conc level
        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
        n_glom = size(this_dff,1);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
%         this_lat(isnan(this_lat)) = Inf;
        for id1_odor = 1:n_odors
            for id2_odor = 1:n_odors
                latency_diff = abs(this_lat(:,id1_odor+1)-this_lat(:,id2_odor+1));
                id_zero = find(isnan(this_lat(:,id1_odor+1)) + isnan(this_lat(:,id2_odor+1))==2);
                latency_diff(id_zero) = 0;
                latency_diff(isnan(latency_diff)) = Inf;
                similarity_index(id_mouse, conc_level,id1_odor,id2_odor) = 1-sum(exp(-latency_diff))/n_glom;
            end
        end
    end
end

figure;
for id_mouse = 1:7
    subplot(2,4,id_mouse); colormap(gray);
    imagesc(squeeze(similarity_index(id_mouse,1,:,:))); caxis([0, .8]);
    xlabel('odor'); ylabel('odor'); colorbar; 
    title(['mouse n. ', num2str(id_mouse)]); 
end
figure;
for id_mouse = 1:7
    subplot(2,4,id_mouse); colormap(gray);
    imagesc(squeeze(similarity_index(id_mouse,2,:,:))); caxis([0, .8]);
    xlabel('odor'); ylabel('odor'); colorbar;
    title(['mouse n. ', num2str(id_mouse)]);
end

figure; colormap(gray)
avg_sim = squeeze(nanmean(similarity_index(:,1,:,:),1));
imagesc(avg_sim); 
xlabel('odor'); ylabel('odor'); colorbar; %caxis([0, 1]);
figure; colormap(gray)
imagesc(squeeze(nanmean(similarity_index(:,2,:,:),1))); 
xlabel('odor'); ylabel('odor'); colorbar; %caxis([0,1]);

%% Neural distance based on rank
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
similarity_index=zeros(id_mouse, conc_level, n_odors, n_odors);
for id_mouse = 1:length(Dff) %iterate over mice
    for conc_level = 1:size(Dff{id_mouse},4) %iterate over conc level
        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
        n_glom = size(this_dff,1);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
        this_rank = LatRank{id_mouse}(:,:,conc_level);
        [~,this_rank_2] = sort(this_rank,1);
        similarity_index(id_mouse, conc_level,:,:) = 1-(0.5+corr(this_rank_2(:,2:end),'Type','Spearman')/2);
    end
end

figure;
for id_mouse = 1:7
    subplot(2,4,id_mouse); colormap(gray);
    imagesc(squeeze(similarity_index(id_mouse,1,:,:))); %caxis([0, .8]);
    xlabel('odor'); ylabel('odor'); colorbar; 
    title(['mouse n. ', num2str(id_mouse)]); 
end
figure;
for id_mouse = 1:7
    subplot(2,4,id_mouse); colormap(gray);
    imagesc(squeeze(similarity_index(id_mouse,2,:,:))); %caxis([0, .8]);
    xlabel('odor'); ylabel('odor'); colorbar;
    title(['mouse n. ', num2str(id_mouse)]);
end

figure; colormap(gray)
avg_sim = squeeze(nanmean(similarity_index(:,1,:,:),1));
imagesc(avg_sim); 
xlabel('odor'); ylabel('odor'); colorbar; %caxis([0, 1]);
figure; colormap(gray)
imagesc(squeeze(nanmean(similarity_index(:,2,:,:),1))); 
xlabel('odor'); ylabel('odor'); colorbar; %caxis([0,1]);

%% Neural distance based on STM
%model parameters
tau_glob = 100*(1e-1);
tau_act = 60*(1e-1);
tau_prim = 200*(1e-1);
teta = pi/2;
beta_glob = 1.1821;
beta_rel = 0.0263;

fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
distance_glob=zeros(id_mouse, conc_level, n_odors, n_odors);
distance_rel=zeros(id_mouse, conc_level, n_odors, n_odors);
% similarity_index=zeros(id_mouse, conc_level, n_odors, n_odors);
for id_mouse = 1:length(Dff) %iterate over mice
    for conc_level = 1:size(Dff{id_mouse},4) %iterate over conc level
        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
        n_glom = size(this_dff,1);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
%         this_lat(isnan(this_lat)) = Inf;
        for id1_odor = 1:n_odors
            for id2_odor = id1_odor:n_odors
                pattern_1 = this_lat(:,id1_odor+1);
                pattern_2 = this_lat(:,id2_odor+1);
                f_1 = from_discrete_to_waveform(pattern_1, tau_act, tau_prim, time);
                f_2 = from_discrete_to_waveform(pattern_2, tau_act, tau_prim, time);
                
                f_1_rot = apply_rotation(f_1, teta);
                CM_1 = compute_CM(f_1,time);
                f_2_rot = apply_rotation(f_2, teta);
                CM_2 = compute_CM(f_2,time);
                
                [dist_glob, dist_rel] = compute_dist(f_1_rot, f_2_rot, CM_1, CM_2, tau_glob, time);
                distance_glob(id_mouse, conc_level,id1_odor,id2_odor) = dist_glob;
                distance_rel(id_mouse, conc_level,id1_odor,id2_odor) = dist_rel;
                distance_glob(id_mouse, conc_level,id2_odor,id1_odor) = dist_glob;
                distance_rel(id_mouse, conc_level,id2_odor,id1_odor) = dist_rel;
%                 similarity_index(id_mouse, conc_level,id1_odor,id2_odor) = beta_glob*dist_glob + beta_rel*dist_rel;
            end
        end
    end
end
%%
% similarity_index = distance_glob;
% similarity_index = distance_rel;
similarity_index = beta_glob*distance_glob + beta_rel*distance_rel;

figure;
for id_mouse = 1:7
    subplot(2,4,id_mouse); colormap(gray);
    imagesc(squeeze(similarity_index(id_mouse,1,:,:))); %caxis([0, .8]);
    xlabel('odor'); ylabel('odor'); colorbar; 
    title(['mouse n. ', num2str(id_mouse)]); 
end
figure;
for id_mouse = 1:7
    subplot(2,4,id_mouse); colormap(gray);
    imagesc(squeeze(similarity_index(id_mouse,2,:,:))); %caxis([0, .8]);
    xlabel('odor'); ylabel('odor'); colorbar;
    title(['mouse n. ', num2str(id_mouse)]);
end

figure; colormap(gray)
avg_sim = squeeze(nanmean(similarity_index(:,1,:,:),1));
imagesc(avg_sim); 
xlabel('odor'); ylabel('odor'); colorbar; %caxis([0, 1]);
figure; colormap(gray)
imagesc(squeeze(nanmean(similarity_index(:,2,:,:),1))); 
xlabel('odor'); ylabel('odor'); colorbar; %caxis([0,1]);

