%% Path to code and data folder
code_path = 'D:\mmoroni\github_repos\Imaging';
working_dir = pwd;
addpath(genpath(code_path));
path_to_dff = 'D:\mmoroni\DATA\olfaction_metric\1P\data\Behavior_odors\Set1';
%Latency
file_name = 'Latency_data_8pures_16mixtures_022220.mat';
behavioral_dist_path = 'D:\mmoroni\Hiro_project\analyses\Behavior\Odor_set1\pooled_distance_matrix.mat';
%% Load data
load(fullfile(path_to_dff, file_name));

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
n_mouse = length(Dff);
n_conc_level = size(Dff{1},4);
% n_conc_level = 1;
distance_glob=zeros(n_mouse, n_conc_level, n_odors, n_odors);
distance_rel=zeros(n_mouse, n_conc_level, n_odors, n_odors);
% similarity_index=zeros(id_mouse, conc_level, n_odors, n_odors);
for id_mouse = 1:length(Dff) %iterate over mice
    for conc_level = 2%:size(Dff{id_mouse},4) %iterate over conc level (1=low, 2=high)
        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
        n_glom = size(this_dff,1);
        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
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
                distance_glob(id_mouse, conc_level, id1_odor, id2_odor) = dist_glob;
                distance_rel(id_mouse, conc_level, id1_odor, id2_odor) = dist_rel;
                distance_glob(id_mouse, conc_level, id2_odor, id1_odor) = dist_glob;
                distance_rel(id_mouse, conc_level, id2_odor, id1_odor) = dist_rel;
%                 similarity_index(id_mouse, conc_level,id1_odor,id2_odor) = beta_glob*dist_glob + beta_rel*dist_rel;
            end
        end
    end
end
%%
% similarity_index = distance_glob;
similarity_index = distance_rel;
% similarity_index = beta_glob*distance_glob + beta_rel*distance_rel;

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

%%
behavioral_data = load(behavioral_dist_path);
behavioral_distance_mat = behavioral_data.distance_mat; %behavioral_data.distance_mat_avg;
BD_lower = behavioral_distance_mat + tril(NaN*ones(n_odors,n_odors),-1);
ND_lower = avg_sim + tril(NaN*ones(n_odors,n_odors),-1);
trials_id = find(1-isnan(BD_lower));
[trials_odor1, trials_odor2] = ind2sub([n_odors, n_odors], trials_id);
trials_pair = [trials_odor1(:), trials_odor2(:)];
BD_lower(isnan(BD_lower)) = [];
ND_lower(isnan(ND_lower)) = [];
corr(BD_lower',ND_lower')
%%
n_trials = n_odors*(n_odors+1)/2;
train_fract = 0.8;
n_train = floor(train_fract*n_trials);
n_folds = 1;
trials_per_fold = floor(n_train/n_folds);
n_train = n_folds*trials_per_fold;
train_fract = n_train/n_trials;

test_fract = 1-train_fract;
n_test = n_trials - n_train;

rng(1)
train_id = randperm(n_trials, n_train);
test_id = setdiff(1:1:n_trials,train_id);
folds_id = [];
for i = 1:n_folds
    folds_id = [folds_id; i*ones(trials_per_fold,1)];
end
folds_id = folds_id(randperm(length(folds_id)));

%%
%model parameters
% tau_glob = 100*(1e-1);
% tau_act = 60*(1e-1);
% tau_prim = 200*(1e-1);
% teta = pi/2;
% beta_glob = 1.1821;
% beta_rel = 0.0263;

tau_glob = [100 200]*(1e-1);
tau_act = [40 80]*(1e-1);
tau_prim = [100 200]*(1e-1);
teta_angle = pi/2;
beta_glob = [0 0.01];% 1.2];
beta_rel = [0.01 0.1 1];% 0.1];

[t_g, t_a, t_p, teta, b_g, b_r] = ndgrid(tau_glob, tau_act, tau_prim, teta_angle, beta_glob, beta_rel);
t_g = t_g(:);t_a = t_a(:); t_p = t_p(:); teta = teta(:); b_g = b_g(:); b_r = b_r(:);



%% SINGLE TRIAL NEURAL RESPONSE FIT
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
n_mouse = length(Dff);
% n_conc_level = size(Dff{1},4);
conc_level_def = [2]; %(1=low, 2=high)
n_conc_level = length(conc_level_def);

% similarity_index=zeros(id_mouse, conc_level, n_odors, n_odors);
corr_train = zeros(length(t_g),n_folds);
corr_val = zeros(n_folds,1);
t_g_val = zeros(n_folds,1); t_a_val = zeros(n_folds,1); t_p_val = zeros(n_folds,1);
teta_val = zeros(n_folds,1); b_g_val = zeros(n_folds,1); b_r_val = zeros(n_folds,1);
if n_folds > 1
    for id_fold = 1:n_folds
        train_trials = train_id(folds_id~=id_fold);
        val_trials = train_id(folds_id==id_fold);
        
        for id_trial = 1:length(train_trials)
            id1_odor = trials_pair(train_trials(id_trial),1);
            id2_odor = trials_pair(train_trials(id_trial),2);
            distance_behav_train(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
        end
        
        for id_param = 1:length(t_g) %%%cycle over parameters to optimize
            disp(['Fold ' num2str(id_fold) ' of ' num2str(n_folds) '. Params: ' num2str(id_param) '/' num2str(length(t_g))]);
            t_g_temp = t_g(id_param); t_a_temp = t_a(id_param);
            t_p_temp = t_p(id_param); teta_temp = teta(id_param);
            b_g_temp = b_g(id_param); b_r_temp = b_r(id_param);
            
            distance_glob_temp=zeros(n_mouse, n_conc_level, length(train_trials));
            distance_rel_temp=zeros(n_mouse, n_conc_level, length(train_trials));
            distance_stm_temp=zeros(n_mouse, n_conc_level, length(train_trials));
            for id_trial = 1:length(train_trials)
                id1_odor = trials_pair(train_trials(id_trial),1);
                id2_odor = trials_pair(train_trials(id_trial),2);
                
                for id_mouse = 1:length(Dff) %iterate over mice
                    for conc_level_id = 1:n_conc_level %iterate over conc level
                        conc_level = conc_level_def(conc_level_id);
                        this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
                        n_glom = size(this_dff,1);
                        this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
                        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
                        
                        pattern_1 = this_lat(:,id1_odor+1);
                        pattern_2 = this_lat(:,id2_odor+1);
                        
                        [dist_glob, dist_rel, dist_stm] = compute_neural_distance_single_pair(pattern_1, pattern_2,...
                            t_a_temp, t_p_temp, t_g_temp, teta_temp, b_g_temp, b_r_temp, time);
                        
                        distance_glob_temp(id_mouse, conc_level_id, id_trial) = dist_glob;
                        distance_rel_temp(id_mouse, conc_level_id, id_trial) = dist_rel;
                        distance_stm_temp(id_mouse, conc_level_id, id_trial) = dist_stm;
                    end
                end
            end
            distance_glob = squeeze(nanmean(nanmean(distance_glob_temp,2),1));
            distance_rel = squeeze(nanmean(nanmean(distance_rel_temp,2),1));
            distance_stm = squeeze(nanmean(nanmean(distance_stm_temp,2),1));
            
            corr_train(id_param, id_fold) = corr(distance_stm(:), distance_behav_train(:));
            if id_param == 1
                t_g_val(id_fold) = t_g_temp; t_a_val(id_fold) = t_a_temp; t_p_val(id_fold) = t_p_temp;
                teta_val(id_fold) = teta_temp; b_g_val(id_fold) = b_g_temp; b_r_val(id_fold) = b_r_temp;
                max_corr = corr_train(id_param, id_fold);
                ind_max = id_param;
            elseif corr_train(id_param, id_fold)>max_corr
                t_g_val(id_fold) = t_g_temp; t_a_val(id_fold) = t_a_temp; t_p_val(id_fold) = t_p_temp;
                teta_val(id_fold) = teta_temp; b_g_val(id_fold) = b_g_temp; b_r_val(id_fold) = b_r_temp;
                max_corr = corr_train(id_params, id_fold);
                ind_max = id_param;
            end
        end
        
        %%%compute validation error (correlations)
        distance_glob_temp=zeros(n_mouse, n_conc_level, length(val_trials));
        distance_rel_temp=zeros(n_mouse, n_conc_level, length(val_trials));
        distance_stm_temp=zeros(n_mouse, n_conc_level, length(val_trials));
        distance_behav_val=zeros(length(val_trials),1);
        for id_trial = 1:length(val_trials)
            id1_odor = trials_pair(val_trials(id_trial),1);
            id2_odor = trials_pair(val_trials(id_trial),2);
            
            distance_behav_val(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
            
            for id_mouse = 1:length(Dff) %iterate over mice
                for conc_level_id = 1:n_conc_level %iterate over conc level
                    conc_level = conc_level_def(conc_level_id);
                    this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
                    n_glom = size(this_dff,1);
                    this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
                    this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
                    
                    pattern_1 = this_lat(:,id1_odor+1);
                    pattern_2 = this_lat(:,id2_odor+1);
                    
                    [dist_glob, dist_rel, dist_stm] = compute_neural_distance_single_pair(pattern_1, pattern_2,...
                        t_a_val(id_fold), t_p_val(id_fold), t_g_val(id_fold), ...
                        teta_val(id_fold), b_g_val(id_fold), b_r_val(id_fold), time);
                    
                    distance_glob_temp(id_mouse, conc_level_id, id_trial) = dist_glob;
                    distance_rel_temp(id_mouse, conc_level_id, id_trial) = dist_rel;
                    distance_stm_temp(id_mouse, conc_level_id, id_trial) = dist_stm;
                    
                end
            end
        end
        distance_glob = squeeze(nanmean(nanmean(distance_glob_temp,2),1));
        distance_rel = squeeze(nanmean(nanmean(distance_glob_temp,2),1));
        distance_stm = squeeze(nanmean(nanmean(distance_glob_temp,2),1));
        
        corr_val(id_fold) = corr(distance_stm(:), distance_behav_val(:));
    end
    
    [corr_val_max, id_best_fold] = max(corr_val);
    t_g_best = t_g_val(id_best_fold); t_a_best = t_a_val(id_best_fold);
    t_p_best = t_p_val(id_best_fold); teta_best = teta_val(id_best_fold);
    b_r_best = b_r_val(id_best_fold); b_g_best = b_g_val(id_best_fold);
    
else
    train_trials = train_id;
    
    for id_trial = 1:length(train_trials)
        id1_odor = trials_pair(train_trials(id_trial),1);
        id2_odor = trials_pair(train_trials(id_trial),2);
        distance_behav_train(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
    end
    
    for id_param = 1:length(t_g) %%%cycle over parameters to optimize
        disp(['Params: ' num2str(id_param) '/' num2str(length(t_g))]);
        t_g_temp = t_g(id_param); t_a_temp = t_a(id_param);
        t_p_temp = t_p(id_param); teta_temp = teta(id_param);
        b_g_temp = b_g(id_param); b_r_temp = b_r(id_param);
        
        distance_glob_temp=zeros(n_mouse, n_conc_level, length(train_trials));
        distance_rel_temp=zeros(n_mouse, n_conc_level, length(train_trials));
        distance_stm_temp=zeros(n_mouse, n_conc_level, length(train_trials));
        for id_trial = 1:length(train_trials)
            id1_odor = trials_pair(train_trials(id_trial),1);
            id2_odor = trials_pair(train_trials(id_trial),2);
            
            for id_mouse = 1:length(Dff) %iterate over mice
                for conc_level_id = 1:n_conc_level %iterate over conc level
                    conc_level = conc_level_def(conc_level_id);
                    this_dff = Dff{id_mouse}(:,:,:,conc_level); %extract df/f
                    n_glom = size(this_dff,1);
                    this_lat = Lat{id_mouse}(:,:,conc_level)+t0; %extract latencies
                    this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
                    
                    pattern_1 = this_lat(:,id1_odor+1);
                    pattern_2 = this_lat(:,id2_odor+1);
                    
                    [dist_glob, dist_rel, dist_stm] = compute_neural_distance_single_pair(pattern_1, pattern_2,...
                        t_a_temp, t_p_temp, t_g_temp, teta_temp, b_g_temp, b_r_temp, time);
                    
                    distance_glob_temp(id_mouse, conc_level_id, id_trial) = dist_glob;
                    distance_rel_temp(id_mouse, conc_level_id, id_trial) = dist_rel;
                    distance_stm_temp(id_mouse, conc_level_id, id_trial) = dist_stm;
                end
            end
        end
        distance_glob = squeeze(nanmean(nanmean(distance_glob_temp,2),1));
        distance_rel = squeeze(nanmean(nanmean(distance_rel_temp,2),1));
        distance_stm = squeeze(nanmean(nanmean(distance_stm_temp,2),1));
        
        corr_train(id_param) = corr(distance_stm(:), distance_behav_train(:));
        if id_param == 1
            t_g_best = t_g_temp; t_a_best = t_a_temp; t_p_best = t_p_temp;
            teta_best = teta_temp; b_g_best = b_g_temp; b_r_best = b_r_temp;
            max_corr = corr_train(id_param);
            ind_max = id_param;
        elseif corr_train(id_param)>max_corr
            t_g_best = t_g_temp; t_a_best = t_a_temp; t_p_best = t_p_temp;
            teta_best = teta_temp; b_g_best = b_g_temp; b_r_best = b_r_temp;
            max_corr = corr_train(id_params);
            ind_max = id_param;
        end
    end
end

%%
%%%find best on validation set and compute correlation on test set