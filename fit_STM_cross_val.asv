% Fit STM model on neural data using average behavioral distance
% The parameters tau_act, tau_prim, tau_global are considered as hyperparameters
% The parameters b_g, b_r are considered as parameters
%% Path to code and data folder
code_path = 'D:\mmoroni\github_repos\Imaging'; addpath(genpath(code_path));
working_dir = pwd;
path_to_dff = 'D:\mmoroni\DATA\olfaction_metric\1P\data\Behavior_odors\Set1';
file_name = 'Latency_data_8pures_16mixtures_022220.mat';
behavioral_dist_path = 'D:\mmoroni\Hiro_project\analyses\Behavior\Odor_set1\pooled_distance_matrix.mat';
%% LOAD PROCESSED IMAGING DATA
load(fullfile(path_to_dff, file_name));
%% IMPORT BEHAVIORAL DISTANCE DATA
behavioral_data = load(behavioral_dist_path);
behavioral_distance_mat = behavioral_data.distance_mat; %behavioral_data.distance_mat_avg; %
behavioral_distance_mat(isnan(behavioral_distance_mat))=0;
%% DEFINE TRAINING, VALIDATION AND TEST TRIALS
rng(1)
n_odors = 24;
train_fract = 0.8;
n_folds = 5;

[~, ~, trials_pair] = lower_half(ones(n_odors, n_odors));
n_trials = n_odors*(n_odors+1)/2;
n_train = floor(train_fract*n_trials);
trials_per_fold = floor(n_train/n_folds);
n_train = n_folds*trials_per_fold;
train_fract = n_train/n_trials;
test_fract = 1-train_fract;
n_test = n_trials - n_train;
train_id = randperm(n_trials, n_train);
test_id = setdiff(1:1:n_trials,train_id);
folds_id = [];
for i = 1:n_folds
    folds_id = [folds_id; i*ones(trials_per_fold,1)];
end
folds_id = folds_id(randperm(length(folds_id)));
%% CHOOSE MODEL AND METRIC
model = 'STM_First_ampl_analytics'; %'STM_First_ampl'; %'STM_CM'; %'STM_First_analytics'; %'STM_First'; %
fit_metric = 'sse'; %sum squared errors
% fit_metric = 'corr'; %pearson correlation
%% DEFINE PARAMETERS GRID
tau_glob = [60:20:120]*(1e-1);%[60:20:120]*(1e-1);
tau_act = [30:10:50]*(1e-1);%[60:10:90]*(1e-1);
tau_prim = [20:20:100]*(1e-1);%[100:20:180]*(1e-1);
teta_angle = pi/2;
beta_glob = [2:2:10];% 1.2];
beta_rel = [500:25:600];%[0.08:0.01:0.11];% 0.1];
beta_0 = [-0.02:0.02:0];%[-0.5:0.1:-0.3];%

[t_g, t_a, t_p, teta] = ndgrid(tau_glob, tau_act, tau_prim, teta_angle); %hyperparameters
t_g = t_g(:);t_a = t_a(:); t_p = t_p(:); teta = teta(:);
[b_g, b_r, b_0] = ndgrid(beta_glob, beta_rel, beta_0); %parameters
b_g = b_g(:); b_r = b_r(:); b_0 = b_0(:);
%% AVERAGE NEURAL RESPONSE FIT
%Compute neural distance matrix separately for each animal, average across
%animals and minimize difference between average neural distance and
%behavioral neural distance
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
n_mouse = length(Dff);
id_keep_mouse = 1:1:n_mouse;
id_keep_conc_level = [2]; %1: low concentration, 2: high concentration

corr_train = zeros(length(b_g),n_folds); corr_val = zeros(length(t_g),n_folds);
t_g_val = zeros(n_folds,1); t_a_val = zeros(n_folds,1); t_p_val = zeros(n_folds,1); teta_val = zeros(n_folds,1);
b_g_val = zeros(n_folds,1); b_r_val = zeros(n_folds,1); b_0_val = zeros(n_folds,1);
min_corr_val = zeros(n_folds,1); ind_val_min = zeros(n_folds,1);
    
for id_fold = 1:n_folds
    train_trials = train_id(folds_id~=id_fold);
    val_trials = train_id(folds_id==id_fold);
    
    distance_behav_train = select_y_true(behavioral_distance_mat, trials_pair, train_trials);
    distance_behav_val = select_y_true(behavioral_distance_mat, trials_pair, val_trials);

    for id_hyper_param = 1:length(t_g)
        disp(['VALIDATION. Fold ' num2str(id_fold) ' of ' num2str(n_folds) '. Params: ' num2str(id_hyper_param) '/' num2str(length(t_g))]);

        for id_param = 1:length(b_g) %%%cycle over parameters to optimize
            
            params = build_params(model, t_a(id_hyper_param), t_p(id_hyper_param),...
                t_g(id_hyper_param), teta(id_hyper_param),...
                b_g(id_param), b_r(id_param), b_0(id_param));
            [~, corr_train(id_param, id_fold)] = compute_predictions_and_error(model, Lat, Dff,...
                trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, t0, frame_period, time, ...
                distance_behav_train(:), params, fit_metric);
           
            if id_param == 1
                b_g_val(id_fold) = b_g(id_param); b_r_val(id_fold) = b_r(id_param); b_0_val(id_fold) = b_0(id_param);
                min_corr = corr_train(id_param, id_fold);
                ind_min = id_param;
            elseif corr_train(id_param, id_fold)<min_corr
                b_g_val(id_fold) = b_g(id_param); b_r_val(id_fold) = b_r(id_param); b_0_val(id_fold) = b_0(id_param);
                min_corr = corr_train(id_param, id_fold);
                ind_min = id_param;
            end
        end
        
        %compute validation error
        params = build_params(model,t_a(id_hyper_param), t_p(id_hyper_param),...
            t_g(id_hyper_param), teta(id_hyper_param),...
            b_g_val(id_fold), b_r_val(id_fold), b_0_val(id_fold));
        [~, corr_val(id_hyper_param, id_fold)] = compute_predictions_and_error(model, Lat, Dff,...
            trials_pair, id_keep_mouse, id_keep_conc_level, val_trials, t0, frame_period, time, ...
            distance_behav_val(:), params, fit_metric);
       
        if id_hyper_param == 1
            t_g_val(id_fold) = t_g(id_hyper_param); t_p_val(id_fold) = t_p(id_hyper_param);
            t_a_val(id_fold) = t_a(id_hyper_param); teta_val(id_fold) = teta(id_hyper_param);
            min_corr_val(id_fold) = corr_val(id_hyper_param, id_fold);
            ind_val_min(id_fold) = id_hyper_param;
        elseif corr_val(id_hyper_param, id_fold)<min_corr_val(id_fold)
            t_g_val(id_fold) = t_g(id_hyper_param); t_p_val(id_fold) = t_p(id_hyper_param);
            t_a_val(id_fold) = t_a(id_hyper_param); teta_val(id_fold) = teta(id_hyper_param);
            min_corr_val(id_fold) = corr_val(id_hyper_param, id_fold);
            ind_val_min(id_fold) = id_hyper_param;
        end
    end
    
end

[~, ind_best] = min(min_corr_val);
t_g_best = t_g_val(ind_best); t_p_best = t_p_val(ind_best);
t_a_best = t_a_val(ind_best); teta_best = teta_val(ind_best);

%re-fit beta_global and beta_rel using the full training set for the best
%cross validated hyperparameters
train_trials = train_id;
distance_behav_train = select_y_true(behavioral_distance_mat, trials_pair, train_trials);

corr_train = zeros(length(b_g),1);
for id_param = 1:length(b_g) %%%cycle over parameters to optimize
    disp(['Params: ' num2str(id_param) '/' num2str(length(b_g))]);
    
    params = build_params(model, t_a_best, t_p_best, t_g_best, teta_best,...
        b_g(id_param), b_r(id_param), b_0(id_param));
    [distance_stm, corr_train(id_param)] = compute_predictions_and_error(model, Lat, Dff,...
        trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, t0, frame_period, time, ...
        distance_behav_train(:), params, fit_metric);

    if id_param == 1
        b_g_best = b_g(id_param); b_r_best = b_r(id_param); b_0_best = b_0(id_param);
        corr_train_def = corr_train(id_param);
        ind_min = id_param;
    elseif corr_train(id_param)<corr_train_def
        b_g_best = b_g(id_param); b_r_best = b_r(id_param); b_0_best = b_0(id_param);
        corr_train_def = corr_train(id_param);
        ind_min = id_param;
    end
end

%% COMPUTE METRIC ON TEST SET AND PLOT RESULTS
test_trials = test_id;
distance_behav_test = select_y_true(behavioral_distance_mat, trials_pair, test_trials);

params = build_params(model, t_a_best, t_p_best, t_g_best, teta_best, b_g_best, b_r_best, b_0_best);
[distance_stm_test, corr_test] = compute_predictions_and_error(model, Lat, Dff,...
    trials_pair, id_keep_mouse, id_keep_conc_level, test_trials, t0, frame_period, time, ...
    distance_behav_test(:), params, fit_metric);

figure;
scatter(distance_stm(:), distance_behav_train(:));
hold on; scatter(distance_stm_test(:), distance_behav_test(:));
hold on; plot(0:.1:1, 0:.1:1,'k--');
xlabel('Neural distance'); ylabel('Behavioral distance'); title([model, '. ', fit_metric, '=', num2str(corr_test)]);

%% SAVE FIT PARAMETERS
save_path = 'D:\mmoroni\Hiro_project\analyses\STM_model\';
save_name = ['bestParams_', model,  '_', fit_metric, '.mat'];
save(fullfile(save_path, save_name),'t_a_best', 't_p_best', 't_g_best', ...
    'teta_best', 'b_g_best', 'b_r_best', 'b_0_best', 'model', 'fit_metric',...
    'distance_stm', 'distance_behav_train', 'distance_stm_test', 'distance_behav_test');

%% REPLICATE BEHAVIORAL DISTANCE TESTS

only_test = 0;

% histogram
p_match = [];
p_A_B = [];
p_A_AB = [];
p_AC_BC = [];
if only_test == 0
    train_trials = train_id;
    trials_id = trials_pair(train_trials,:);  
    for i_trial = 1:length(train_trials)
        odor_components = trials_pair(train_trials(i_trial),:);
        if odor_components(1) == odor_components(2)
            p_match = [p_match; distance_stm(i_trial)];
        elseif (odor_components(1)<=8 ) && (odor_components(2)<=8 )
            p_A_B = [p_A_B; distance_stm(i_trial)];
        elseif (odor_components(1)<=8 ) && (odor_components(2)>8 )
            p_A_AB = [p_A_AB; distance_stm(i_trial)];
        elseif (odor_components(1)>8 ) && (odor_components(2)>8 )
            p_AC_BC = [p_AC_BC; distance_stm(i_trial)];
        end
    end
end

test_trials = test_id;
trials_id = trials_pair(test_trials,:);
for i_trial = 1:length(test_trials)
    odor_components = trials_pair(test_trials(i_trial),:);
    if odor_components(1) == odor_components(2)
        p_match = [p_match; distance_stm_test(i_trial)];
    elseif (odor_components(1)<=8 ) && (odor_components(2)<=8 ) 
        p_A_B = [p_A_B; distance_stm_test(i_trial)];
    elseif (odor_components(1)<=8 ) && (odor_components(2)>8 ) 
        p_A_AB = [p_A_AB; distance_stm_test(i_trial)];
    elseif (odor_components(1)>8 ) && (odor_components(2)>8 )
        p_AC_BC = [p_AC_BC; distance_stm_test(i_trial)];
    end
end

p_mean = [nanmean(p_match) nanmean(p_A_B) nanmean(p_AC_BC) nanmean(p_A_AB)];
p_std = [nanstd(p_match) nanstd(p_A_B) nanstd(p_AC_BC) nanstd(p_A_AB)];
p_num = [length(p_match) length(p_A_B) length(p_AC_BC) length(p_A_AB)];

figure;
errorbar([1,2,3,4], p_mean,p_std./sqrt(p_num));
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'match','A-B','AC-BC','A-AB'});

%build correspondence for scatter pairs
odors_noAir = odors(2:end);
scatter_A_AB_pair = [];
scatter_AC_BC_pair = [];
for i1 = 1:8
    odor1 = odors_noAir{i1};
    for i2 = i1+1:8
        odor2 = odors_noAir{i2};
        for i3 = 9:length(odors_noAir)
            odor3 = odors_noAir{i3};
            if contains(odor3, odor1) && contains(odor3, odor2)
%                 disp([odor1, ', ' odor2, ', ', odor3])
                scatter_A_AB_pair = [scatter_A_AB_pair; i1, i2, i3];
            end
        end
        for i4 = 1:8
            if (i4 ~= i1) && (i4~=i2)
                odor4 = odors_noAir{i4};
                for i3 = 9:length(odors_noAir)
                    odor3 = odors_noAir{i3};
                    for i5 = 9:length(odors_noAir)
                        odor5 = odors_noAir{i5};
                        if contains(odor3, odor1) && contains(odor3, odor4) && contains(odor5, odor2) && contains(odor5, odor4)
%                             disp([odor1, ', ' odor2, ', ', odor4, ', ', odor3, ', ', odor5])
                            scatter_AC_BC_pair = [scatter_AC_BC_pair; i1, i2, i3, i5];
                        end
                    end
                end
            end
        end
    end
end

%scatter plot
scatter_A_B = [];
scatter_A_AB = [];

if only_test == 0
    train_trials = train_id;
    trials_id = trials_pair(train_trials,:);
    for i_A_AB = 1:size(scatter_A_AB_pair,1)
        A_B_temp = scatter_A_AB_pair(i_A_AB,[1 2]);
        A_AB_temp1 = scatter_A_AB_pair(i_A_AB,[1 3]);
        A_AB_temp2 = scatter_A_AB_pair(i_A_AB,[2 3]);
        
        id1 = find(sum(abs(trials_id-A_B_temp),2)==0);
        id2 = find(sum(abs(trials_id-A_AB_temp1),2)==0);
        id3 = find(sum(abs(trials_id-A_AB_temp2),2)==0);
        if ~isempty(id1) && ~isempty(id2)
            scatter_A_B = [scatter_A_B; distance_stm(id1)];
            scatter_A_AB = [scatter_A_AB; distance_stm(id2)];
        end
        if ~isempty(id1) && ~isempty(id3)
            scatter_A_B = [scatter_A_B; distance_stm(id1)];
            scatter_A_AB = [scatter_A_AB; distance_stm(id3)];
        end
    end
end

test_trials = test_id;
trials_id = trials_pair(test_trials,:);
for i_A_AB = 1:size(scatter_A_AB_pair,1)
    A_B_temp = scatter_A_AB_pair(i_A_AB,[1 2]);
    A_AB_temp1 = scatter_A_AB_pair(i_A_AB,[1 3]);
    A_AB_temp2 = scatter_A_AB_pair(i_A_AB,[2 3]);
    
    id1 = find(sum(abs(trials_id-A_B_temp),2)==0);
    id2 = find(sum(abs(trials_id-A_AB_temp1),2)==0);
    id3 = find(sum(abs(trials_id-A_AB_temp2),2)==0);
    if ~isempty(id1) && ~isempty(id2)
        scatter_A_B = [scatter_A_B; distance_stm_test(id1)];
        scatter_A_AB = [scatter_A_AB; distance_stm_test(id2)];
    end
    if ~isempty(id1) && ~isempty(id3)
        scatter_A_B = [scatter_A_B; distance_stm_test(id1)];
        scatter_A_AB = [scatter_A_AB; distance_stm_test(id3)];
    end
end

figure; scatter(scatter_A_B, scatter_A_AB);
hold on; plot([0:.1:1],[0:.1:1],'k--');
xlabel('A-B'); ylabel('A-AB');

scatter_A_B = [];
scatter_AC_BC = [];

if only_test == 0
    train_trials = train_id;
    trials_id = trials_pair(train_trials,:);
    for i_AC_BC = 1:size(scatter_AC_BC_pair,1)
        A_B_temp = scatter_AC_BC_pair(i_AC_BC,[1 2]);
        AC_BC_temp = scatter_AC_BC_pair(i_AC_BC,[3 4]);
        
        id1 = find(sum(abs(trials_id-A_B_temp),2)==0);
        id2 = find(sum(abs(trials_id-AC_BC_temp),2)==0);
        if ~isempty(id1) && ~isempty(id2)
            scatter_A_B = [scatter_A_B; distance_stm(id1)];
            scatter_AC_BC = [scatter_AC_BC; distance_stm(id2)];
        end
    end
end

test_trials = test_id;
trials_id = trials_pair(test_trials,:);

for i_AC_BC = 1:size(scatter_AC_BC_pair,1)
    A_B_temp = scatter_AC_BC_pair(i_AC_BC,[1 2]);
    AC_BC_temp = scatter_AC_BC_pair(i_AC_BC,[3 4]);
    
    id1 = find(sum(abs(trials_id-A_B_temp),2)==0);
    id2 = find(sum(abs(trials_id-AC_BC_temp),2)==0);
    if ~isempty(id1) && ~isempty(id2)
        scatter_A_B = [scatter_A_B; distance_stm_test(id1)];
        scatter_AC_BC = [scatter_AC_BC; distance_stm_test(id2)];
    end
end

figure; scatter(scatter_A_B, scatter_AC_BC);
hold on; plot([0:.1:1],[0:.1:1],'k--');
xlabel('A-B'); ylabel('AC-BC');

%% IMAGING DISTANCE MATRIX
imaging_distance_mat = zeros(size(behavioral_distance_mat));
train_trials = train_id;
trials_id = trials_pair(train_trials,:);
for i = 1:length(train_trials)
    imaging_distance_mat(trials_id(i,1), trials_id(i,2)) =  distance_stm(i);
    imaging_distance_mat(trials_id(i,2), trials_id(i,1)) =  distance_stm(i);
end
test_trials = test_id;
trials_id = trials_pair(test_trials,:);
for i = 1:length(test_trials)
    imaging_distance_mat(trials_id(i,1), trials_id(i,2)) =  distance_stm_test(i);
    imaging_distance_mat(trials_id(i,2), trials_id(i,1)) =  distance_stm_test(i);
end
figure;
subplot(1,2,1); imagesc(behavioral_distance_mat); colorbar; colormap(gray);
subplot(1,2,2); imagesc(imaging_distance_mat); colorbar; colormap(gray);

