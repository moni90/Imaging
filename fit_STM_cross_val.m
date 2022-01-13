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
behavioral_distance_mat = behavioral_data.distance_mat_avg; %behavioral_data.distance_mat; %
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
model = 'STM_CM'; %'STM_First'; %
fit_metric = 'sse'; %sum squared errors
% fit_metric = 'corr'; %pearson correlation
%% DEFINE PARAMETERS GRID
tau_glob = [80]*(1e-1);
tau_act = [60:5:90]*(1e-1);
tau_prim = [120:10:180]*(1e-1);
teta_angle = pi/2;
beta_glob = [0];% 1.2];
beta_rel = [0.08:0.01:0.12];% 0.1];
beta_0 = [-0.6:0.05:-0.4];

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
    
%     distance_behav_train = zeros(length(train_trials),1);
%     for id_trial = 1:length(train_trials)
%         id1_odor = trials_pair(train_trials(id_trial),1);
%         id2_odor = trials_pair(train_trials(id_trial),2);
%         distance_behav_train(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
%     end
%     distance_behav_val = zeros(length(val_trials),1);
%     for id_trial = 1:length(val_trials)
%         id1_odor = trials_pair(val_trials(id_trial),1);
%         id2_odor = trials_pair(val_trials(id_trial),2);
%         distance_behav_val(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
%     end
    
    for id_hyper_param = 1:length(t_g)
        disp(['VALIDATION. Fold ' num2str(id_fold) ' of ' num2str(n_folds) '. Params: ' num2str(id_hyper_param) '/' num2str(length(t_g))]);
%         t_g_temp = t_g(id_hyper_param); t_a_temp = t_a(id_hyper_param);
%         t_p_temp = t_p(id_hyper_param); teta_temp = teta(id_hyper_param);

        for id_param = 1:length(b_g) %%%cycle over parameters to optimize
            
            params = build_params(model, t_a(id_hyper_param), t_p(id_hyper_param),...
                t_g(id_hyper_param), teta(id_hyper_param),...
                b_g(id_param), b_r(id_param), b_0(id_param));
            [~, corr_train(id_param, id_fold)] = compute_predictions_and_error(model, Lat, ...
                trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, t0, frame_period, time, ...
                distance_behav_train(:), params, fit_metric);

%             b_g_temp = b_g(id_param); b_r_temp = b_r(id_param); b_0_temp = b_0(id_param);
%             [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
%                 Lat, trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, ...
%                 t_a_temp, t_p_temp, t_g_temp, teta_temp, b_g_temp, b_r_temp, b_0_temp, t0, frame_period, time);
%             corr_train(id_param, id_fold) = prediction_error(distance_behav_train(:), distance_stm(:), fit_metric);
            
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
        [~, corr_val(id_hyper_param, id_fold)] = compute_predictions_and_error(model, Lat, ...
            trials_pair, id_keep_mouse, id_keep_conc_level, val_trials, t0, frame_period, time, ...
            distance_behav_val(:), params, fit_metric);
        
%         [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
%                 Lat, trials_pair, id_keep_mouse, id_keep_conc_level, val_trials, ...
%                 t_a_temp, t_p_temp, t_g_temp, teta_temp, b_g_val(id_fold),...
%                 b_r_val(id_fold), b_0_val(id_fold), t0, frame_period, time);       
%         corr_val(id_hyper_param, id_fold) = prediction_error(distance_behav_val(:), distance_stm(:), fit_metric);
        
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
% distance_behav_train = zeros(length(train_trials),1);
% for id_trial = 1:length(train_trials)
%     id1_odor = trials_pair(train_trials(id_trial),1);
%     id2_odor = trials_pair(train_trials(id_trial),2);
%     distance_behav_train(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
% end

corr_train = zeros(length(b_g),1);
for id_param = 1:length(b_g) %%%cycle over parameters to optimize
    disp(['Params: ' num2str(id_param) '/' num2str(length(b_g))]);
    
    params = build_params(model, t_a_best, t_p_best, t_g_best, teta_best,...
        b_g(id_param), b_r(id_param), b_0(id_param));
    [distance_stm, corr_train(id_param)] = compute_predictions_and_error(model, Lat, ...
        trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, t0, frame_period, time, ...
        distance_behav_train(:), params, fit_metric);
        
%     b_g_temp = b_g(id_param); b_r_temp = b_r(id_param); b_0_temp = b_0(id_param);
%     [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
%         Lat, trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, ...
%         t_a_best, t_p_best, t_g_best, teta_best, b_g_temp, b_r_temp, b_0_temp, t0, frame_period, time);
%     corr_train(id_param) = prediction_error(distance_behav_train(:), distance_stm(:), fit_metric);
    
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
% distance_behav_test = zeros(length(test_trials),1);
% for id_trial = 1:length(test_trials)
%     id1_odor = trials_pair(test_trials(id_trial),1);
%     id2_odor = trials_pair(test_trials(id_trial),2);
%     distance_behav_test(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
% end

params = build_params(model, t_a_best, t_p_best, t_g_best, teta_best, b_g_best, b_r_best, b_0_best);
[distance_stm_test, corr_test] = compute_predictions_and_error(model, Lat, ...
    trials_pair, id_keep_mouse, id_keep_conc_level, test_trials, t0, frame_period, time, ...
    distance_behav_test(:), params, fit_metric);
    
% [distance_glob, distance_rel, distance_stm_test] = compute_avg_neural_distance(...
%     Lat, trials_pair, id_keep_mouse, id_keep_conc_level, test_trials, ...
%     t_a_best, t_p_best, t_g_best, teta_best, b_g_best, b_r_best, b_0_best, t0, frame_period, time);
% corr_test = prediction_error(distance_behav_test(:), distance_stm_test(:), fit_metric);

figure;
scatter(distance_stm(:), distance_behav_train(:));
hold on; scatter(distance_stm_test(:), distance_behav_test(:));
hold on; plot(0:.1:1, 0:.1:1,'k--');
xlabel('Neural distance'); ylabel('Behavioral distance'); title([model, '. ', fit_metric, '=', num2str(corr_test)]);