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

%% SOME FIXED VARIABLES
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
n_mouse = length(Dff);
id_keep_mouse = 1:1:n_mouse;
id_keep_conc_level = [2];
%% PERMUTATION OF BEHAVIORAL DISTANCE VALUES
n_permutations = 100;
test_error = zeros(n_permutations,1);
for id_permute = 1:n_permutations
    behavioral_distance_mat_sh = shuffle_symmetric_matrix(behavioral_distance_mat);

    %AVERAGE NEURAL RESPONSE FIT
    corr_train = zeros(length(b_g),n_folds); corr_val = zeros(length(t_g),n_folds);
    t_g_val = zeros(n_folds,1); t_a_val = zeros(n_folds,1); t_p_val = zeros(n_folds,1); teta_val = zeros(n_folds,1);
    b_g_val = zeros(n_folds,1); b_r_val = zeros(n_folds,1); b_0_val = zeros(n_folds,1);
    min_corr_val = zeros(n_folds,1); ind_val_min = zeros(n_folds,1);
    
    for id_fold = 1:n_folds
        train_trials = train_id(folds_id~=id_fold);
        val_trials = train_id(folds_id==id_fold);
        
        distance_behav_train = select_y_true(behavioral_distance_mat_sh, trials_pair, train_trials);
        distance_behav_val = select_y_true(behavioral_distance_mat_sh, trials_pair, val_trials);
        
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
    distance_behav_train = select_y_true(behavioral_distance_mat_sh, trials_pair, train_trials);
    
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
    
    % COMPUTE METRIC ON TEST SET AND PLOT RESULTS
    test_trials = test_id;
    distance_behav_test = select_y_true(behavioral_distance_mat_sh, trials_pair, test_trials);
    
    params = build_params(model, t_a_best, t_p_best, t_g_best, teta_best, b_g_best, b_r_best, b_0_best);
    [distance_stm_test, corr_test] = compute_predictions_and_error(model, Lat, Dff,...
        trials_pair, id_keep_mouse, id_keep_conc_level, test_trials, t0, frame_period, time, ...
        distance_behav_test(:), params, fit_metric);
    test_error(id_permute) = corr_test;
end
