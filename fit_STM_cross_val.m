% Fit STM model on neural data using average behavioral distance
% The parameters tau_act, tau_prim, tau_global are considered as
% hyperparameters
% The parameters b_g, b_r are considered as parameters

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

%% IMPORT BEHAVIORAL DISTANCE DATA
behavioral_data = load(behavioral_dist_path);
behavioral_distance_mat = behavioral_data.distance_mat_avg; %behavioral_data.distance_mat; %
behavioral_distance_mat(isnan(behavioral_distance_mat))=0;
%% DEFINE TRAINING, VALIDATION AND TEST TRIALS
n_odors = 24;
[~, ~, trials_pair] = lower_half(ones(n_odors, n_odors));

n_trials = n_odors*(n_odors+1)/2;
train_fract = 0.8;
n_train = floor(train_fract*n_trials);
n_folds = 5;
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

%% DEFINE PARAMETERS GRID
%model parameters
% tau_glob = 100*(1e-1); tau_act = 60*(1e-1); tau_prim = 200*(1e-1);
% teta = pi/2; beta_glob = 1.1821; beta_rel = 0.0263;

tau_glob = [50 100 150]*(1e-1);
tau_act = [20 40 60 80]*(1e-1);
tau_prim = [50 100 150 200]*(1e-1);
teta_angle = pi/2;
beta_glob = [0 0.01 0.05 1];% 1.2];
beta_rel = [0.01 0.1 1];% 0.1];

[t_g, t_a, t_p, teta] = ndgrid(tau_glob, tau_act, tau_prim, teta_angle); %hyperparameters
t_g = t_g(:);t_a = t_a(:); t_p = t_p(:); teta = teta(:);
[b_g, b_r] = ndgrid(beta_glob, beta_rel); %parameters
b_g = b_g(:); b_r = b_r(:);
%% AVERAGE NEURAL RESPONSE FIT
%Compure neural distance matrix separately for each animal, average across
%animals and maximize correlation between average neural distance and
%behavioral neural distance
fps = 100;
frame_period = 1/fps;
time = -200:1:1899;
time = time*frame_period;
t0 = find(time==0);
n_odors = 24;
n_mouse = length(Dff);
id_keep_mouse = 1:1:n_mouse;
% n_conc_level = size(Dff{1},4);
id_keep_conc_level = [2];

corr_train = zeros(length(b_g),n_folds);
corr_val = zeros(length(t_g),n_folds);
t_g_val = zeros(n_folds,1); t_a_val = zeros(n_folds,1); t_p_val = zeros(n_folds,1);
teta_val = zeros(n_folds,1); b_g_val = zeros(n_folds,1); b_r_val = zeros(n_folds,1);
max_corr_val = zeros(n_folds,1); ind_val_max = zeros(n_folds,1);
    
for id_fold = 1:n_folds
    train_trials = train_id(folds_id~=id_fold);
    val_trials = train_id(folds_id==id_fold);
    
    distance_behav_train = zeros(length(train_trials),1);
    for id_trial = 1:length(train_trials)
        id1_odor = trials_pair(train_trials(id_trial),1);
        id2_odor = trials_pair(train_trials(id_trial),2);
        distance_behav_train(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
    end
    distance_behav_val = zeros(length(val_trials),1);
    for id_trial = 1:length(val_trials)
        id1_odor = trials_pair(val_trials(id_trial),1);
        id2_odor = trials_pair(val_trials(id_trial),2);
        distance_behav_val(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
    end
    
    for id_hyper_param = 1:length(t_g)
        disp(['VALIDATION. Fold ' num2str(id_fold) ' of ' num2str(n_folds) '. Params: ' num2str(id_hyper_param) '/' num2str(length(t_g))]);
        t_g_temp = t_g(id_hyper_param); t_a_temp = t_a(id_hyper_param);
        t_p_temp = t_p(id_hyper_param); teta_temp = teta(id_hyper_param);

        for id_param = 1:length(b_g) %%%cycle over parameters to optimize
%             disp(['TRAINING. Fold ' num2str(id_fold) ' of ' num2str(n_folds) '. Params: ' num2str(id_param) '/' num2str(length(b_g))]);

            b_g_temp = b_g(id_param); b_r_temp = b_r(id_param);

            [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
                Lat, trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, ...
                t_a_temp, t_p_temp, t_g_temp, teta_temp, b_g_temp, b_r_temp, t0, frame_period, time);
            
            corr_train(id_param, id_fold) = corr(distance_stm(:), distance_behav_train(:));
            if id_param == 1
                b_g_val(id_fold) = b_g_temp; b_r_val(id_fold) = b_r_temp;
                max_corr = corr_train(id_param, id_fold);
                ind_max = id_param;
            elseif corr_train(id_param, id_fold)>max_corr
                b_g_val(id_fold) = b_g_temp; b_r_val(id_fold) = b_r_temp;
                max_corr = corr_train(id_params, id_fold);
                ind_max = id_param;
            end
        end
        
        %compute validation error (correlations)
        [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
                Lat, trials_pair, id_keep_mouse, id_keep_conc_level, val_trials, ...
                t_a_temp, t_p_temp, t_g_temp, teta_temp, b_g_val(id_fold), b_r_val(id_fold), t0, frame_period, time);
        
        corr_val(id_hyper_param, id_fold) = corr(distance_stm(:), distance_behav_val(:));
        if id_hyper_param == 1
            t_g_val(id_fold) = t_g_temp; t_p_val(id_fold) = t_p_temp;
            t_a_val(id_fold) = t_a_temp; teta_val(id_fold) = teta_temp;
            max_corr_val(id_fold) = corr_val(id_hyper_param, id_fold);
            ind_val_max(id_fold) = id_hyper_param;
        elseif corr_val(id_hyper_param, id_fold)>max_corr_val
            t_g_val(id_fold) = t_g_temp; t_p_val(id_fold) = t_p_temp;
            t_a_val(id_fold) = t_a_temp; teta_val(id_fold) = teta_temp;
            max_corr_val(id_fold) = corr_val(id_hyper_param, id_fold);
            ind_val_max(id_fold) = id_hyper_param;
        end
    end
    
end

[~, ind_best] = max(max_corr_val);
t_g_best = t_g_val(ind_best); t_p_best = t_p_val(ind_best);
t_a_best = t_a_val(ind_best); teta_best = teta_val(ind_best);

%re-fit beta_global and beta_rel using the full training set for the best
%cross validated hyperparameters
train_trials = train_id;
distance_behav_train = zeros(length(train_trials),1);
for id_trial = 1:length(train_trials)
    id1_odor = trials_pair(train_trials(id_trial),1);
    id2_odor = trials_pair(train_trials(id_trial),2);
    distance_behav_train(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
end

corr_train = zeros(length(b_g),1);
for id_param = 1:length(b_g) %%%cycle over parameters to optimize
    disp(['Params: ' num2str(id_param) '/' num2str(length(b_g))]);
    b_g_temp = b_g(id_param); b_r_temp = b_r(id_param);
    
    [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
        Lat, trials_pair, id_keep_mouse, id_keep_conc_level, train_trials, ...
        t_a_best, t_p_best, t_g_best, teta_best, b_g_temp, b_r_temp, t0, frame_period, time);
    
    corr_train(id_param) = corr(distance_stm(:), distance_behav_train(:));
    if id_param == 1
        b_g_best = b_g_temp; b_r_best = b_r_temp;
        corr_train_def = corr_train(id_param);
        ind_max = id_param;
    elseif corr_train(id_param)>corr_train_def
        b_g_best = b_g_temp; b_r_best = b_r_temp;
        corr_train_def = corr_train(id_params);
        ind_max = id_param;
    end
end

%%
%%%find best on validation set and compute correlation on test set
test_trials = test_id;
distance_behav_test = zeros(length(test_trials),1);
for id_trial = 1:length(test_trials)
    id1_odor = trials_pair(test_trials(id_trial),1);
    id2_odor = trials_pair(test_trials(id_trial),2);
    distance_behav_test(id_trial) = behavioral_distance_mat(id1_odor,id2_odor);
end

[distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance(...
    Lat, trials_pair, id_keep_mouse, id_keep_conc_level, test_trials, ...
    t_a_best, t_p_best, t_g_best, teta_best, b_g_best, b_r_best, t0, frame_period, time);

corr_test = corr(distance_stm(:), distance_behav_test(:));
disp(corr_test);