function [distance_glob, distance_rel, distance_stm] = compute_avg_neural_distance_analytics(...
    Lat, trials_pair, id_keep_mouse, id_keep_conc_level, id_keep_trials, align_mode, ...
    t_a, t_p, t_g, teta, b_g, b_r, b_0, t0, frame_period, time)

n_mouse = length(id_keep_mouse);
n_conc_level = length(id_keep_conc_level);
n_trials = length(id_keep_trials);

distance_glob_temp=zeros(n_mouse, n_conc_level, n_trials);
distance_rel_temp=zeros(n_mouse, n_conc_level, n_trials);
distance_stm_temp=zeros(n_mouse, n_conc_level, n_trials);

for i_mouse = 1:n_mouse %iterate over mice
    id_mouse = id_keep_mouse(i_mouse);
    
    for i_conc_level = 1:n_conc_level %iterate over conc level
        id_conc_level = id_keep_conc_level(i_conc_level);
        this_lat = Lat{id_mouse}(:,:,id_conc_level)+t0; %extract latencies
        this_lat = (this_lat-200)*frame_period; %convert to seconds from odor presentations
        
        parfor id_trial = 1:n_trials
            id1_odor = trials_pair(id_keep_trials(id_trial),1);
            id2_odor = trials_pair(id_keep_trials(id_trial),2);
            
            pattern_1 = this_lat(:,id1_odor+1);
            pattern_2 = this_lat(:,id2_odor+1);
            
            [dist_glob, dist_rel, dist_stm] = compute_neural_distance_single_pair_analytics(pattern_1, pattern_2, align_mode,...
                t_a, t_p, t_g, teta, b_g, b_r, b_0);
            
            distance_glob_temp(i_mouse, i_conc_level, id_trial) = dist_glob;
            distance_rel_temp(i_mouse, i_conc_level, id_trial) = dist_rel;
            distance_stm_temp(i_mouse, i_conc_level, id_trial) = dist_stm;
        end
    end
end
distance_glob = squeeze(nanmean(nanmean(distance_glob_temp,2),1));
distance_rel = squeeze(nanmean(nanmean(distance_rel_temp,2),1));
distance_stm = squeeze(nanmean(nanmean(distance_stm_temp,2),1));