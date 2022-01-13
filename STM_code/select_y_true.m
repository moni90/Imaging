function y_true = select_y_true(y_all, trials_pair, selected_trials)

y_true = zeros(length(selected_trials),1);
for id_trial = 1:length(selected_trials)
    id1_odor = trials_pair(selected_trials(id_trial),1);
    id2_odor = trials_pair(selected_trials(id_trial),2);
    y_true(id_trial) = y_all(id1_odor,id2_odor);
end