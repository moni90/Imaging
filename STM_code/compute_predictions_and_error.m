function [y_pred, err_pred] = compute_predictions_and_error(model, ...
    Lat, Dff, trials_pair, id_keep_mouse, id_keep_conc_level, id_trials, t0, frame_period, time,...
    y_true, params, fit_metric)

% compute y_pred and error for selected trials
if strcmp(model, 'STM_CM')
    t_a = params.t_a; t_p = params.t_p; t_g = params.t_g; teta = params.teta;
    b_g = params.b_g; b_r = params.b_r; b_0 = params.b_0;
    align_mode = 'CM';

    [distance_glob, distance_rel, y_pred] = compute_avg_neural_distance(...
        Lat, trials_pair, id_keep_mouse, id_keep_conc_level, id_trials, align_mode, ...
        t_a, t_p, t_g, teta, b_g, b_r, b_0, t0, frame_period, time);
    
    err_pred = prediction_error( y_true(:), y_pred(:), fit_metric);
    
elseif strcmp(model, 'STM_First')
    t_a = params.t_a; t_p = params.t_p; t_g = params.t_g; teta = params.teta;
    b_g = params.b_g; b_r = params.b_r; b_0 = params.b_0;
    align_mode = 'First';

    [distance_glob, distance_rel, y_pred] = compute_avg_neural_distance(...
        Lat, trials_pair, id_keep_mouse, id_keep_conc_level, id_trials, align_mode, ...
        t_a, t_p, t_g, teta, b_g, b_r, b_0, t0, frame_period, time);
    
    err_pred = prediction_error( y_true(:), y_pred(:), fit_metric);
    
elseif strcmp(model, 'STM_First_analytics')
    t_a = params.t_a; t_p = params.t_p; t_g = params.t_g; teta = params.teta;
    b_g = params.b_g; b_r = params.b_r; b_0 = params.b_0;
    align_mode = 'First';

    [distance_glob, distance_rel, y_pred] = compute_avg_neural_distance_analytics(...
        Lat, trials_pair, id_keep_mouse, id_keep_conc_level, id_trials, align_mode, ...
        t_a, t_p, t_g, teta, b_g, b_r, b_0, t0, frame_period, time);
    
    err_pred = prediction_error( y_true(:), y_pred(:), fit_metric);
    
elseif strcmp(model, 'STM_First_ampl_analytics')
    t_a = params.t_a; t_p = params.t_p; t_g = params.t_g; teta = params.teta;
    b_g = params.b_g; b_r = params.b_r; b_0 = params.b_0;
    align_mode = 'First';

    [distance_glob, distance_rel, y_pred] = compute_avg_neural_distance_with_amplitude_analytics(...
        Lat, Dff, trials_pair, id_keep_mouse, id_keep_conc_level, id_trials, align_mode, ...
        t_a, t_p, t_g, teta, b_g, b_r, b_0, t0, frame_period, time);
    
    err_pred = prediction_error( y_true(:), y_pred(:), fit_metric);
    
elseif strcmp(model, 'STM_First_ampl')
    t_a = params.t_a; t_p = params.t_p; t_g = params.t_g; teta = params.teta;
    b_g = params.b_g; b_r = params.b_r; b_0 = params.b_0;
    align_mode = 'First';

    [distance_glob, distance_rel, y_pred] = compute_avg_neural_distance_with_amplitude(...
        Lat, Dff, trials_pair, id_keep_mouse, id_keep_conc_level, id_trials, align_mode, ...
        t_a, t_p, t_g, teta, b_g, b_r, b_0, t0, frame_period, time);
    
    err_pred = prediction_error( y_true(:), y_pred(:), fit_metric);
end