function params = build_params(model, t_a, t_p, t_g, teta, b_g, b_r, b_0)

if strcmp(model, 'STM_CM') || strcmp(model, 'STM_First') || strcmp(model, 'STM_First_analytics') || strcmp(model, 'STM_First_ampl_analytics') || strcmp(model, 'STM_First_ampl')
    params.t_a = t_a; params.t_p = t_p; params.t_g = t_g; params.teta = teta;
    params.b_g = b_g; params.b_r = b_r; params.b_0 = b_0;
end