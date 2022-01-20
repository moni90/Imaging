function [dist_glob, dist_rel, dist_stm] = compute_neural_distance_with_amplitude_single_pair_analytics(odor1, odor2, A_odor1, A_odor2, align_mode, tau_act, tau_prim, tau_glob, teta, beta_glob, beta_rel, beta_0)

    t_min_odor1 = nanmin(odor1);
    t_min_odor2 = nanmin(odor2);
    t_norm_odor1 = odor1 - t_min_odor1;
    t_norm_odor2 = odor2 - t_min_odor2;
%     C_odor1 = exp(-t_norm_odor1/tau_prim + (odor1 + t_min_odor1)/tau_act);
%     C_odor2 = exp(-t_norm_odor2/tau_prim + (odor2 + t_min_odor2)/tau_act);
    %in case response amplitude is considered
    C_odor1 = A_odor1.*exp(-t_norm_odor1/tau_prim + (odor1 + t_min_odor1)/tau_act);
    C_odor2 = A_odor2.*exp(-t_norm_odor2/tau_prim + (odor2 + t_min_odor2)/tau_act);
    
    
    ind_odor1 = 1*(t_norm_odor1<=t_norm_odor2);
    ind_odor2 = 1 - ind_odor1;
    
    t_min = ind_odor1.*(t_norm_odor1) + ind_odor2.*(t_norm_odor2);
    t_max = ind_odor1.*(t_norm_odor2) + ind_odor2.*(t_norm_odor1);
    C_min = ind_odor1.*(C_odor1) + ind_odor2.*(C_odor2);
    
    dist_rel = nanmean((C_min.^2).*exp(-2*t_min/tau_act) + ( ( (C_odor1-C_odor2).^2 - C_min.^2 ).*exp(-2*t_max/tau_act)  ));
    dist_glob = (1-exp(-abs(t_min_odor2-t_min_odor1)/tau_glob));

    dist_stm = beta_glob*dist_glob + beta_rel*dist_rel;
    dist_stm = 1./(1+exp(-log(max(0,dist_stm+beta_0))));%

end