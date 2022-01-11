function [dist_glob, dist_rel] = compute_dist(f_target_rot, f_probe_rot, CM_target, CM_probe, tau_glob, tt)

    dt_CM = (CM_target-CM_probe);
    dist_glob = (1-exp(-abs(dt_CM)/tau_glob));
    
    %compute relative lag between patterns and align to CM
    time_step = tt(2)-tt(1);
    shift_rel = dt_CM;
    shift_rel = round(shift_rel/time_step); 
    if shift_rel >= 0 %CM target later than CM spot. Shift CM_spot later.
        f_probe_zero_pad = [zeros(size(f_probe_rot,1),abs(shift_rel)) f_probe_rot];
        f_target_zero_pad = [f_target_rot zeros(size(f_target_rot,1),abs(shift_rel)) ];
    else
        f_target_zero_pad = [zeros(size(f_target_rot,1),abs(shift_rel)) f_target_rot];
        f_probe_zero_pad = [f_probe_rot zeros(size(f_probe_rot,1),abs(shift_rel)) ];
    end
    f_diff = (f_target_zero_pad-f_probe_zero_pad).^2;
    diff_channels = cumtrapz(time_step,f_diff,2);
    dist_rel = sum(diff_channels(:,end));
    

end