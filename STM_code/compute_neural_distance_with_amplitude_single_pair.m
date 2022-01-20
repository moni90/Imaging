function [dist_glob, dist_rel, dist_stm] = compute_neural_distance_with_amplitude_single_pair(odor1, odor2, amplitude_1, amplitude_2, align_mode, tau_act, tau_prim, tau_glob, teta, beta_glob, beta_rel, beta_0, tt)

    amplitude_1 = amplitude_1/(nanmax(amplitude_1(:)));
    amplitude_2 = amplitude_2/(nanmax(amplitude_2(:)));
    
    f_odor1 = from_discrete_to_waveform_with_amplitude(odor1, amplitude_1, tau_act, tau_prim, tt);
    f_odor2 = from_discrete_to_waveform_with_amplitude(odor2, amplitude_2, tau_act, tau_prim, tt);
    %apply rotation to glomerular id
    f_odor1_rot = apply_rotation(f_odor1, teta);    
    f_odor2_rot = apply_rotation(f_odor2, teta);
    %compute alignment_reference
%     CM_odor1 = compute_CM(f_odor1,tt);
%     CM_odor2 = compute_CM(f_odor2,tt);
    CM_odor1 = compute_align_ref(f_odor1,tt, align_mode);
    CM_odor2 = compute_align_ref(f_odor2,tt, align_mode);
    
    dt_CM = (CM_odor1-CM_odor2);
    dist_glob = (1-exp(-abs(dt_CM)/tau_glob));
    
    %compute relative lag between patterns and align to CM
    time_step = tt(2)-tt(1);
    shift_rel = dt_CM;
    shift_rel = round(shift_rel/time_step); 
    if shift_rel >= 0 %CM target later than CM spot. Shift CM_spot later.
        f_odor2_zero_pad = [zeros(size(f_odor2_rot,1),abs(shift_rel)) f_odor2_rot];
        f_odor1_zero_pad = [f_odor1_rot zeros(size(f_odor1_rot,1),abs(shift_rel)) ];
    else
        f_odor1_zero_pad = [zeros(size(f_odor1_rot,1),abs(shift_rel)) f_odor1_rot];
        f_odor2_zero_pad = [f_odor2_rot zeros(size(f_odor2_rot,1),abs(shift_rel)) ];
    end
    f_diff = (f_odor1_zero_pad-f_odor2_zero_pad).^2;
    diff_channels = cumtrapz(time_step,f_diff,2);
    dist_rel = sum(diff_channels(:,end));
    
    dist_stm = beta_glob*dist_glob + beta_rel*dist_rel;
    dist_stm = 1./(1+exp(-log(max(0,dist_stm+beta_0))));%

end