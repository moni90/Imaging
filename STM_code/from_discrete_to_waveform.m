function f_pattern = from_discrete_to_waveform(pattern, tau_act, tau_prim, tt)

    %activation of single channels
    f_pattern = single_exp(pattern, tau_act, tt);
    %primacy modulation
    onset_primacy_pattern = nanmin(pattern);
    primacy_kernel_pattern = exp(-(tt-onset_primacy_pattern)/tau_prim).*((tt-onset_primacy_pattern)>=0);
    f_pattern = f_pattern.*repmat(primacy_kernel_pattern,size(f_pattern,1),1);
       
end
