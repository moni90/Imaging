function align_ref = compute_align_ref(f_target,time_pattern, align_mode)

if strcmp(align_mode,'CM')
    %sum activity across channels
    f_sum = sum(f_target,1);
    %compute covered area
    covered_area1_cum = cumtrapz(time_pattern,f_sum,2);
    %find time point where half of the total area is covered
    [~,ind_min] = min((covered_area1_cum-covered_area1_cum(end)/2).^2);
    align_ref = time_pattern(ind_min);
elseif strcmp(align_mode,'First')
    %sum activity across channels
    f_sum = sum(f_target,1);
    align_ref = time_pattern(find(f_sum>0,1,'first'));

end