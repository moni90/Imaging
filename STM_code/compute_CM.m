function CM = compute_CM(f_target,time_pattern)
%sum activity across channels
f_sum = sum(f_target,1);
%compute covered area
covered_area1_cum = cumtrapz(time_pattern,f_sum,2);
%find time point where half of the total area is covered
[~,ind_min] = min((covered_area1_cum-covered_area1_cum(end)/2).^2);
CM = time_pattern(ind_min);

end