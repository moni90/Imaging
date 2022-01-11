%% set fixed params
%target patterns
target = [10:40:210 NaN*ones(1,6)];
%some example of probe trials
probe_T = [30 50 80 130 170 150 NaN*ones(1,6)];
probe_G = target + 30;
probe_S = [10 NaN 90 130 170 NaN NaN 50 NaN NaN NaN 210];
probes = [probe_T; probe_G; probe_S];
%model parameters
tau_glob = 100;
tau_act = 60;
tau_prim = 200;
teta = pi/2;

tt = 0:0.1:800;

%%

f_target = from_discrete_to_waveform(target, tau_act, tau_prim, tt);
f_target_rot = apply_rotation(f_target, teta);
CM_target = compute_CM(f_target,tt);

for i_trial = 1:size(probes,1)
    
    pattern = probes(i_trial,:);
    f_probe = from_discrete_to_waveform(pattern, tau_act, tau_prim, tt);
    figure;
    for i_spot = 1:size(f_probe,1)
        subplot(2,6,i_spot);
        plot(tt,f_probe(i_spot,:),'r','linewidth',2);
        hold on; plot(tt,f_target(i_spot,:),'k','linewidth',2);
        ylim([0 1]); xlim([0 500])
        xlabel('time (ms)'); title(['Channel ' num2str(i_spot)])
        if i_spot == 6
            legend('probe','target');
        end
    end
    f_probe_rot = apply_rotation(f_probe, teta);
    CM_probe = compute_CM(f_probe,tt);   
%     figure;
%     for i_spot = 1:size(f_target_rot,1)
%         subplot(2,6,i_spot);
%         plot(tt,f_probe_rot(i_spot,:),'r','linewidth',2);
%         hold on; plot(tt,f_target_rot(i_spot,:),'k','linewidth',2);
%         ylim([0 1]);
%         xlabel('time (ms)'); title(['Channel ' num2str(i_spot)])
%     end

    [dist_glob, dist_rel] = compute_dist(f_target_rot, f_probe_rot, CM_target, CM_probe, tau_glob, tt);
    disp(['global distance = ' num2str(dist_glob)]);
    disp(['relative distance = ' num2str(dist_rel)]);
end
