function f = single_exp(t_onset, tau, tt)
%single exponential function
% t_onset(t_onset==0)=NaN;
n_trials = length(t_onset);
if size(t_onset,2)>1 %transform t_onset into a column vector
    t_onset = t_onset';
end
if size(tt,1)>1 %transform tt into a row vector
    tt = tt';
end
tt = repmat(tt,n_trials,1);
f = exp(-(tt-t_onset)/tau).*((tt-t_onset)>=0);
f(isnan(f))=0;