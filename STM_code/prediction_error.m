function pred_error = prediction_error(y_true, y_predict, fit_metric)

%I want a metric to be minimized
if strcmp(fit_metric, 'sse')
    pred_error = mean((y_true(:)-y_predict(:)).^2);
elseif strcmp(fit_metric, 'corr')
    pred_error = -corr(y_true(:),y_predict(:));
end
