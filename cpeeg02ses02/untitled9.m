%% =========================================================================
% BASELINE CORRECTION (PER TRIAL)
% =========================================================================

cfg = [];
cfg.demean          = 'yes';

% ---- DEFINE BASELINE WINDOW ----
% Example: pre-stim baseline (adjust based on your epoch timing)
cfg.baselinewindow  = [-0.5 0];   % change if needed

data_bl = ft_preprocessing(cfg, data_clean);

fprintf('Baseline correction applied: [%0.2f %0.2f] sec\n', cfg.baselinewindow);