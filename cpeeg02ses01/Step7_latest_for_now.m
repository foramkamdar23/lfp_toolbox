
fieldtrip_path = "C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507";

addpath(fieldtrip_path);
ft_defaults;

%load("C:\Users\fkamdar\Desktop\repos\lfp_toolbox\derivatives\sub-CPEEG02_ses-01_data_clean.mat")

% Plot layout: left column = LEFT hemisphere, right column = RIGHT.
%plot_layout = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT', 'ONE_THREE_RIGHT'};
plot_layout = { 'ZERO_TWO_LEFT'; 'ZERO_TWO_RIGHT'; 'ONE_THREE_RIGHT'};
label_map = containers.Map();
label_map('ZERO_TWO_LEFT')   = 'ventral ch - left hemisphere';
label_map('ZERO_TWO_RIGHT')  = 'ventral ch - right hemisphere';
label_map('ONE_THREE_RIGHT') = 'dorsal ch - right hemisphere';

% fprintf('Loaded %d clean trials, channels: %s\n', ...
%         numel(data_clean.trial), strjoin(data_clean.label, ', '));

%% Cluster permutation stats on full time x frequency: FEEL Pos vs Neu
chan_name = plot_layout{1};   % ventral channel right hemisphere

%Split freq_clean_bl by condition, keeping full freq x time (no averaging)
cfg = [];
cfg.trials = feel_pos_idx;
freq_feelpos = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = feel_neu_idx;
freq_feelneu = ft_selectdata(cfg, freq_clean_bl);

%Stats config
cfg = [];
cfg.channel           = {chan_name};
cfg.frequency        = [2 29];      % frequency range of interest
%cfg.frequency         = 'all';
cfg.latency           = [0 2];       % full image window
cfg.method            = 'montecarlo';
cfg.statistic         = 'indepsamplesT';
cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.05; 0.3
cfg.clusterstatistic  = 'maxsum';
cfg.tail              = 0;
cfg.clustertail       = 0;
cfg.alpha             = 0.025;     0.05  % two-sided
cfg.numrandomization  = 10000;

n_pos = numel(freq_feelpos.trialinfo);
n_neu = numel(freq_feelneu.trialinfo);
cfg.design = [ones(1, size(freq_feelpos.powspctrm,1)), ...
              ones(1, size(freq_feelneu.powspctrm,1)) * 2];
cfg.ivar   = 1;

stat_freq_feelpos_vs_feelneu_ZTR = ft_freqstatistics(cfg, freq_feelpos, freq_feelneu);

%Inspect
% fprintf('Positive clusters:\n'); disp(stat_freq_feelpos_vs_feelneu_ZTR.posclusters)
% fprintf('Negative clusters:\n'); disp(stat_freq_feelpos_vs_feelneu_ZTR.negclusters)
% disp([stat_freq_feelpos_vs_feelneu_ZTR.negclusters.prob])


%% Cluster permutation stats on full time x frequency: FEEL Neg vs Neu


%Split freq_clean_bl by condition, keeping full freq x time (no averaging)
cfg = [];
cfg.trials = feel_neg_idx;
freq_feelneg = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = feel_neu_idx;
freq_feelneu = ft_selectdata(cfg, freq_clean_bl);

%Stats config
cfg = [];
cfg.channel           = {chan_name};
cfg.frequency        = [2 29];      % frequency range of interest
cfg.latency           = [0 3];       % full image window
cfg.method            = 'montecarlo';
cfg.statistic         = 'indepsamplesT';
cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.05;
cfg.clusterstatistic  = 'maxsum';
cfg.tail              = 0;
cfg.clustertail       = 0;
cfg.alpha             = 0.025;       % two-sided
cfg.numrandomization  = 10000;

n_neg = numel(freq_feelneg.trialinfo);
n_neu = numel(freq_feelneu.trialinfo);
cfg.design = [ones(1, size(freq_feelneg.powspctrm,1)), ...
              ones(1, size(freq_feelneu.powspctrm,1)) * 2];
cfg.ivar   = 1;

stat_freq_feelneg_vs_feelneu_ZTR = ft_freqstatistics(cfg, freq_feelneg, freq_feelneu);

%Inspect
% fprintf('Positive clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
% fprintf('Negative clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)
% 
% disp([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob])
% disp([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob])


%% Plot stats as TFR: t-values with cluster outline (Pos vs Neu, Neg vs Neu)

sig_alpha = 0.025;

figure('Name', sprintf('%s: Cluster stats TFR', label_map(chan_name)));

% --- Pos vs Neu ---
subplot(1,2,1)
imagesc(stat_freq_feelpos_vs_feelneu_ZTR.time, ...
        stat_freq_feelpos_vs_feelneu_ZTR.freq, ...
        squeeze(stat_freq_feelpos_vs_feelneu_ZTR.stat));
axis xy;
xline(0, 'k--');
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim_val = max(abs(squeeze(stat_freq_feelpos_vs_feelneu_ZTR.stat)), [], 'all');
clim([-clim_val clim_val]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Emotion Perception Condition: Positive vs Neutral');

% Outline significant clusters (pos)
hold on;
if ~isempty(stat_freq_feelpos_vs_feelneu_ZTR.posclusters)
    sig_k = find([stat_freq_feelpos_vs_feelneu_ZTR.posclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelpos_vs_feelneu_ZTR.posclusterslabelmat == k);
        contour(stat_freq_feelpos_vs_feelneu_ZTR.time, ...
                stat_freq_feelpos_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
% Outline significant clusters (neg)
if ~isempty(stat_freq_feelpos_vs_feelneu_ZTR.negclusters)
    sig_k = find([stat_freq_feelpos_vs_feelneu_ZTR.negclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelpos_vs_feelneu_ZTR.negclusterslabelmat == k);
        contour(stat_freq_feelpos_vs_feelneu_ZTR.time, ...
                stat_freq_feelpos_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
hold off;

% --- Neg vs Neu ---
subplot(1,2,2)
imagesc(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
        stat_freq_feelneg_vs_feelneu_ZTR.freq, ...
        squeeze(stat_freq_feelneg_vs_feelneu_ZTR.stat));
axis xy;
xline(0, 'k--');
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim_val = max(abs(squeeze(stat_freq_feelneg_vs_feelneu_ZTR.stat)), [], 'all');
clim([-clim_val clim_val]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Emotion Perception Condition: Negative vs Neutral');

hold on;
if ~isempty(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
    sig_k = find([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelneg_vs_feelneu_ZTR.posclusterslabelmat == k);
        contour(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
                stat_freq_feelneg_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
if ~isempty(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)
    sig_k = find([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelneg_vs_feelneu_ZTR.negclusterslabelmat == k);
        contour(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
                stat_freq_feelneg_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
hold off;

%sgtitle(sprintf('%s: t-statistic TFR)', label_map(chan_name), sig_alpha), 'Interpreter', 'none');
sgtitle(sprintf('%s: t-statistic TFR', label_map(chan_name)), 'Interpreter', 'none');


fprintf('Positive clusters:\n'); disp(stat_freq_feelpos_vs_feelneu_ZTR.posclusters)
fprintf('Negative clusters:\n'); disp(stat_freq_feelpos_vs_feelneu_ZTR.negclusters)
disp([stat_freq_feelpos_vs_feelneu_ZTR.negclusters.prob])

fprintf('Positive clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
fprintf('Negative clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)

disp([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob])
disp([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob])