clear;
close all;
fieldtrip_path = "C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507";

addpath(fieldtrip_path);
ft_defaults;

load("C:\Users\fkamdar\Desktop\repos\lfp_toolbox\derivatives\sub-CPEEG02_ses-01_data_clean_neu4-6.mat")

% Plot layout: left column = LEFT hemisphere, right column = RIGHT.
%plot_layout = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT', 'ONE_THREE_RIGHT'};
plot_layout = { 'ZERO_TWO_LEFT'; 'ZERO_TWO_RIGHT'; 'ONE_THREE_RIGHT'};
label_map = containers.Map();
label_map('ZERO_TWO_LEFT')   = 'ventral ch - left hemisphere';
label_map('ZERO_TWO_RIGHT')  = 'ventral ch - right hemisphere';
label_map('ONE_THREE_RIGHT') = 'dorsal ch - right hemisphere';

fprintf('Loaded %d clean trials, channels: %s\n', ...
        numel(data_clean.trial), strjoin(data_clean.label, ', '));

%% TFR
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = logspace(log10(2), log10(200));    % cfg.foi = 1:1:99;   % current: linear
cfg.toi = -0.5:0.001:3;
cfg.keeptrials = 'yes';
freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
tfr_baseline = [-0.5 -0.2];     % s, pre-image baseline window
tfr_basetype = 'relchange';     % 'relchange' (fractional) or 'db' % change

cfg = [];
cfg.baseline = tfr_baseline;
cfg.baselinetype = tfr_basetype; 
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);



%% Cluster permutation stats on full time x frequency:  Neg vs Neu
chan_name = plot_layout{1};   % 
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');

%Split freq_clean_bl by condition, keeping full freq x time (no averaging)
cfg = [];
cfg.trials = neg_idx;
freq_neg = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = neu_idx;
freq_neu = ft_selectdata(cfg, freq_clean_bl);

%Stats config
cfg = [];
cfg.channel           = {chan_name};
cfg.frequency        = 'all';     
cfg.latency           = 'all';       
cfg.method            = 'montecarlo';
cfg.statistic         = 'indepsamplesT';
cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.3; 
cfg.clusterstatistic  = 'maxsum';
cfg.tail              = 0;
cfg.clustertail       = 0;
cfg.alpha             = 0.05;      
cfg.numrandomization  = 1000;

n_neg = numel(freq_neg.trialinfo);
n_neu = numel(freq_neu.trialinfo);
cfg.design = [ones(1, size(freq_neg.powspctrm,1)), ...
              ones(1, size(freq_neu.powspctrm,1)) * 2];
cfg.ivar   = 1;

stat_freq_neg_vs_neu = ft_freqstatistics(cfg, freq_neg, freq_neu);

%Inspect
% fprintf('Positive clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
% fprintf('Negative clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)
% 
% disp([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob])
% disp([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob])


fprintf('Positive clusters:\n'); disp(stat_freq_neg_vs_neu.posclusters)
fprintf('Negative clusters:\n'); disp(stat_freq_neg_vs_neu.negclusters)

disp([stat_freq_neg_vs_neu.negclusters.prob])
disp([stat_freq_neg_vs_neu.posclusters.prob])

%% Cluster permutation stats on full time x frequency:  Pos vs Neu
chan_name = plot_layout{3};   % 
pos_idx = strcmp(trial_info_clean.Valence, 'Pos');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');

%Split freq_clean_bl by condition, keeping full freq x time (no averaging)
cfg = [];
cfg.trials = pos_idx;
freq_pos = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = neu_idx;
freq_neu = ft_selectdata(cfg, freq_clean_bl);

%Stats config
cfg = [];
cfg.channel           = {chan_name};
cfg.frequency        = [2 29];     
cfg.latency           = [0 2];       
cfg.method            = 'montecarlo';
cfg.statistic         = 'indepsamplesT';
cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.3; 
cfg.clusterstatistic  = 'maxsum';
cfg.tail              = 0;
cfg.clustertail       = 0;
cfg.alpha             = 0.05;      
cfg.numrandomization  = 10000;

n_pos = numel(freq_pos.trialinfo);
n_neu = numel(freq_neu.trialinfo);
cfg.design = [ones(1, size(freq_pos.powspctrm,1)), ...
              ones(1, size(freq_neu.powspctrm,1)) * 2];
cfg.ivar   = 1;

stat_freq_pos_vs_neu = ft_freqstatistics(cfg, freq_pos, freq_neu);

fprintf('Positive clusters:\n'); disp(stat_freq_pos_vs_neu.posclusters)
fprintf('Negative clusters:\n'); disp(stat_freq_pos_vs_neu.negclusters)

disp([stat_freq_pos_vs_neu.negclusters.prob])
disp([stat_freq_pos_vs_neu.posclusters.prob])