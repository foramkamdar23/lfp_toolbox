
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



%% Cluster permutation stats: FEEL&Pos vs FEEL&Neu
chan_name = plot_layout{3}; % ventral channel right hemisphere
 
frequency_range = [1 40];   % full range; narrow later if you have a prior

Stats_structure.cfg        = [];
Stats_structure.individual = [];
Stats_structure.time       = [];
Stats_structure.label      = [];
Stats_structure.dimord     = 'subj_chan_time';

% FEEL & Pos
tl_feelpos = Stats_structure;
tl_feelpos.time  = freq_axis;
tl_feelpos.label = {chan_name};
tl_feelpos.individual(:,1,:) = pow_feelpos.(chan_name);

% FEEL & Neu
tl_feelneu = Stats_structure;
tl_feelneu.time  = freq_axis;
tl_feelneu.label = {chan_name};
tl_feelneu.individual(:,1,:) = pow_feelneu.(chan_name);

% Stats config
cfg = [];
cfg.channel          = {chan_name};
cfg.neighbours       = [];
cfg.latency          = frequency_range;
cfg.method           = 'montecarlo';
cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;       % two-sided
cfg.numrandomization = 10000;

cfg.design = [ones(1, size(tl_feelpos.individual,1)), ...
              ones(1, size(tl_feelneu.individual,1)) * 2];
cfg.ivar   = 1;

stat_feelpos_vs_feelneu_ZTR = ft_timelockstatistics(cfg, tl_feelpos, tl_feelneu);

%save('Stats_FeelPosVsNeu_ZERO_TWO_RIGHT', 'stat_feelpos_vs_feelneu_ZTR')



%% Cluster permutation stats: FEEL&Neg vs FEEL&Neu
chan_name = plot_layout{3}; % ventral channel right hemisphere
 
frequency_range = [1 40];   % full range; narrow later if you have a prior

Stats_structure.cfg        = [];
Stats_structure.individual = [];
Stats_structure.time       = [];
Stats_structure.label      = [];
Stats_structure.dimord     = 'subj_chan_time';

% FEEL & Pos
tl_feelneg = Stats_structure;
tl_feelneg.time  = freq_axis;s
tl_feelneg.label = {chan_name};
tl_feelneg.individual(:,1,:) = pow_feelneg.(chan_name);

% % FEEL & Neu
% tl_feelneu = Stats_structure;
% tl_feelneu.time  = freq_axis;
% tl_feelneu.label = {chan_name};
% tl_feelneu.individual(:,1,:) = pow_feelneu.(chan_name);

% Stats config
cfg = [];
cfg.channel          = {chan_name};
cfg.neighbours       = [];
cfg.latency          = frequency_range;
cfg.method           = 'montecarlo';
cfg.statistic        = 'indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;       % two-sided
cfg.numrandomization = 10000;

cfg.design = [ones(1, size(tl_feelneg.individual,1)), ...
              ones(1, size(tl_feelneu.individual,1)) * 2];
cfg.ivar   = 1;

stat_feelneg_vs_feelneu = ft_timelockstatistics(cfg, tl_feelneg, tl_feelneu);

%save('Stats_FeelPosVsNeu_ZERO_TWO_RIGHT', 'stat_feelpos_vs_feelneu_ZTR')
fprintf('Positive clusters:\n'); disp(stat_feelneg_vs_feelneu_ZTR.posclusters)
fprintf('Negative clusters:\n'); disp(stat_feelneg_vs_feelneu_ZTR.negclusters)
fprintf('Negative cluster freq bins:\n'); disp(freq_axis(find(stat_feelneg_vs_feelneu_ZTR.negclusterslabelmat == 1)))