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
cfg.foi = logspace(log10(2), log10(29));    % cfg.foi = 1:1:99;   % current: linear
cfg.toi = -0.5:0.01:3;
cfg.keeptrials = 'yes';
freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
tfr_baseline = [-0.5 -0.2];     % s, pre-image baseline window
tfr_basetype = 'relchange';     % 'relchange' (fractional) or 'db' % change

cfg = [];
cfg.baseline = tfr_baseline;
cfg.baselinetype = tfr_basetype; 
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

%avg
cfg = [];
cfg.avgoverrpt = 'yes';
freq_avg = ft_selectdata(cfg, freq_clean_bl);

% PLOT TFR
subj_id = 'CPEEG02';
sess_cond = 'off stim';

tfr_clim = [-2 2];

plot_tfr_grid(freq_avg, plot_layout, tfr_clim, ...
              sprintf('%s, %s', subj_id, sess_cond), label_map);

% function plot_tfr_grid(freq, layout, clim, ttl)
%     % Per-channel TFR grid (LEFT col = left hemisphere, RIGHT col = right).
%     % freq.powspctrm is [channels x freq x time] after averaging over trials.
%     [nr, nc] = size(layout);
%     figure('Name', ttl);
%     for r = 1:nr
%         for c = 1:nc
%             lab = layout{r,c};
%             ch = find(strcmp(freq.label, lab), 1);
%             subplot(nr, nc, (r-1)*nc + c);
%             imagesc(freq.time, freq.freq, squeeze(freq.powspctrm(ch,:,:)), clim);
%             axis xy;
%             xline(0, 'k--');                 % image onset
%             title(lab, 'Interpreter', 'none');
%             xlabel('Time (s)'); ylabel('Frequency (Hz)');
%             colorbar;
%         end
%     end
%     colormap jet
%     sgtitle(ttl);
% end
% 

%% TFR BY instruction
feel_idx = strcmp(trial_info_clean.Instruction, 'FEEL');
tone_idx = strcmp(trial_info_clean.Instruction, 'TONE');

%Average separately by instruction
cfg = [];
cfg.trials = feel_idx;
cfg.avgoverrpt = 'yes';
freq_avg_feel = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = tone_idx;
cfg.avgoverrpt = 'yes';
freq_avg_tone = ft_selectdata(cfg, freq_clean_bl);

% Plot
plot_tfr_grid(freq_avg_feel, plot_layout, tfr_clim, ...
              sprintf('%s, %s, FEEL', subj_id, sess_cond), label_map);

plot_tfr_grid(freq_avg_tone, plot_layout, tfr_clim, ...
              sprintf('%s, %s, TONE', subj_id, sess_cond), label_map);

%% TFR POS AND Instruction
pos_idx = strcmp(trial_info_clean.Valence, 'Pos');

feel_pos_idx = feel_idx & pos_idx;
tone_pos_idx = tone_idx & pos_idx;


%Average separately by instruction, positive valence only
cfg = [];
cfg.trials = feel_pos_idx;
cfg.avgoverrpt = 'yes';
freq_avg_feel_pos = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = tone_pos_idx;
cfg.avgoverrpt = 'yes';
freq_avg_tone_pos = ft_selectdata(cfg, freq_clean_bl);

%plot
plot_tfr_grid(freq_avg_feel_pos, plot_layout, tfr_clim, ...
              sprintf('%s, %s, FEEL, Pos', subj_id, sess_cond), label_map);

plot_tfr_grid(freq_avg_tone_pos, plot_layout, tfr_clim, ...
              sprintf('%s, %s, TONE, Pos', subj_id, sess_cond), label_map);

%% TFR Neu AND Instruction
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');

feel_neu_idx = feel_idx & neu_idx;
tone_neu_idx = tone_idx & neu_idx;


%Average separately by instruction, positive valence only
cfg = [];
cfg.trials = feel_neu_idx;
cfg.avgoverrpt = 'yes';
freq_avg_feel_neu = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = tone_neu_idx;
cfg.avgoverrpt = 'yes';
freq_avg_tone_neu = ft_selectdata(cfg, freq_clean_bl);

%plot
plot_tfr_grid(freq_avg_feel_neu, plot_layout, tfr_clim, ...
              sprintf('%s, %s, FEEL, Neu', subj_id, sess_cond), label_map);

plot_tfr_grid(freq_avg_tone_neu, plot_layout, tfr_clim, ...
              sprintf('%s, %s, TONE, Neu', subj_id, sess_cond), label_map);

%% TFR Neg AND Instuction
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');


feel_neg_idx = feel_idx & neg_idx;
tone_neg_idx = tone_idx & neg_idx;


%Average separately by instruction, positive valence only
cfg = [];
cfg.trials = feel_neg_idx;
cfg.avgoverrpt = 'yes';
freq_avg_feel_neg = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = tone_neg_idx;
cfg.avgoverrpt = 'yes';
freq_avg_tone_neg = ft_selectdata(cfg, freq_clean_bl);

%plot
plot_tfr_grid(freq_avg_feel_neg, plot_layout, tfr_clim, ...
              sprintf('%s, %s, FEEL, Neg', subj_id, sess_cond), label_map);

plot_tfr_grid(freq_avg_tone_neg, plot_layout, tfr_clim, ...
              sprintf('%s, %s, TONE, Neg', subj_id, sess_cond), label_map);

%% 
%% Extract trial-level power per condition, per channel (for stats)

toi_stats = [0 3];   % adjust based on your hypothesis (e.g. [0 1] for early effects)
tidx = freq_clean_bl.time >= toi_stats(1) & freq_clean_bl.time <= toi_stats(2);

freq_axis = freq_clean_bl.freq; 


for c = 1:numel(plot_layout)
    chan_name = plot_layout{c};
    ch = find(strcmp(freq_clean_bl.label, chan_name), 1);

    % powspctrm dims: [trial x channel x freq x time]

    % FEEL & Pos
    pow_feelpos.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(feel_pos_idx, ch, :, tidx), 4));   % -> [n_trials x n_freqs]

    % FEEL & Neu
    pow_feelneu.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(feel_neu_idx, ch, :, tidx), 4));

     % FEEL & Neu
    pow_feelneg.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(feel_neg_idx, ch, :, tidx), 4));

end


for c = 1:numel(plot_layout)
    chan_name = plot_layout{c};
    ch = find(strcmp(freq_clean_bl.label, chan_name), 1);

    % powspctrm dims: [trial x channel x freq x time]

    % FEEL & Pos
    pow_tonepos.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(tone_pos_idx, ch, :, tidx), 4));   % -> [n_trials x n_freqs]

    % FEEL & Neu
    pow_toneneu.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(tone_neu_idx, ch, :, tidx), 4));

     % FEEL & Neu
    pow_toneneg.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(tone_neg_idx, ch, :, tidx), 4));

end

%% Extract trial-level power per condition, per channel (for stats)

toi_stats = [0 3];   % adjust based on your hypothesis (e.g. [0 1] for early effects)
tidx = freq_clean_bl.time >= toi_stats(1) & freq_clean_bl.time <= toi_stats(2);

freq_axis = freq_clean_bl.freq; 


for c = 1:numel(plot_layout)
    chan_name = plot_layout{c};
    ch = find(strcmp(freq_clean_bl.label, chan_name), 1);

    % powspctrm dims: [trial x channel x freq x time]

    % FEEL & Pos
    pow_feelpos.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(feel_pos_idx, ch, :, tidx), 4));   % -> [n_trials x n_freqs]

    % FEEL & Neu
    pow_feelneu.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(feel_neu_idx, ch, :, tidx), 4));

     % FEEL & Neu
    pow_feelneg.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(feel_neg_idx, ch, :, tidx), 4));

end
%%

toi_stats = [0 2];   % adjust based on your hypothesis (e.g. [0 1] for early effects)
tidx = freq_clean_bl.time >= toi_stats(1) & freq_clean_bl.time <= toi_stats(2);

freq_axis = freq_clean_bl.freq; 

for c = 1:numel(plot_layout)
    chan_name = plot_layout{c};
    ch = find(strcmp(freq_clean_bl.label, chan_name), 1);

    % powspctrm dims: [trial x channel x freq x time]

    % Pos
    pow_pos.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(pos_idx, ch, :, tidx), 4));   % -> [n_trials x n_freqs]

    % Neu
    pow_neu.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(neu_idx, ch, :, tidx), 4));

     % Neu
    pow_neg.(chan_name) = squeeze(nanmean( ...
        freq_clean_bl.powspctrm(neg_idx, ch, :, tidx), 4));

end