clear;
close all;
fieldtrip_path = "C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507";

addpath(fieldtrip_path);
ft_defaults;

load("C:\Users\fkamdar\Desktop\repos\lfp_toolbox\derivatives\sub-CPEEG02_ses-01_data_clean.mat")

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
cfg.foi = 1:1:99;
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

%% TFR POS AND Instuction
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

%% TFR Neu AND Instuction
pos_idx = strcmp(trial_info_clean.Valence, 'Neu');

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
              sprintf('%s, %s, FEEL, Neu', subj_id, sess_cond), label_map);

plot_tfr_grid(freq_avg_tone_pos, plot_layout, tfr_clim, ...
              sprintf('%s, %s, TONE, Neu', subj_id, sess_cond), label_map);

%% TFR Neg AND Instuction
pos_idx = strcmp(trial_info_clean.Valence, 'Neg');

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
              sprintf('%s, %s, FEEL, Neg', subj_id, sess_cond), label_map);

plot_tfr_grid(freq_avg_tone_pos, plot_layout, tfr_clim, ...
              sprintf('%s, %s, TONE, Neg', subj_id, sess_cond), label_map);