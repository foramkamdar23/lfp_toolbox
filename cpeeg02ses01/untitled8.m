%% =========================================================================
% FULL PIPELINE — PSD: REST vs FEEL vs TONE
% =========================================================================

clear; clc; close all;

addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

%% =========================================================================
% LOAD TASK LFP
% =========================================================================
load('lfp_firsthalf.mat');
data_first = data_percept_epoched;

load('lfp_secondhalf.mat');
data_second = data_percept_epoched;

data_all = data_first;

data_all.trial = [data_first.trial, data_second.trial];
data_all.time  = [data_first.time,  data_second.time];

if isfield(data_first, 'trialinfo') && isfield(data_second, 'trialinfo')
    data_all.trialinfo = [data_first.trialinfo; data_second.trialinfo];
end

if isfield(data_first, 'sampleinfo') && isfield(data_second, 'sampleinfo')
    data_all.sampleinfo = [data_first.sampleinfo; data_second.sampleinfo];
end

nTrials = length(data_all.trial);
fprintf('Combined LFP: %d trials\n', nTrials);

%% =========================================================================
% LOAD CONDITION CSV
% =========================================================================
T = readtable('C:\Users\fkamdar\Desktop\repos\cp_nonmotor\stimuli\blocks\cpeeg02\cpeeg02_block_96.csv');
T(1,:) = [];

if height(T) ~= nTrials
    error('Mismatch between CSV and trials');
end

condition = zeros(nTrials,1);

for i = 1:nTrials
    if strcmp(T.condition{i}, 'FEEL')
        condition(i) = 1;
    elseif strcmp(T.condition{i}, 'TONE')
        condition(i) = 2;
    end
end

%% =========================================================================
% REMOVE NaN TRIALS
% =========================================================================
nan_any = false(nTrials,1);

for tr = 1:nTrials
    if any(isnan(data_all.trial{tr}(:)))
        nan_any(tr) = true;
    end
end

clean_idx = find(~nan_any);

cfg = [];
cfg.trials = clean_idx;

data_clean = ft_selectdata(cfg, data_all);
condition_clean = condition(clean_idx);

fprintf('Clean trials: %d\n', length(clean_idx));

%% =========================================================================
% %%% NEW: KEEP ONLY TRIALS WITH FULL [-0.5 1] WINDOW
% =========================================================================
baseline_win = [-0.5 0];
active_win   = [0 1];

good_trials = false(length(data_clean.trial),1);

for tr = 1:length(data_clean.trial)
    t = data_clean.time{tr};

    if min(t) <= baseline_win(1) && max(t) >= active_win(2)
        good_trials(tr) = true;
    end
end

fprintf('Keeping %d / %d trials with full window\n', ...
    sum(good_trials), length(good_trials));

cfg = [];
cfg.trials = find(good_trials);

data_clean = ft_selectdata(cfg, data_clean);
condition_clean = condition_clean(good_trials);

%% =========================================================================
% SPLIT FEEL vs TONE
% =========================================================================
cfg = [];
cfg.trials = find(condition_clean == 1);
data_feel = ft_selectdata(cfg, data_clean);

cfg.trials = find(condition_clean == 2);
data_tone = ft_selectdata(cfg, data_clean);

fprintf('FEEL: %d | TONE: %d\n', ...
    length(data_feel.trial), length(data_tone.trial));

%% =========================================================================
% %%% NEW: SPLIT BASELINE + ACTIVE WINDOWS
% =========================================================================
cfg = [];
cfg.latency = baseline_win;
data_feel_base = ft_selectdata(cfg, data_feel);
data_tone_base = ft_selectdata(cfg, data_tone);

cfg.latency = active_win;
data_feel_act = ft_selectdata(cfg, data_feel);
data_tone_act = ft_selectdata(cfg, data_tone);

%% =========================================================================
% LOAD REST (UNCHANGED)
% =========================================================================
load("C:\Users\fkamdar\Desktop\repos\WorkingMemory_CPEEG\cpeeg02\cpeeg02_b01_offmedoffstim_rest_percept.mat");

nChannels = length(processed_data);
Fs = processed_data(1).fs;
signal_length = length(processed_data(1).signal);

percept_all = nan(nChannels, signal_length);

for ch = 1:nChannels
    percept_all(ch,:) = processed_data(ch).signal;
end

t = (0:signal_length-1)/Fs;

data_rest_ft = [];
data_rest_ft.label = {processed_data.channel}';
data_rest_ft.fsample = Fs;
data_rest_ft.trial = {percept_all};
data_rest_ft.time = {t};
data_rest_ft.sampleinfo = [1 signal_length];

data_rest_ft = ft_checkdata(data_rest_ft, 'datatype','raw');

%% =========================================================================
% REST EPOCHING (UNCHANGED)
% =========================================================================
clean_start = 20000;
clean_end   = 70000;

PRE = 0.5;
POST = 0.5;

pre_samples  = round(PRE*Fs);
post_samples = round(POST*Fs);
epoch_len    = pre_samples + post_samples;

trl = [];
step = round(epoch_len/2);
start_idx = clean_start;

while (start_idx + epoch_len - 1) <= clean_end
    beg = start_idx;
    en  = start_idx + epoch_len - 1;
    trl(end+1,:) = [beg en -pre_samples]; %#ok<AGROW>
    start_idx = start_idx + step;
end

cfg = [];
cfg.trl = trl;

data_rest = ft_redefinetrial(cfg, data_rest_ft);

%% =========================================================================
% PSD CONFIG
% =========================================================================
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'pow';
cfg.taper      = 'dpss';
cfg.foilim     = [1 100];
cfg.tapsmofrq  = 1.5;
cfg.keeptrials = 'no';

%% =========================================================================
% %%% NEW: COMPUTE PSD WITH BASELINE NORMALIZATION
% =========================================================================
freq_rest = ft_freqanalysis(cfg, data_rest);

freq_feel_base = ft_freqanalysis(cfg, data_feel_base);
freq_tone_base = ft_freqanalysis(cfg, data_tone_base);

freq_feel_act  = ft_freqanalysis(cfg, data_feel_act);
freq_tone_act  = ft_freqanalysis(cfg, data_tone_act);

% Normalize
freq_feel = freq_feel_act;
freq_tone = freq_tone_act;

freq_feel.powspctrm = ...
    (freq_feel_act.powspctrm - freq_feel_base.powspctrm) ./ ...
    (freq_feel_base.powspctrm + 1e-10);

freq_tone.powspctrm = ...
    (freq_tone_act.powspctrm - freq_tone_base.powspctrm) ./ ...
    (freq_tone_base.powspctrm + 1e-10);

fprintf('PSD done\n');

%% =========================================================================
% PLOT (UNCHANGED except ylabel)
% =========================================================================
nChannels = length(freq_rest.label);

figure('Color','w','Position',[100 100 1200 800]);

colors.rest = [0.00 0.62 0.65];
colors.feel = [0.00 0.00 0.90];
colors.tone = [0.90 0.20 0.20];

for ch = 1:nChannels

    subplot(ceil(nChannels/2),2,ch); hold on;

    plot(freq_rest.freq, squeeze(freq_rest.powspctrm(ch,:)), ...
        'Color', colors.rest, 'LineWidth', 1.8);

    plot(freq_feel.freq, squeeze(freq_feel.powspctrm(ch,:)), ...
        'Color', colors.feel, 'LineWidth', 1.8);

    plot(freq_tone.freq, squeeze(freq_tone.powspctrm(ch,:)), ...
        'Color', colors.tone, 'LineWidth', 1.8);

    set(gca,'YScale','log','FontSize',12);

    xlabel('Frequency (Hz)');
    ylabel('Normalized Power');

    title(freq_rest.label{ch}, 'Interpreter','none');

    xlim([1 100]);

    legend({'Rest','Perception','Regulation'}, ...
        'Location','northeast','Box','off');

    grid on;

end

sgtitle('PSD (0–1s, baseline -0.5–0)');