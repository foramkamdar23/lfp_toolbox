%% =========================================================================
% FULL PIPELINE — PSD: REST vs FEEL vs TONE
% =========================================================================

clear; clc; close all;

addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

%% =========================================================================
% LOAD TASK LFP (ALREADY EPOCHED — 96 trials)
% =========================================================================
% --- Load first half ---
load('lfp_firsthalf.mat');   % gives data_percept_epoched
data_first = data_percept_epoched;

% --- Load second half ---
load('lfp_secondhalf.mat');  % same variable name
data_second = data_percept_epoched;

% --- Combine ---
data_all = data_first;

data_all.trial = [data_first.trial, data_second.trial];
data_all.time  = [data_first.time,  data_second.time];

% Combine trialinfo if present
if isfield(data_first, 'trialinfo') && isfield(data_second, 'trialinfo')
    data_all.trialinfo = [data_first.trialinfo; data_second.trialinfo];
end

% Update sampleinfo if present (optional but good practice)
if isfield(data_first, 'sampleinfo') && isfield(data_second, 'sampleinfo')
    data_all.sampleinfo = [data_first.sampleinfo; data_second.sampleinfo];
end

% --- Summary ---
nTrials = length(data_all.trial);
fprintf('Combined LFP: %d trials\n', nTrials);
%% =========================================================================
% LOAD CONDITION CSV (FEEL / TONE)
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
    else
        error('Unknown condition');
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
% LOAD REST PERCEPT DATA
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

fprintf('REST loaded: %d samples\n', signal_length);

%% =========================================================================
% CREATE CLEAN REST EPOCHS (MANUAL GAP — NO TENS)
% =========================================================================

% >>> FROM YOUR PLOT <<<
clean_start = 20000;   % adjust slightly if needed
clean_end   = 70000;   % adjust slightly if needed

fprintf('Using manual clean window: %.1f – %.1f sec\n', ...
    clean_start/Fs, clean_end/Fs);

PRE = 0.5;
POST = 0.5;

pre_samples  = round(PRE*Fs);
post_samples = round(POST*Fs);
epoch_len    = pre_samples + post_samples;

trl = [];
step = round(epoch_len/2);

start_idx = clean_start;
trial_id = 1;

while (start_idx + epoch_len - 1) <= clean_end

    beg = start_idx;
    en  = start_idx + epoch_len - 1;

    trl(end+1,:) = [beg en -pre_samples trial_id]; %#ok<AGROW>

    start_idx = start_idx + step;
    trial_id  = trial_id + 1;
end

cfg = [];
cfg.trl = trl;

data_rest = ft_redefinetrial(cfg, data_rest_ft);

fprintf('Clean REST epochs: %d\n', length(data_rest.trial));

%% =========================================================================
% REMOVE NaNs FROM REST
% =========================================================================
nTrials_rest = length(data_rest.trial);
nan_mask = false(nTrials_rest,1);

for tr = 1:nTrials_rest
    if any(isnan(data_rest.trial{tr}(:)))
        nan_mask(tr) = true;
    end
end

cfg = [];
cfg.trials = find(~nan_mask);

data_rest = ft_selectdata(cfg, data_rest);

fprintf('Clean REST trials: %d\n', length(data_rest.trial));

%% =========================================================================
% PSD COMPUTATION
% =========================================================================
cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'pow';
cfg.taper      = 'dpss';
cfg.foilim     = [1 100];
cfg.tapsmofrq  = 1.5;
cfg.keeptrials = 'no';

freq_rest = ft_freqanalysis(cfg, data_rest);
freq_feel = ft_freqanalysis(cfg, data_feel);
freq_tone = ft_freqanalysis(cfg, data_tone);

fprintf('PSD done\n');

%% =========================================================================
% PLOT PSD
% =========================================================================
nChannels = length(freq_rest.label);
figure()
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
    ylabel('Power');

    title(freq_rest.label{ch}, 'Interpreter','none');

    xlim([1 100]);

    legend({'Rest','Perception','Regulation'}, ...
        'Location','northeast','Box','off');

    grid on;

end

sgtitle('Power Spectra — REST vs FEEL vs TONE');



%%
%% =========================================================================
% VALENCE-BASED SPLIT (NEG / NEU / POS) + PSD
% =========================================================================

% --- Extract valence ---
valence = T.Valence;

% Match cleaned trials (IMPORTANT)
valence_clean = valence(clean_idx);

% --- Define bins ---
neg_idx = find(valence_clean >= 1 & valence_clean <= 4);
neu_idx = find(valence_clean >= 4 & valence_clean <= 6);
pos_idx = find(valence_clean >= 6 & valence_clean <= 9);

fprintf('NEG: %d | NEU: %d | POS: %d\n', ...
    length(neg_idx), length(neu_idx), length(pos_idx));

% --- Select data ---
cfg = [];

cfg.trials = neg_idx;
data_neg = ft_selectdata(cfg, data_clean);

cfg.trials = neu_idx;
data_neu = ft_selectdata(cfg, data_clean);

cfg.trials = pos_idx;
data_pos = ft_selectdata(cfg, data_clean);

%% =========================================================================
% PSD FOR VALENCE CONDITIONS
% =========================================================================

cfg = [];
cfg.method     = 'mtmfft';
cfg.output     = 'pow';
cfg.taper      = 'dpss';
cfg.foilim     = [1 100];
cfg.tapsmofrq  = 1.5;
cfg.keeptrials = 'no';

freq_neg = ft_freqanalysis(cfg, data_neg);
freq_neu = ft_freqanalysis(cfg, data_neu);
freq_pos = ft_freqanalysis(cfg, data_pos);

fprintf('Valence PSD done\n');

%% =========================================================================
% PLOT VALENCE PSD
% =========================================================================

nChannels = length(freq_neg.label);

figure()
figure('Color','w','Position',[100 100 1200 800]);

colors.neg = [0.85 0.10 0.10];
colors.neu = [0.50 0.50 0.50];
colors.pos = [0.10 0.70 0.10];

for ch = 1:nChannels

    subplot(ceil(nChannels/2),2,ch); hold on;

    plot(freq_neg.freq, squeeze(freq_neg.powspctrm(ch,:)), ...
        'Color', colors.neg, 'LineWidth', 1.8);

    plot(freq_neu.freq, squeeze(freq_neu.powspctrm(ch,:)), ...
        'Color', colors.neu, 'LineWidth', 1.8);

    plot(freq_pos.freq, squeeze(freq_pos.powspctrm(ch,:)), ...
        'Color', colors.pos, 'LineWidth', 1.8);

    set(gca,'YScale','log','FontSize',12);

    xlabel('Frequency (Hz)');
    ylabel('Power');

    title(freq_neg.label{ch}, 'Interpreter','none');

    xlim([1 100]);

    legend({'Negative','Neutral','Positive'}, ...
        'Location','northeast','Box','off');

    grid on;

end

sgtitle('Power Spectra — NEG vs NEU vs POS');