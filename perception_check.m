%% =========================================================================
% VISUALIZE ALL PERCEPTION (FEEL) TRIALS
% =========================================================================
%% =========================================================================
% FULL PIPELINE — PSD: REST vs FEEL vs TONE
% =========================================================================

clear; clc; close all;

addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

%% =========================================================================
% LOAD TASK LFP (ALREADY EPOCHED — 96 trials)
% =========================================================================
load('lfp_firsthalf.mat');   % contains data_percept_epoched (96 trials)

data_all = data_percept_epoched;
nTrials  = length(data_all.trial);

fprintf('Loaded task LFP: %d trials\n', nTrials);

%% =========================================================================
% LOAD CONDITION CSV (FEEL / TONE)
% =========================================================================
T = readtable('C:\Users\fkamdar\Desktop\repos\cp_nonmotor\stimuli\blocks\cpeeg02\cpeeg02_block_96.csv');

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

time_ax = data_feel.time{1};
nTrials = length(data_feel.trial);
nChannels = length(data_feel.label);

figure('Color','w','Position',[100 100 1400 900]);

for ch = 1:nChannels

    subplot(ceil(nChannels/2),2,ch); hold on;

    for tr = 1:nTrials
        signal = data_feel.trial{tr}(ch,:);
        
        plot(time_ax, signal, ...
            'Color',[0 0 1 0.2], ...   % transparent blue
            'LineWidth',0.8);
    end

    % mark event time
    xline(0,'k--','LineWidth',1.2);

    xlabel('Time (s)');
    ylabel('Amplitude');
    title(data_feel.label{ch}, 'Interpreter','none');

    grid on; box off;

end

sgtitle(sprintf('All PERCEPTION trials (%d)', nTrials), ...
    'FontSize',13,'FontWeight','bold');