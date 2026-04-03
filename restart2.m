%CPEEG02
clear all;
close all;
clc;

%% Cell 1
data1 = load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124017_firsthalf_1437.mat");
data2 = load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124034_secondhalf_1453.mat");
first_half = data1.data.IndefiniteStreaming;
second_half = data2.data.IndefiniteStreaming;

%% Cell 2 load lfp
ch_id = 1;
first_half_ch1 = first_half(ch_id);
second_half_ch1 = second_half(ch_id);


sig1 = first_half_ch1.TimeDomainData;
sig2 = second_half_ch1.TimeDomainData;
fs_lfp = first_half_ch1.SampleRateInHz;

t_lfp1 = (0:length(sig1)-1)/fs_lfp;
t_lfp2  = (0:length(sig2)-1)/fs_lfp;


figure;
subplot(2,1,1); plot(t_lfp1, sig1); title('LFP1');
subplot(2,1,2); plot(t_lfp2, sig2); title('LFP2');


% load EX 7 to align
% fieldtrip 
addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

cfg = [];
cfg.dataset = 'C:\Users\fkamdar\Desktop\repos\data\eeg-selected\cpeeg02_b01_ert.bdf';
%cfg.channel = {'EXG7', 'Status'};
cfg.channel = {'EXG7'};

exg7_data = ft_preprocessing(cfg);

exg7_eeg = exg7_data.trial{1};
t_eeg = exg7_data.time{1};
fs_eeg = exg7_data.fsample;

figure;
plot(t_eeg, exg7_eeg);
title('EXG7');

% down sampling
exg7_ds = resample(exg7_eeg, fs_lfp, fs_eeg);
t_exg7_ds = (0:length(exg7_ds)-1)/fs_lfp;

%qc
figure;
subplot(2,1,1); plot(t_eeg, exg7_eeg); title('Original EEG');
subplot(2,1,2); plot(t_exg7_ds, exg7_ds); title('Downsampled EEG (250 Hz)');


% filter for enevelope around 80hz tens burst
[b,a] = butter(4, [75 85]/(fs_lfp/2), 'bandpass');

exg7_ds_filt = filtfilt(b,a, exg7_ds);
lfp1_filt = filtfilt(b,a, sig1);
lfp2_filt = filtfilt(b,a, sig2);


figure;
subplot(3,1,1); plot(exg7_ds_filt); title('EEG');
subplot(3,1,2); plot(lfp1_filt); title('LFP1');
subplot(3,1,3); plot(lfp2_filt); title('LFP2');


%env created and normalized
exg7_env = abs(hilbert(exg7_ds_filt));
lfp1_env = abs(hilbert(lfp1_filt));
lfp2_env = abs(hilbert(lfp2_filt));

exg7_env = zscore(exg7_env); %(exg7_env - mean(exg7_env))/std(exg7_env); 
lfp1_env = zscore(lfp1_env); %(lfp_env - mean(lfp_env))/std(lfp_env);
lfp2_env = zscore(lfp2_env); 

%cross correlation
[xc1, lags1] = xcorr(exg7_env, lfp1_env);
[~, idx1] = max(xc1); % find peak
lag_samples1 = lags1(idx1);

[xc2, lags2] = xcorr(exg7_env, lfp2_env);
[~, idx2] = max(xc2); % find peak
lag_samples2 = lags2(idx2);


% convert to seconds
lag_sec1 = lag_samples1 / fs_lfp;
lag_sec2 = lag_samples2 / fs_lfp;


%align here
t_lfp1_aligned = t_lfp1 + lag_sec1;
t_lfp2_aligned = t_lfp2 + lag_sec2;

figure;
plot(t_exg7_ds, exg7_env); hold on;
plot(t_lfp1_aligned, lfp1_env); hold on;
plot(t_lfp2_aligned, lfp2_env);
legend('EEG ENV','LFP ENV (aligned)');

fprintf("done with alinging")

%% Combine LFP1 & LPF2

t_master = t_exg7_ds;
lfp_merged = nan(size(t_master)); % master timeline


% mapping lfp1 to master timeline
sample_idx1 = round(t_lfp1_aligned * fs_lfp) + 1;
valid1 = sample_idx1 > 0 & sample_idx1 <= length(lfp_merged);
lfp_merged(sample_idx1(valid1)) = sig1(valid1);

sample_idx2 = round(t_lfp2_aligned * fs_lfp) + 1;
valid2 = sample_idx2 > 0 & sample_idx2 <= length(lfp_merged);
lfp_merged(sample_idx2(valid2)) = sig2(valid2);

figure;
plot(t_master, exg7_ds_filt); hold on;
plot(t_master, lfp_merged);
legend('EEG','Merged LFP');


%plot nans as gray shaded areafigure; hold on;

% plot LFP (MATLAB will break at NaNs automatically)
plot(t_master, lfp_merged, 'b');

% highlight NaNs as gray background
nan_mask = isnan(lfp_merged);
d = diff([0 nan_mask 0]);
start_idx = find(d == 1);
end_idx   = find(d == -1) - 1;

y_limits = ylim;

for i = 1:length(start_idx)
    patch([t_master(start_idx(i)) t_master(end_idx(i)) t_master(end_idx(i)) t_master(start_idx(i))], ...
          [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
          [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

legend('LFP','NaN regions');
title('Merged LFP with NaN gaps');

%% Biosemi Triggers

% Lets look into the triggers now
events = ft_read_event('C:\Users\fkamdar\Desktop\repos\data\eeg-selected\cpeeg02_b01_ert.bdf');
events(1:3) = [];
trig = [events.value]'; 
trig = trig-512;
trig = double(bitand(int32(trig),255)); % we have triggers now

% lets do time of the triggers
trig_samples = [events.sample]';
trig_time = trig_samples / exg7_data.hdr.Fs;

trig_table = table(trig,trig_time, 'VariableNames', {'TriggerCode', 'TriggerTime'});
trig_table = trig_table(trig_table.TriggerCode ~= 0, :);

trig_map = containers.Map( ...
    [1 2 11 20 12 30 13 40 99 200 201 202], ...
    {'intro1','intro2','fix1','cue','fix2','image','fix3','valence_onset','end_screen','break_start','resync_start','resync_end'} ...
);
labels = strings(height(trig_table),1);

for i = 1:height(trig_table)
    code = trig_table.TriggerCode(i);
    
    if isKey(trig_map, code)
        labels(i) = trig_map(code);
        
    elseif code >= 51 && code <= 59   % valence responses
        rating = code - 50;
        labels(i) = "valence_resp_" + string(rating);
        
    else
        labels(i) = "unknown";
    end
end

trig_table.Event = labels;



figure();
plot(t_master, lfp_merged); hold on;

% filter only fix1 events
fix1_idx = trig_table.Event == "fix1";
fix1_times = trig_table.TriggerTime(fix1_idx);

% plot only fix1 triggers
for i = 1:length(fix1_times)
    xline(fix1_times(i), 'r');
end

title('LFP + fix1 Triggers');
ylabel('LFP');

%% MERGE ALL CHANNELS
n_ch = 6;
lfp_merged_all = nan(n_ch, length(t_master));
lfp1_all = data1.data.IndefiniteStreaming;
lfp2_all = data2.data.IndefiniteStreaming;


for ch = 1:n_ch
    
    sig1 = lfp1_all(ch).TimeDomainData;
    sig2 = lfp2_all(ch).TimeDomainData;
    
    % align time
    t1_aligned = t_lfp1 + lag_sec1;
    t2_aligned = t_lfp2 + lag_sec2;
    
    % indices
    idx1 = round(t1_aligned * fs_lfp) + 1;
    idx2 = round(t2_aligned * fs_lfp) + 1;
    
    % valid masks
    valid1 = idx1 > 0 & idx1 <= length(t_master);
    valid2 = idx2 > 0 & idx2 <= length(t_master);
    
    % insert into merged matrix
    lfp_merged_all(ch, idx1(valid1)) = sig1(valid1);
    lfp_merged_all(ch, idx2(valid2)) = sig2(valid2);
end

%lfp_merged_all


%% EPOCHING

%  Convert LFP to fieldtrip

% channel labels (simple for now)
channel_names = cell(n_ch,1);
for ch = 1:n_ch
    channel_names{ch} = first_half(ch).Channel;
end

data_ft = [];
data_ft.label   = channel_names;
data_ft.fsample = fs_lfp;

% continuous data
data_ft.trial = {lfp_merged_all};     % [channels x time]
data_ft.time  = {t_master};

% required bookkeeping
data_ft.sampleinfo = [1 length(t_master)];

% header (important for some FT functions)
data_ft.hdr = [];
data_ft.hdr.label     = channel_names;
data_ft.hdr.Fs        = fs_lfp;
data_ft.hdr.nChans    = n_ch;
data_ft.hdr.nSamples  = length(t_master);
data_ft.hdr.nTrials   = 1;

% validate
data_ft = ft_checkdata(data_ft, 'datatype','raw','hassampleinfo','yes');

fprintf("FieldTrip structure created\n");


% Define trials
%% ================= DEFINE TRIALS =================

image_idx = trig_table.Event == "image";
image_samples = round(trig_table.TriggerTime(image_idx) * fs_lfp);

pre  = round(4 * fs_lfp);
post = round(5 * fs_lfp);

trl = [];

for i = 1:length(image_samples)
    
    s = image_samples(i);
    
    start_sample = s - pre;
    end_sample   = s + post;
    offset       = -pre;
    
    if start_sample > 0 && end_sample <= length(t_master)
        trl = [trl; start_sample end_sample offset];
    end
end

%
cfg = [];
cfg.trl = trl;

data_epoched = ft_redefinetrial(cfg, data_ft);

fprintf("Epoching complete\n");

%% Reject bad trials 

% check for TENS contamination:

[b,a] = butter(4, [75 85]/(fs_lfp/2), 'bandpass');

n_trials = length(data_epoched.trial);
tens_bad = false(n_trials,1);
power_80hz = zeros(n_trials,1);

for i = 1:n_trials
    
    trial_data = data_epoched.trial{i};   % [channels x time]
    
    % --- filter ---
    trial_filt = filtfilt(b,a, trial_data')';  % filter per channel
    
    % --- envelope ---
    env = abs(hilbert(trial_filt')');  % [channels x time]
    
    % --- collapse across channels ---
    env_mean = mean(env,1);
    
    % --- quantify energy ---
    power_80hz(i) = mean(env_mean);
    
end

% --- threshold (adaptive) ---
threshold = mean(power_80hz) + 2*std(power_80hz);

tens_bad = power_80hz > threshold;

bad_trials = find(tens_bad);


% ================= REMOVE NaN TRIALS =================

nan_bad = false(n_trials,1);

for i = 1:n_trials
    trial_data = data_epoched.trial{i};
    
    if any(isnan(trial_data), 'all')
        nan_bad(i) = true;
    end
end

nan_trials = find(nan_bad);

fprintf('NaN trials (%d):\n', length(nan_trials));
disp(nan_trials);


% ALL bad:
bad_all = tens_bad | nan_bad;

bad_trials_all = find(bad_all);

fprintf('Total bad trials (%d):\n', length(bad_trials_all));
disp(bad_trials_all);


%================= TRACKING STRUCT =================

trial_info = table( ...
    (1:n_trials)', ...
    tens_bad, ...
    nan_bad, ...
    bad_all, ...
    'VariableNames', {'Trial','TENS_Bad','NaN_Bad','Bad_Overall'} ...
);

disp(trial_info);

% clean data
%% ================= REMOVE BAD TRIALS =================

good_idx = find(~bad_all);

cfg = [];
cfg.trials = good_idx;

data_clean = ft_selectdata(cfg, data_epoched);



%% -------- CLEAN TABLE --------
trial_info_clean = trial_info(~bad_all, :);

fprintf('Labeling complete\n');

%
figure;

n_ch = length(data_clean.label);

for ch = 1:n_ch
    
    subplot(3,2,ch); hold on;
    
    % --- plot all trials ---
    for i = 1:length(data_clean.trial)
        plot(data_clean.time{i}, data_clean.trial{i}(ch,:), ...
             'Color', [0.8 0.8 0.8]);
    end
    
    % --- compute mean ---
    n_trials = length(data_clean.trial);
    n_time = length(data_clean.time{1});
    
    all_trials = zeros(n_trials, n_time);
    
    for i = 1:n_trials
        all_trials(i,:) = data_clean.trial{i}(ch,:);
    end
    
    mean_sig = mean(all_trials, 1);
    
    % --- plot mean ---
    plot(data_clean.time{1}, mean_sig, 'k', 'LineWidth', 2);
    
    xline(0,'r--');
    
    title(data_clean.label{ch});
    xlabel('Time (s)');
    ylabel('LFP');
end

% plot withou

%% Baseline correction:
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.5 0];

data_bl = ft_preprocessing(cfg, data_clean);

% Plot baseline plot
figure;

n_ch = length(data_bl.label);

for ch = 1:n_ch
    
    subplot(3,2,ch); hold on;
    
    % --- plot all trials ---
    for i = 1:length(data_bl.trial)
        plot(data_bl.time{i}, data_bl.trial{i}(ch,:), ...
             'Color', [0.8 0.8 0.8]);
    end
    
    % --- compute mean ---
    n_trials = length(data_bl.trial);
    n_time = length(data_bl.time{1});
    
    all_trials = zeros(n_trials, n_time);
    
    for i = 1:n_trials
        all_trials(i,:) = data_bl.trial{i}(ch,:);
    end
    
    mean_sig = mean(all_trials, 1);
    
    % --- plot mean ---
    plot(data_bl.time{1}, mean_sig, 'k', 'LineWidth', 2);
    
    xline(0,'r--');
    
    title(data_bl.label{ch});
    xlabel('Time (s)');
    ylabel('LFP');
end

%%

% %% LABEL Trials
csv_data = readtable('C:\Users\fkamdar\Desktop\repos\lfp_toolbox\cpeeg02_block_96.csv');

trial_info.Condition = string(csv_data.condition);
trial_info_clean = trial_info(~bad_all, :);

% Now define condition indices
feel_idx = strcmp(trial_info_clean.Condition, 'FEEL');
tone_idx = strcmp(trial_info_clean.Condition, 'TONE');

cfg = [];
cfg.trials = find(feel_idx);
data_feel = ft_selectdata(cfg, data_clean);

cfg = [];
cfg.trials = find(tone_idx);
data_tone = ft_selectdata(cfg, data_clean);

% VALENCE (NEG / NEU / POS) 
trial_info.Valence = string(csv_data.Valence_group);
trial_info_clean = trial_info(~bad_all, :);

pos_idx = strcmp(trial_info_clean.Valence, 'Pos');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');

cfg = [];
cfg.trials = find(pos_idx);
data_pos = ft_selectdata(cfg, data_clean);

cfg = [];
cfg.trials = find(neu_idx);
data_neu = ft_selectdata(cfg, data_clean);

cfg = [];
cfg.trials = find(neg_idx);
data_neg = ft_selectdata(cfg, data_clean);

%% PLOT MEAN 
cfg = [];
avg_feel = ft_timelockanalysis(cfg, data_feel);

cfg = [];
avg_tone = ft_timelockanalysis(cfg, data_tone);

figure;

n_ch = length(avg_feel.label);

for ch = 1:n_ch
    
    subplot(3,2,ch); hold on;
    
    plot(avg_feel.time, avg_feel.avg(ch,:), 'b', 'LineWidth',2);
    plot(avg_tone.time, avg_tone.avg(ch,:), 'r', 'LineWidth',2);
    
    xline(0,'k--'); % image onset
    
    title(avg_feel.label{ch});
    xlabel('Time (s)');
    ylabel('LFP');
end

legend('FEEL','TONE');


%%
cfg = [];
avg_neg = ft_timelockanalysis(cfg, data_neg);

cfg = [];
avg_neu = ft_timelockanalysis(cfg, data_neu);

cfg = [];
avg_pos = ft_timelockanalysis(cfg, data_pos);

figure;

n_ch = length(avg_neg.label);

for ch = 1:n_ch
    
    subplot(3,2,ch); hold on;
    
    plot(avg_neg.time, avg_neg.avg(ch,:), 'r', 'LineWidth',2);
    plot(avg_neu.time, avg_neu.avg(ch,:), 'g', 'LineWidth',2);
    plot(avg_pos.time, avg_pos.avg(ch,:), 'b', 'LineWidth',2);
    
    xline(0,'k--'); % image onset
    
    title(avg_neg.label{ch});
    xlabel('Time (s)');
    ylabel('LFP');
end

legend('NEG','NEU','POS');

%% TFR
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

%avg
cfg = [];
cfg.avgoverrpt = 'yes';

freq_avg = ft_selectdata(cfg, freq_clean_bl);

n_ch = length(freq_avg.label);

%plot
figure;
for ch = 1:n_ch
    
    subplot(ceil(n_ch/2), 2, ch);
    
    imagesc(freq_avg.time, freq_avg.freq, ...
        squeeze(freq_avg.powspctrm(ch,:,:)), [-2 2]); 
    
    axis xy;
    title(freq_avg.label{ch});
    
    xline(0,'w--');
    colorbar;

    colormap jet
end

xlabel('Time (s)');
ylabel('Frequency (Hz)');


data = freq_avg.powspctrm;

neg_vals = data(data < 0);
length(neg_vals)
pos_vals = data(data > 0);
length(pos_vals)
