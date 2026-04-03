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

% %% save data for python for later
% fprintf("Saving data for Python...\n");
% 
% % --- 1. Ensure correct shapes ---
% t_master = t_master(:);  % column vector
% 
% % --- 2. Trigger samples in LFP space ---
% trig_samples_lfp = round(trig_table.TriggerTime * fs_lfp);
% 
% % --- 3. Channel labels ---
% channel_names = strings(n_ch,1);
% for ch = 1:n_ch
%     channel_names(ch) = "ch_" + string(ch);
% end
% 
% % --- 4. Save everything ---
% save_path = "C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\cpeeg02_lfp_merged.mat";
% 
% save(save_path, ...
%     'lfp_merged_all', ...
%     't_master', ...
%     'fs_lfp', ...
%     'trig_table', ...
%     'trig_samples_lfp', ...
%     'channel_names', ...
%     '-v7.3');   % use v7.3 for large data
% 
% fprintf("Saved to: %s\n", save_path);