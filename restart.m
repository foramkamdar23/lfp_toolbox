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
figure;
subplot(2,1,1); plot(exg7_ds_filt); title('EEG 75–85 Hz');
subplot(2,1,2); plot(lfp1_filt); title('LFP1 75–85 Hz');

%env created and normalized
exg7_env = abs(hilbert(exg7_ds_filt));
lfp_env = abs(hilbert(lfp1_filt));
exg7_env = zscore(exg7_env); %(exg7_env - mean(exg7_env))/std(exg7_env); 
lfp_env = zscore(lfp_env); %(lfp_env - mean(lfp_env))/std(lfp_env);
n = min(length(exg7_env), length(lfp_env));
exg7_env = exg7_env(1:n);
lfp_env  = lfp_env(1:n);
t_exg7_ds = t_exg7_ds(1:n);


%cross correlation
max_lag_sec = 30; % search window
max_lag = fs_lfp * max_lag_sec;

[xc, lags] = xcorr(exg7_env, lfp_env, max_lag, 'coeff');

% find peak
[~, idx] = max(xc);
lag_samples = lags(idx);

% convert to seconds
lag_sec = lag_samples / fs_lfp;


%align here
t_lfp1_aligned = t_lfp1 + lag_sec;



figure;
plot(t_exg7_ds, exg7_env); hold on;
plot(t_lfp1_aligned, lfp_env);
legend('EEG ENV','LFP ENV (aligned)');
xlim([0 200]);
title('Envelope alignment');
xlim([t_lfp1_aligned(end)-100 t_lfp1_aligned(end)]);




























% filter at 80 hz
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [75 85];   % adjust if needed
cfg.channel = {'EXG7'};
exg7_data = ft_preprocessing(cfg, sync_data);



% Lets look into the triggers now
events = ft_read_event('C:\Users\fkamdar\Desktop\repos\data\eeg-selected\cpeeg02_b01_ert.bdf');
events(1:3) = [];
trig = [events.value]'; 
trig = trig-512;
trig = double(bitand(int32(trig),255)); % we have triggers now

% lets do time of the triggers
trig_samples = [events.sample]';
trig_time = trig_samples / exg7_data.hdr.Fs;

trig_table = table(trig,trig_time, 'VariableNames', {'TriggerCode', 'Trigger Time'});
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