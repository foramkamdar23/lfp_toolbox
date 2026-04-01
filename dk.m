%% =========================================================================
% =========================================================================
clear; clc; close all;

addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

%% =========================================================================
% LOAD PERCEPT DATA
% =========================================================================
load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124017_firsthalf_1437.mat");

new_datastruct = data.IndefiniteStreaming;
fs = new_datastruct(1).SampleRateInHz;

processed_data = struct();

for i = 1:length(new_datastruct)

    [signal_with_nan, ~] = detect_packet_loss_percept( ...
        new_datastruct(i).TimeDomainData, ...
        new_datastruct(i).GlobalSequences, ...
        new_datastruct(i).GlobalPacketSizes, ...
        fs);

    processed_data(i).channel = new_datastruct(i).Channel;
    processed_data(i).signal  = signal_with_nan;
    processed_data(i).fs      = fs;
end

%% =========================================================================
% BUILD FT STRUCTURE
% =========================================================================
nChannels = length(processed_data);
signal_length = length(processed_data(1).signal);

percept_all = nan(nChannels, signal_length);

for ch = 1:nChannels
    percept_all(ch,:) = processed_data(ch).signal;
end

t = (0:signal_length-1)/fs;

data_percept_ft = [];
data_percept_ft.label = {processed_data.channel}';
data_percept_ft.fsample = fs;
data_percept_ft.trial = {percept_all};
data_percept_ft.time = {t};
data_percept_ft.sampleinfo = [1 signal_length];

data_percept_ft = ft_checkdata(data_percept_ft,'datatype','raw');

%% =========================================================================
% LOAD EEG EVENTS
% =========================================================================
dataDir = 'C:\Users\fkamdar\Desktop\repos\data\eeg-selected';
eegFile = 'cpeeg02_b01_ert.bdf';
eegFullPath = fullfile(dataDir, eegFile);

events_raw    = ft_read_event(eegFullPath);
status_events = events_raw(strcmp({events_raw.type}, 'STATUS'));

if ischar(status_events(1).value) || isstring(status_events(1).value)
    raw_values = str2double({status_events.value});
else
    raw_values = double([status_events.value]);
end

corrected_values = bitand(uint16(raw_values) - 512, 255);

%% =========================================================================
% SETTINGS
% =========================================================================
EPOCH_CODE = 30;
POST_STIM_S = 1.0;
BASELINE_S  = 0.5;

Fs_percept = fs;
Fs_bdf     = 2048;

post_samples     = round(POST_STIM_S * Fs_percept);
baseline_samples = round(BASELINE_S  * Fs_percept);

%% =========================================================================
% EVENT SAMPLES
% =========================================================================
epoch_events  = status_events(corrected_values == EPOCH_CODE);
epoch_samples = [epoch_events.sample];

epoch_samples_ds = round(epoch_samples * (Fs_percept / Fs_bdf));

%% =========================================================================
% BUILD TRIALS (0–1 sec)
% =========================================================================
trl = [];

for k = 1:length(epoch_samples_ds)
    beg = epoch_samples_ds(k);
    en  = epoch_samples_ds(k) + post_samples;
    trl(end+1,:) = [beg en 0 30]; %#ok<AGROW>
end

% remove invalid
nSamples = size(data_percept_ft.trial{1},2);
valid = trl(:,1) > baseline_samples & trl(:,2) <= nSamples;

trl = trl(valid,:);
epoch_samples_ds = epoch_samples_ds(valid);

%% =========================================================================
% BASELINE EXTRACTION
% =========================================================================
baseline_values = nan(size(trl,1), nChannels);

for k = 1:size(trl,1)
    base_beg = epoch_samples_ds(k) - baseline_samples;
    base_end = epoch_samples_ds(k) - 1;

    seg = data_percept_ft.trial{1}(:, base_beg:base_end);
    baseline_values(k,:) = nanmean(seg,2);
end

%% =========================================================================
% EPOCH DATA
% =========================================================================
cfg = [];
cfg.trl = trl;
data_epoched = ft_redefinetrial(cfg, data_percept_ft);

%% =========================================================================
% BASELINE NORMALIZATION
% =========================================================================
for tr = 1:length(data_epoched.trial)
    for ch = 1:nChannels
        baseline = baseline_values(tr,ch);
        if ~isnan(baseline)
            data_epoched.trial{tr}(ch,:) = ...
                data_epoched.trial{tr}(ch,:) - baseline;
        end
    end
end

%% =========================================================================
% COMPUTE TENS ENVELOPE (80 Hz)
% =========================================================================
signal = processed_data(1).signal;
signal(isnan(signal)) = interp1(find(~isnan(signal)), signal(~isnan(signal)), find(isnan(signal)),'linear','extrap');

[b,a] = butter(4,[75 85]/(Fs_percept/2),'bandpass');
signal_filt = filtfilt(b,a,signal);
signal_percept_env = zscore(abs(signal_filt));

%% =========================================================================
% REMOVE TENS TRIALS
% =========================================================================
THRESH_Z = 2.5;
THRESH_FRAC = 0.1;

good_trials = true(length(data_epoched.trial),1);

for tr = 1:length(data_epoched.trial)

    beg = trl(tr,1);
    en  = trl(tr,2);

    seg = signal_percept_env(beg:en);

    if mean(seg > THRESH_Z) > THRESH_FRAC
        good_trials(tr) = false;
    end
end

removed_trials = find(~good_trials);

fprintf('\n--- TENS REMOVAL ---\n');
fprintf('Removed %d trials\n', length(removed_trials));
disp('Removed trial IDs:');
disp(removed_trials');

cfg = [];
cfg.trials = find(good_trials);

data_clean = ft_selectdata(cfg, data_epoched);

%% =========================================================================
% KEEP TRIALS 2–48
% =========================================================================
data_out = data_clean;

startIdx = 2;
endIdx   = min(48, length(data_out.trial));

data_out.trial = data_out.trial(startIdx:endIdx);
data_out.time  = data_out.time(startIdx:endIdx);

fprintf('Saved trials %d → %d\n', startIdx, endIdx);

%% =========================================================================
% SAVE
% =========================================================================
save('lfp_firsthalf_clean.mat','data_out');

%% =========================================================================
% QUICK CHECK
% =========================================================================
figure;
plot(data_out.time{1}, data_out.trial{1}(1,:));
title('Example Clean Trial');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;