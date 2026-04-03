

% CPEEG02 MULTI-CHANNEL ALIGN + MERGE
clear; close all; clc;

%% ================= LOAD DATA =================
data1 = load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124017_firsthalf_1437.mat");
data2 = load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124034_secondhalf_1453.mat");

lfp1_all = data1.data.IndefiniteStreaming;
lfp2_all = data2.data.IndefiniteStreaming;

n_ch = length(lfp1_all);   % should be 6
fs_lfp = lfp1_all(1).SampleRateInHz;

%% ================= TIME VECTORS =================
sig1 = lfp1_all(1).TimeDomainData;
sig2 = lfp2_all(1).TimeDomainData;

t_lfp1 = (0:length(sig1)-1)/fs_lfp;
t_lfp2 = (0:length(sig2)-1)/fs_lfp;

%% ================= LOAD EEG (ALIGNMENT ONLY) =================
addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

cfg = [];
cfg.dataset = 'C:\Users\fkamdar\Desktop\repos\data\eeg-selected\cpeeg02_b01_ert.bdf';
cfg.channel = {'EXG7'};
exg7_data = ft_preprocessing(cfg);

exg7 = exg7_data.trial{1};
fs_eeg = exg7_data.fsample;

% downsample ONLY for alignment
exg7_ds = resample(exg7, fs_lfp, fs_eeg);
t_master = (0:length(exg7_ds)-1)/fs_lfp;

%% ================= ALIGNMENT (CHANNEL 1 ONLY) =================
[b,a] = butter(4, [75 85]/(fs_lfp/2), 'bandpass');

exg7_filt = filtfilt(b,a, exg7_ds);
lfp1_filt = filtfilt(b,a, sig1);
lfp2_filt = filtfilt(b,a, sig2);

exg7_env = zscore(abs(hilbert(exg7_filt)));
lfp1_env = zscore(abs(hilbert(lfp1_filt)));
lfp2_env = zscore(abs(hilbert(lfp2_filt)));

% cross-correlation
[xc1, lags1] = xcorr(exg7_env, lfp1_env);
[~, idx1] = max(xc1);
lag_sec1 = lags1(idx1)/fs_lfp;

[xc2, lags2] = xcorr(exg7_env, lfp2_env);
[~, idx2] = max(xc2);
lag_sec2 = lags2(idx2)/fs_lfp;

fprintf("Alignment done: LFP1=%.3fs | LFP2=%.3fs\n", lag_sec1, lag_sec2);

%% ================= MULTI-CHANNEL MERGE =================
lfp_merged_all = nan(n_ch, length(t_master));

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

fprintf("All channels merged\n");

%% ================= QC PLOT =================
figure;
subplot(2,1,1)
plot(t_master, lfp_merged_all(1,:))
title('Merged LFP - Channel 1')

% subplot(2,1,2)
% imagesc(lfp_merged_all)
% title('All channels (rows)'); xlabel('Time'); ylabel('Channel')
% colormap jet; colorbar;

%% ================= NaN VISUALIZATION =================
figure; hold on;

plot(t_master, lfp_merged_all(1,:), 'b');

nan_mask = isnan(lfp_merged_all(1,:));
d = diff([0 nan_mask 0]);
start_idx = find(d == 1);
end_idx   = find(d == -1) - 1;

yl = ylim;

for i = 1:length(start_idx)
    patch([t_master(start_idx(i)) t_master(end_idx(i)) t_master(end_idx(i)) t_master(start_idx(i))], ...
          [yl(1) yl(1) yl(2) yl(2)], ...
          [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor','none');
end

title('NaN gaps (Channel 1)');