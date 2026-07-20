%PDNM2_s1
clear all;
close all;
clc;

%% Cell 1 LOAD LFP
data1 = load("C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm091\session04\Report_Json_Session_Report_20260709T173126_WM_BLOCK02.mat");
first_half = data1.data.IndefiniteStreaming;

ch_id = 1;
ch1 = first_half(1);
ch2 = first_half(2);
ch3 = first_half(3);
ch4 = first_half(4);
ch5 = first_half(5);
ch6 = first_half(6);

sig1 = ch1.TimeDomainData;
sig2 = ch2.TimeDomainData;
sig3 = ch3.TimeDomainData;
sig4 = ch4.TimeDomainData;
sig5 = ch5.TimeDomainData;
sig6 = ch6.TimeDomainData;


fs_lfp = ch1.SampleRateInHz;

t_lfp = (0:length(sig1)-1)/fs_lfp;

sigs = {sig1, sig2, sig3, sig4, sig5, sig6};

figure;
for i = 1:6
    sig = sigs{i};
    subplot(6,1,i);
    plot(t_lfp, sig);
    title(['LFP Channel ' num2str(i)]);
end



%% Cell 2 LOAD CLAVICLE EEG

% fieldtrip 
addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

vhdr = 'C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm091\session04\PDNM091_0709026_s4_stim_OFF_medON_WM_b02_1311.vhdr';
vmrk = 'C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm091\session04\PDNM091_0709026_s4_stim_OFF_medON_WM_b02_1311.vmrk';
eeg_data = 'C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm091\session04\PDNM091_0709026_s4_stim_OFF_medON_WM_b02_1311.eeg';

cfg.dataset = vhdr;  % .vhdr is the header file to start with
cfg.channel = 'Right Clavicle'; % Replace with your channel name, e.g., 'EOG' or 'EEG001'
right_clavicle_data = ft_preprocessing(cfg);


right_clavicle_eeg = right_clavicle_data.trial{1};
t_eeg = right_clavicle_data.time{1};
fs_eeg = right_clavicle_data.fsample;

figure;
plot(t_eeg, right_clavicle_eeg);
title('Right Clavicle');

% down sampling
right_clavicle_ds = resample(right_clavicle_eeg, fs_lfp, fs_eeg);
t_right_clavicle_ds = (0:length(right_clavicle_ds)-1)/fs_lfp;

%qc
figure;
subplot(2,1,1); plot(t_eeg, right_clavicle_eeg); title('Original EEG');
subplot(2,1,2); plot(t_right_clavicle_ds, right_clavicle_ds); title('Downsampled EEG (250 Hz)');

close all



% filter for enevelope around 80hz tens burst
[b,a] = butter(4, [78 82]/(fs_lfp/2), 'bandpass');

right_clavicle_ds_filt = filtfilt(b,a, right_clavicle_ds);
ch1_filt = filtfilt(b,a, sig1);
ch2_filt = filtfilt(b,a, sig2);
ch3_filt = filtfilt(b,a, sig3);
ch4_filt = filtfilt(b,a, sig4);
ch5_filt = filtfilt(b,a, sig5);
ch6_filt = filtfilt(b,a, sig6);

chs_filt = {ch1_filt, ch2_filt, ch3_filt, ch4_filt, ch5_filt, ch6_filt};

figure;

for i = 1:6
    ch_filt = chs_filt{i};
    subplot(7,1,i);
    plot(ch_filt);
    title(['LFP Channel ' num2str(i)]);
end
subplot(7,1,7); plot(right_clavicle_ds_filt); title('EEG'); 




subplot(3,1,2); plot(ch1_filt); title('ch1');


% for loop ?

%env created and normalized
right_clavicle_env = abs(hilbert(right_clavicle_ds_filt));
ch1_env = abs(hilbert(ch1_filt));

right_clavicle_env = zscore(right_clavicle_env); 
ch1_env = zscore(ch1_env); 

%cross correlation
[xc1, lags1] = xcorr(right_clavicle_env, ch1_env);
[~, idx1] = max(xc1); % find peak
lag_samples1 = lags1(idx1);


% convert to seconds
lag_sec1 = lag_samples1 / fs_lfp;


%align here
t_lfp_aligned = t_lfp + lag_sec1;

figure;
plot(t_right_clavicle_ds, right_clavicle_env); hold on;
plot(t_lfp_aligned, ch1_env);
legend('EEG ENV','LFP ENV (aligned)');

fprintf("done with alinging")


figure('Name','Envelope alignment');
plot(t_right_clavicle_ds, right_clavicle_env); hold on;
plot(t_lfp_aligned, ch1_env);
legend('EEG env','LFP1 env (aligned)');
title('TENS envelope alignment'); xlabel('Time (s)');
