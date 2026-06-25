%Try to load .mat file and LFP data over time

clear all
clc
close all

filepath = "C:\Users\fkamdar\Desktop\percept\code\testing\Report_Json_Session_Report_20250903T170810.mat";
load(filepath)


% for ch_id = 1:length(data.BrainSenseTimeDomain)
%     bs = data.BrainSenseTimeDomain(ch_id);
%     fprintf('Channel: %s | Fs: %d Hz\n', bs.Channel, bs.SampleRateInHz);
% end

%% Lets do one channel at a time ( out of the 7)
ch_id = 1; 
bs = data.BrainSenseTimeDomain(ch_id);   % choose channel


% Check Packet Drops  % Does not work
seqDiff = diff(bs.GlobalSequences);
missingPackets = seqDiff(seqDiff > 1) - 1;
totalLostPackets = sum(missingPackets) % It should be clean ideally

% Use this to check for missing packets
packetSizes = str2double(split(bs.GlobalPacketSizes, ','));
packetSizes = packetSizes(~isnan(packetSizes)); % remove empty entries
packetSize = unique(packetSizes);    
% expectedSamples and actualSamples should match sum(packetSizes) == length(TimeDomainData)

expectedSamples = sum(packetSizes);
actualSamples   = length(bs.TimeDomainData);

% Sampling frequency and time vectors 
fs = bs.SampleRateInHz;
n = length(bs.TimeDomainData);
duration_secs = n/fs;
lfp = double(bs.TimeDomainData);   % make sure numeric
t = (0:length(lfp)-1) / fs;


% Plot Raw LFP  %typical lfp std (confirm) 
figure(1);
plot(t, lfp)
xlabel('Time (seconds)')
ylabel('Amplitude (\muV)')
title(['BrainSense LFP - ' bs.Channel])


% Delete outlier peaks 
mu = mean(lfp); sd = std(lfp);
scale_k = 10;

isArt   = abs(lfp-mu) > scale_k *sd;  %is Artifacts
isnotArt =  ~isArt;

lfp_clean = lfp;

lfp_clean(isArt) = interp1(t(isnotArt), lfp(isnotArt), t(isArt), 'linear',  'extrap');


% Plot LFP with outlier peak removed
figure(2);
plot(t, lfp_clean)
xlabel('Time (seconds)')
ylabel('Amplitude (\muV)')
title(['BrainSense LFP - ' bs.Channel])

%Plot power spectrum
figure(3);
nyquist_freq = fs/2;
windowLength = 2*fs;   % 2-second window (500 samples)
noverlap     = [];     % MATLAB default (50%)
nfft         = [];     % Let MATLAB choose

[pxx, freq] = pwelch(lfp_clean, windowLength, noverlap, nfft, fs);

figure(3);
plot(freq,10*log10(pxx))
xlim([0 nyquist_freq])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title(['BrainSense LFP - ' bs.Channel])


