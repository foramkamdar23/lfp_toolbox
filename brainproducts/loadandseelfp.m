%PDNM2_s1
clear all;
close all;
clc;
%% Cell 1
data2 = load("C:\Users\fkamdar\Desktop\repos\data\Report_Json_Session_Report_20260618T132944.mat");
brainsense = data2.data.BrainSenseTimeDomain;
indefinite = data2.data.IndefiniteStreaming;
%% Cell 2 load and filter lfp across channels
n_bs    = numel(brainsense);   % expected 4
n_indef = numel(indefinite);   % expected 6

% preallocate
sig1      = cell(n_bs,1);
t_lfp1    = cell(n_bs,1);
lfp1_filt = cell(n_bs,1);

sig2      = cell(n_indef,1);
t_lfp2    = cell(n_indef,1);
lfp2_filt = cell(n_indef,1);

% --- BrainSense channels ---
for ch = 1:n_bs
    bs_ch    = brainsense(ch);
    sig1{ch} = bs_ch.TimeDomainData;
    fs_bs    = bs_ch.SampleRateInHz;

    t_lfp1{ch} = (0:length(sig1{ch})-1)/fs_bs;

    [b,a] = butter(4, [75 85]/(fs_bs/2), 'bandpass');
    lfp1_filt{ch} = filtfilt(b,a, sig1{ch});

    figure;
    subplot(2,1,1); plot(t_lfp1{ch}, sig1{ch});      title(sprintf('BrainSense raw - ch %d', ch));
    subplot(2,1,2); plot(lfp1_filt{ch});             title(sprintf('BrainSense filt - ch %d', ch));
end

% --- Indefinite channels ---
for ch = 1:n_indef
    indef_ch = indefinite(ch);
    sig2{ch} = indef_ch.TimeDomainData;
    fs_indef = indef_ch.SampleRateInHz;

    t_lfp2{ch} = (0:length(sig2{ch})-1)/fs_indef;

    [b,a] = butter(4, [78 82]/(fs_indef/2), 'bandpass');
    lfp2_filt{ch} = filtfilt(b,a, sig2{ch});

    figure;
    subplot(2,1,1); plot(t_lfp2{ch}, sig2{ch});      title(sprintf('Indefinite raw - ch %d', ch));
    subplot(2,1,2); plot(lfp2_filt{ch});             title(sprintf('Indefinite filt - ch %d', ch));
end