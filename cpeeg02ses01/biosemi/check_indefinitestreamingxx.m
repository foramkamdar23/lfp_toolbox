% patient64 session 2

clear all;
clc;
close all;

sub64_s02_filepath = "C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm064\session02\Report_Json_Session_Report_1201.mat";
sub64_s03_filepath = "C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm064\session03\Report_Json_Session_Report_1332.mat"; 
sub91_s02_filepath = "C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm091\session02\Report_Json_Session_Report_1433.mat";

data1 = load(sub64_s02_filepath);
data2 = load(sub64_s03_filepath);
data3 = load(sub91_s02_filepath);

sub64_s02 = data1.data.IndefiniteStreaming;
sub64_s03 = data2.data.IndefiniteStreaming;
sub91_s02 = data3.data.IndefiniteStreaming;

% lets plot data1 first 
fs = sub64_s02.SampleRateInHz;
%ch_ids = length(sub64_s02);

for fig = 1:2
    figure;
    ch_start = (fig-1)*6 + 1;
    ch_end   = fig*6;
    for i = ch_start : ch_end
        subplot(6, 1, i - ch_start + 1);
        sig = sub64_s02(i).TimeDomainData;
        t = (0:length(sig)-1) / fs;
        plot(t, sig);
        title(sub64_s02(i).Channel);
    end
end

% lets plot data2
figure;
ch_start = 1;
ch_end   = 6;
for i = ch_start : ch_end
    subplot(6, 1, i - ch_start + 1);
    sig = sub64_s03(i).TimeDomainData;
    t = (0:length(sig)-1) / fs;
    plot(t, sig);
    title(sub64_s03(i).Channel);
end

% lets plot data3
figure;
ch_start = 1;
ch_end   = 6;
for i = ch_start : ch_end
    subplot(6, 1, i - ch_start + 1);
    sig = sub91_s02(i).TimeDomainData;
    t = (0:length(sig)-1) / fs;
    plot(t, sig);
    title(sub91_s02(i).Channel);
end

%%
% patient64 session 2
clear all;
clc;
close all;

sub64_s02_filepath = "C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm064\session02\Report_Json_Session_Report_1201.mat";
sub64_s03_filepath = "C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm064\session03\Report_Json_Session_Report_1332.mat";
sub91_s02_filepath = "C:\Users\fkamdar\Desktop\repos\data\pdnm\pdnm091\session02\Report_Json_Session_Report_1433.mat";

data1 = load(sub64_s02_filepath);
data2 = load(sub64_s03_filepath);
data3 = load(sub91_s02_filepath);

sub64_s02 = data1.data.IndefiniteStreaming;
sub64_s03 = data2.data.IndefiniteStreaming;
sub91_s02 = data3.data.IndefiniteStreaming;

fs = sub64_s02.SampleRateInHz;

% Design 80 Hz bandpass filter (75-85 Hz), shared across all plots
[b, a] = butter(4, [75 85] / (fs/2), 'bandpass');

% --- Plot data1 (sub64_s02) ---
for fig = 1:2
    figure;
    ch_start = (fig-1)*6 + 1;
    ch_end   = fig*6;
    for i = ch_start:ch_end
        subplot(6, 1, i - ch_start + 1);
        sig = sub64_s02(i).TimeDomainData;
        sig_filt = filtfilt(b, a, sig);
        t = (0:length(sig_filt)-1) / fs;
        plot(t, sig_filt);
        title(sprintf('sub64\\_s02 | %s | 75-85 Hz filtered', sub64_s02(i).Channel));
    end
end

% --- Plot data2 (sub64_s03) ---
figure;
for i = 1:6
    subplot(6, 1, i);
    sig = sub64_s03(i).TimeDomainData;
    sig_filt = filtfilt(b, a, sig);
    t = (0:length(sig_filt)-1) / fs;
    plot(t, sig_filt);
    title(sprintf('sub64\\_s03 | %s | 75-85 Hz filtered', sub64_s03(i).Channel));
end

% --- Plot data3 (sub91_s02) ---
figure;
for i = 1:6
    subplot(6, 1, i);
    sig = sub91_s02(i).TimeDomainData;
    sig_filt = filtfilt(b, a, sig);
    t = (0:length(sig_filt)-1) / fs;
    plot(t, sig_filt);
    title(sprintf('sub91\\_s02 | %s | 75-85 Hz filtered', sub91_s02(i).Channel));
end

%%
