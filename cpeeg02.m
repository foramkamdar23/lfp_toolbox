%CPEEG02

load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124017_firsthalf_1437.mat");

new_datastruct = data.IndefiniteStreaming;

fs = new_datastruct(1).SampleRateInHz;
rec_start_datetime = new_datastruct(1).FirstPacketDateTime_PDT;

processed_data = struct();
for i = 1:length(new_datastruct)

    % Reconstruct signal with NaNs
    [signal_with_nan, info] = detect_packet_loss_percept( ...
        new_datastruct(i).TimeDomainData, ...
        new_datastruct(i).GlobalSequences, ...
        new_datastruct(i).GlobalPacketSizes, ...
        fs);

    % Time vector in datetime
    t = rec_start_datetime + seconds((0:length(signal_with_nan)-1)/fs);

    % Populate structure
    processed_data(i).channel       = new_datastruct(i).Channel;
    processed_data(i).signal        = signal_with_nan;
    processed_data(i).time          = t;
    processed_data(i).fs            = fs;
    processed_data(i).gap_info      = info.gap_info;
    processed_data(i).packet_loss_rate    = info.packet_loss_rate;
    processed_data(i).bad_packet_indices  = info.bad_packet_indices;

    % Optional: plotting (can be removed if not needed)
    figure(i)
    plot(t, signal_with_nan)
    hold on
    nan_mask = isnan(signal_with_nan);
    d = diff([0; nan_mask; 0]);
    start_idx = find(d == 1);
    end_idx   = find(d == -1) - 1;
    yl = ylim;
    for k = 1:length(start_idx)
        patch([t(start_idx(k)) t(end_idx(k)) t(end_idx(k)) t(start_idx(k))], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 0 0], 'EdgeColor','none');
    end
    title(sprintf('%s | Packet loss: %.3f%%', strrep(new_datastruct(i).Channel,"_"," "), ...
        info.packet_loss_rate*100))
    xlabel('Time of day')
    ylabel('Amplitude')
end



% fieldtrip 
addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;

%load eeg data
dataDir = 'C:\Users\saosorio\Projects\WorkingMemory_CPEEG';
eegFile = 'cpeeg01_b02_offstimoffmed_wm.bdf';
eegFullPath = fullfile(dataDir, eegFile);

cfg = [];
cfg.dataset = eegFullPath;
cfg.continuous = 'yes';

data = ft_preprocessing(cfg);
