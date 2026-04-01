%CPEEG02
clear all;
close all;
clc;

load("C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected\Report_Json_Session_Report_20260324T124034_secondhalf_1453.mat");

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
dataDir = 'C:\Users\fkamdar\Desktop\repos\data\eeg-selected';
eegFile = 'cpeeg02_b01_ert.bdf';
eegFullPath = fullfile(dataDir, eegFile);

cfg = [];
cfg.dataset = eegFullPath;
cfg.continuous = 'yes';

data = ft_preprocessing(cfg);


%% -------------------------------------------------------------------------
% BLOCK 3 — Select EEG + EXG channels
% -------------------------------------------------------------------------
cfg = [];
cfg.channel = ft_channelselection({'EEG*','EXG*'}, data.label);

data_clean = ft_selectdata(cfg, data);


%% -------------------------------------------------------------------------
% BLOCK 4 — Downsample EEG to match Percept sampling rate
% -------------------------------------------------------------------------
cfg = [];
cfg.resamplefs = 250;
cfg.detrend = 'no';

data_ds = ft_resampledata(cfg, data_clean);

Fs_eeg = data_ds.fsample;

%% -------------------------------------------------------------------------
% BLOCK 5 — Extract EXG7 channel (contains TENS artifact)
% -------------------------------------------------------------------------
channelName = 'EXG7';

cfg = [];
cfg.channel = channelName;

data_exg = ft_selectdata(cfg, data_ds);

%% -------------------------------------------------------------------------
% BLOCK 6 — Filter EXG7 and compute envelope
% -------------------------------------------------------------------------
cfg = [];
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [75 85];
cfg.bpfilttype = 'but';
cfg.bpfiltord  = 4;

% Hilbert transform returns the analytic signal envelope
cfg.hilbert = 'abs';

data_exg_filt = ft_preprocessing(cfg, data_exg);

timeVec_eeg = data_exg_filt.time{1};
signal_eeg_env = data_exg_filt.trial{1}(1,:);


%% -------------------------------------------------------------------------
% BLOCK 7 — Normalize EEG signal
% -------------------------------------------------------------------------
signal_eeg_env = zscore(signal_eeg_env);

%% -------------------------------------------------------------------------
% BLOCK 8 — Plot filtered EEG artifact
% -------------------------------------------------------------------------
figure('Color','w');

plot(timeVec_eeg, signal_eeg_env,'LineWidth',1.2)

xlabel('Time (s)')
ylabel('Z-scored amplitude')
title('EXG7 TENS Artifact (Filtered + Envelope)')
grid on

%% -------------------------------------------------------------------------
% BLOCK 9 — Percept intracranial data
% -------------------------------------------------------------------------

signal_orig = processed_data(1).signal;
Fs_percept = processed_data(1).fs;

%% -------------------------------------------------------------------------
% BLOCK 10 — Create synchronization copy and interpolate NaNs
% -------------------------------------------------------------------------
signal_sync = signal_orig;
nanIdx = isnan(signal_sync);
signal_sync(nanIdx) = interp1( ...
    find(~nanIdx), ...
    signal_sync(~nanIdx), ...
    find(nanIdx), ...
    'linear','extrap');

%% -------------------------------------------------------------------------
% BLOCK 11 — Bandpass filter Percept signal
% -------------------------------------------------------------------------
lowCut = 75;
highCut = 85;
order = 4;

[b,a] = butter(order,[lowCut highCut]/(Fs_percept/2),'bandpass');

signal_percept_filt = filtfilt(b,a,signal_sync);

%% -------------------------------------------------------------------------
% BLOCK 12 — Compute envelope and normalize Percept signal
% -------------------------------------------------------------------------
signal_percept_env = abs(signal_percept_filt);
signal_percept_env = zscore(signal_percept_env);


%% -------------------------------------------------------------------------
% BLOCK 13 — Visual comparison of signals
% -------------------------------------------------------------------------
figure('Color','w')

plot(signal_percept_env,'b')
hold on
plot(signal_eeg_env,'r', 'LineWidth',1)

legend('Percept','EEG EXG7')
title('Z-scored TENS artifacts')
xlabel('Samples')
ylabel('Amplitude')
grid on

%% -------------------------------------------------------------------------
% BLOCK 14 — Estimate synchronization lag
% -------------------------------------------------------------------------
[xc,lags] = xcorr(signal_percept_env,signal_eeg_env);
xc = xc ./ max(abs(xc));
[~,idx] = max(xc);
lag_samples = lags(idx);
lag_seconds = lag_samples / Fs_percept;
fprintf('Estimated lag: %d samples (%.3f seconds)\n',lag_samples,lag_seconds)

%% -------------------------------------------------------------------------
% BLOCK 15 — Plot cross-correlation
% -------------------------------------------------------------------------

figure('Color','w')
plot(lags/Fs_percept,xc,'LineWidth',1.2)
xlabel('Lag (seconds)')
ylabel('Correlation')
title('Cross-correlation between Percept and EEG artifacts')
grid on

%% -------------------------------------------------------------------------
% BLOCK 16 — Align Percept to EEG using estimated lag
% -------------------------------------------------------------------------
% Align both ORIGINAL signal (with NaNs) AND envelope signal
lag_abs = abs(lag_samples);

if lag_samples > 0
    % Percept is ahead → shift LEFT
    signal_percept_aligned     = [signal_orig; nan(lag_abs,1)];
    signal_percept_env_aligned = [signal_percept_env; nan(lag_abs,1)]; 
    
elseif lag_samples < 0
    % Percept is behind → shift RIGHT
    signal_percept_aligned     = [nan(lag_abs,1); signal_orig]; 
    signal_percept_env_aligned = [nan(lag_abs,1); signal_percept_env];
    
else
    % No lag
    signal_percept_aligned     = signal_orig;
    signal_percept_env_aligned = signal_percept_env;
end
%% -------------------------------------------------------------------------
% BLOCK 17 — Create common time base
% -------------------------------------------------------------------------
t_eeg = (0:length(signal_eeg_env)-1) / Fs_eeg;
t_percept = (0:length(signal_percept_aligned)-1) / Fs_percept;
%% -------------------------------------------------------------------------
% BLOCK 18a — Plot EEG vs Percept envelope (signal_percept_env_aligned)
% -------------------------------------------------------------------------
figure('Color','w'); hold on;

% Plot signals first
plot(t_percept, signal_percept_env_aligned, 'b')
plot(t_eeg, zscore(signal_eeg_env), 'r', 'LineWidth',1)

% Identify NaN segments (based on original aligned signal)
nan_mask = isnan(signal_percept_aligned);
nan_mask = nan_mask(:)';  % ensure row
d = diff([0 nan_mask 0]);
nan_start = find(d == 1);
nan_end   = find(d == -1) - 1;

% Get current y-limits
yl = ylim;

% Plot NaN patches
for i = 1:length(nan_start)
    x_patch = [t_percept(nan_start(i)) t_percept(nan_end(i)) t_percept(nan_end(i)) t_percept(nan_start(i))];
    y_patch = [yl(1) yl(1) yl(2) yl(2)];
    patch(x_patch, y_patch, 'k', 'FaceAlpha', 0.5, 'EdgeColor','none');
end

% Bring signals to top
plot(t_percept, signal_percept_env_aligned, 'b')
plot(t_eeg, zscore(signal_eeg_env), 'r', 'LineWidth',1)

legend('Percept envelope','EEG EXG7','Missing Packets')
xlabel('Time (s)')
ylabel('Amplitude (z-score)')
title('Aligned EEG and Percept envelope (NaNs highlighted)')
grid on


%% -------------------------------------------------------------------------
% BLOCK 18b — Plot EEG vs z-scored original Percept (signal_percept_z)
% -------------------------------------------------------------------------
% Compute z-scored aligned Percept signal
signal_percept_z = (signal_percept_aligned - nanmean(signal_percept_aligned(~isnan(signal_percept_aligned)))) ./ ...
                   nanstd(signal_percept_aligned(~isnan(signal_percept_aligned)));

figure('Color','w'); hold on;

% Plot signals first
plot(t_percept, signal_percept_z, 'b')
plot(t_eeg, zscore(signal_eeg_env), 'r', 'LineWidth',1)

% Get y-limits
yl = ylim;

% Plot NaN patches
for i = 1:length(nan_start)
    x_patch = [t_percept(nan_start(i)) t_percept(nan_end(i)) t_percept(nan_end(i)) t_percept(nan_start(i))];
    y_patch = [yl(1) yl(1) yl(2) yl(2)];
    patch(x_patch, y_patch, 'k', 'FaceAlpha', 0.3, 'EdgeColor','none');
end

% Bring signals to top
plot(t_percept, signal_percept_z, 'b')
plot(t_eeg, zscore(signal_eeg_env), 'r', 'LineWidth',1)

legend('NaN segments','Percept (z-scored)','EEG EXG7')
xlabel('Time (s)')
ylabel('Amplitude (z-score)')
title('Aligned EEG and Percept z-scored signal (NaNs highlighted)')
grid on



%% -------------------------------------------------------------------------
% BLOCK 19 — Build continuous FieldTrip structure with aligned Percept LFP
% -------------------------------------------------------------------------
% Use t_eeg as the time base: it matches the data matrix length exactly,
% avoiding sampleinfo mismatches in ft_checkdata. The Percept clock is
% preserved through the lag-shifted trl indices in Block 20.

nChannels = length(processed_data);
nSamples  = length(t_eeg);             % EEG is the length reference

% --- 19a: Align ALL Percept channels with the lag from Block 14 -----------
percept_aligned_all = nan(nChannels, nSamples);

for ch = 1:nChannels

    sig_ch = processed_data(ch).signal;     % original signal, NaNs intact

    if lag_samples > 0
        sig_shifted = [sig_ch; nan(lag_abs, 1)];    % Percept ahead → pad end
    elseif lag_samples < 0
        sig_shifted = [nan(lag_abs, 1); sig_ch];    % Percept behind → pad start
    else
        sig_shifted = sig_ch;
    end

    % Trim or pad to exactly nSamples
    len = length(sig_shifted);
    if len >= nSamples
        percept_aligned_all(ch, :) = sig_shifted(1:nSamples);
    else
        percept_aligned_all(ch, 1:len) = sig_shifted;
    end

end

% --- 19b: Channel metadata ------------------------------------------------
chan_labels = {processed_data.channel}';    % Nx1 cell, e.g. 'ZERO_THREE_LEFT'
chan_type   = repmat({'lfp'}, nChannels, 1);
chan_units  = repmat({'uV'},  nChannels, 1);

% --- 19c: Assemble the raw continuous structure ---------------------------
data_percept_ft              = [];
data_percept_ft.label        = chan_labels;
data_percept_ft.fsample      = Fs_percept;          % 250 Hz
data_percept_ft.trial        = {percept_aligned_all};
data_percept_ft.time         = {t_eeg};             % length matches data matrix
data_percept_ft.sampleinfo   = [1, nSamples];

% Header — required by some ft_ functions
data_percept_ft.hdr.label        = chan_labels;
data_percept_ft.hdr.chantype     = chan_type;
data_percept_ft.hdr.chanunit     = chan_units;
data_percept_ft.hdr.Fs           = Fs_percept;
data_percept_ft.hdr.nChans       = nChannels;
data_percept_ft.hdr.nSamples     = nSamples;
data_percept_ft.hdr.nTrials      = 1;
data_percept_ft.hdr.nSamplesPre  = 0;

% Percept quality metadata — survives ft_checkdata as a cfg subfield
data_percept_ft.cfg.percept_gap_info    = {processed_data.gap_info};
data_percept_ft.cfg.percept_packet_loss = [processed_data.packet_loss_rate];
data_percept_ft.cfg.percept_bad_packets = {processed_data.bad_packet_indices};
data_percept_ft.cfg.sync_lag_samples    = lag_samples;
data_percept_ft.cfg.sync_lag_seconds    = lag_seconds;

% Validate
data_percept_ft = ft_checkdata(data_percept_ft, 'datatype', 'raw', ...
                               'feedback', 'yes', 'hassampleinfo', 'yes');

fprintf('Percept FT structure ready: %d channels x %d samples at %d Hz\n', ...
        nChannels, nSamples, Fs_percept);


%% -------------------------------------------------------------------------
% BLOCK 20 — Define trials from BDF events and epoch Percept LFP
%            + REMOVE TENS-CONTAMINATED TRIALS
% -------------------------------------------------------------------------

% --- 20a: Read and correct BioSemi STATUS events --------------------------
events_raw    = ft_read_event(eegFullPath);
status_events = events_raw(strcmp({events_raw.type}, 'STATUS'));

% Handle both numeric and string value fields
if ischar(status_events(1).value) || isstring(status_events(1).value)
    raw_values = str2double({status_events.value});
else
    raw_values = double([status_events.value]);
end

corrected_values = bitand(uint16(raw_values) - 512, 255);

% Print summary
fprintf('\nCorrected BioSemi event codes:\n');
unique_vals = unique(corrected_values);
for i = 1:length(unique_vals)
    v     = unique_vals(i);
    n     = sum(corrected_values == v);
    first = status_events(find(corrected_values == v, 1)).sample;
    fprintf('Code %d | Count %d | First sample %d\n', v, n, first);
end

% --- Epoch settings -------------------------------------------------------
PRE_STIM_S  = 0;
POST_STIM_S = 1.0;
EPOCH_CODE  = 30;   % image onset

pre_samples  = round(PRE_STIM_S  * Fs_percept);
post_samples = round(POST_STIM_S * Fs_percept);

Fs_bdf = data.fsample;

% --- Extract event samples ------------------------------------------------
epoch_events  = status_events(corrected_values == EPOCH_CODE);
epoch_samples = [epoch_events.sample];

% Convert BDF → Percept sample space
epoch_samples_ds = round(epoch_samples * (Fs_percept / Fs_bdf));

% --- Build initial trl ----------------------------------------------------
trl_percept = [];

for k = 1:length(epoch_samples_ds)

    beg = epoch_samples_ds(k) - pre_samples  + lag_samples;
    en  = epoch_samples_ds(k) + post_samples + lag_samples;

    trl_percept(end+1,:) = [beg, en, -pre_samples, 30]; %#ok<AGROW>
end

% --- Remove out-of-bounds trials -----------------------------------------
nSamples = size(data_percept_ft.trial{1}, 2);

valid = trl_percept(:,1) >= 1 & trl_percept(:,2) <= nSamples;

fprintf('%d / %d trials removed (out of bounds)\n', ...
    sum(~valid), size(trl_percept,1));

trl_percept = trl_percept(valid,:);

%% -------------------------------------------------------------------------
% 🔥 NEW: REMOVE TENS-CONTAMINATED TRIALS USING ENVELOPE
% -------------------------------------------------------------------------

threshold = 3;   % tune (2.5–4 works well)

good_trials = true(size(trl_percept,1),1);

for k = 1:size(trl_percept,1)

    beg = trl_percept(k,1);
    en  = trl_percept(k,2);

    % safety
    if beg < 1 || en > length(signal_percept_env_aligned)
        good_trials(k) = false;
        continue;
    end

    seg = signal_percept_env_aligned(beg:en);

    % robust detection (better than max)
    if mean(seg > 2.5) > 0.1   % >10% contaminated
        good_trials(k) = false;
    end
end

fprintf('Removed %d / %d trials due to TENS\n', ...
    sum(~good_trials), length(good_trials));

trl_percept = trl_percept(good_trials,:);


%% -------------------------------------------------------------------------
% DEBUG: WHICH TRIALS WERE REMOVED
% -------------------------------------------------------------------------

removed_idx = find(~good_trials);
kept_idx    = find(good_trials);

fprintf('\n--- TENS REJECTION SUMMARY ---\n');
fprintf('Total trials: %d\n', length(good_trials));
fprintf('Kept trials:  %d\n', length(kept_idx));
fprintf('Removed:      %d\n\n', length(removed_idx));

fprintf('Removed trial indices:\n');
disp(removed_idx');

%% -------------------------------------------------------------------------
% Epoch using CLEAN trials
% -------------------------------------------------------------------------

cfg_epoch     = [];
cfg_epoch.trl = trl_percept;

data_percept_epoched = ft_redefinetrial(cfg_epoch, data_percept_ft);

% store event code
data_percept_epoched.trialinfo = trl_percept(:,4);

fprintf('\nFINAL: %d trials after TENS removal\n', ...
    length(data_percept_epoched.trial));

%% -------------------------------------------------------------------------
% BLOCK 21 — Visual sanity check: all epochs per channel
% -------------------------------------------------------------------------

time_ax = data_percept_epoched.time{1};
nTrials = length(data_percept_epoched.trial);
nChannels = length(data_percept_epoched.label);

% Use unique trial types (e.g., 30 for image)
trial_types = unique(data_percept_epoched.trialinfo);
colors = lines(length(trial_types));

%% --- 21a: Identify NaN-contaminated trials -------------------------------
nan_trial_mask = false(nTrials, nChannels);

for tr = 1:nTrials
    for ch = 1:nChannels
        nan_trial_mask(tr, ch) = any(isnan(data_percept_epoched.trial{tr}(ch, :)));
    end
end

% Trial rejected if ANY channel has NaN
nan_any_ch = any(nan_trial_mask, 2);

n_nan   = sum(nan_any_ch);
pct_nan = 100 * n_nan / nTrials;

%% --- Print summary -------------------------------------------------------
fprintf('\n--- NaN trial report ---\n');
fprintf('Total trials:    %d\n', nTrials);
fprintf('Trials with NaN: %d (%.1f%%)\n', n_nan, pct_nan);
fprintf('Clean trials:    %d (%.1f%%)\n', nTrials - n_nan, 100 - pct_nan);

% Per-channel breakdown
fprintf('\nPer-channel breakdown:\n');
fprintf('%-20s  %-10s  %s\n', 'Channel', 'NaN trials', '%');
fprintf('%s\n', repmat('-', 1, 40));

for ch = 1:nChannels
    n_ch   = sum(nan_trial_mask(:, ch));
    pct_ch = 100 * n_ch / nTrials;

    fprintf('%-20s  %-10d  %.1f%%\n', ...
        data_percept_epoched.label{ch}, n_ch, pct_ch);
end

% Per-condition breakdown
fprintf('\nPer-condition breakdown:\n');
fprintf('%-10s  %-10s  %s\n', 'Code', 'NaN trials', '%');
fprintf('%s\n', repmat('-', 1, 35));

for i = 1:length(trial_types)
    code = trial_types(i);

    mask     = data_percept_epoched.trialinfo == code;
    n_code   = sum(nan_any_ch & mask);
    tot_code = sum(mask);

    fprintf('%-10d  %-10d  %.1f%%\n', ...
        code, n_code, 100 * n_code / tot_code);
end

% %% --- 21b: Plot -----------------------------------------------------------
% figure('Color','w', 'Position', [100 100 1400 900]);
% 
% for ch = 1:nChannels
% 
%     subplot(ceil(nChannels/2), 2, ch); hold on;
% 
%     for tr = 1:nTrials
% 
%         code   = data_percept_epoched.trialinfo(tr);
%         signal = data_percept_epoched.trial{tr}(ch, :);
% 
%         color_idx = find(trial_types == code);
% 
%         if nan_trial_mask(tr, ch)
%             % NaN trials → dashed red
%             plot(time_ax, signal, ...
%                 'Color', [0.85 0.15 0.15 0.3], ...
%                 'LineWidth', 0.6, ...
%                 'LineStyle', '--');
%         else
%             % Clean trials → colored by condition
%             plot(time_ax, signal, ...
%                 'Color', [colors(color_idx,:), 0.4], ...
%                 'LineWidth', 0.8);
%         end
%     end
% 
%     xline(0, 'k--', 'LineWidth', 1.2);
% 
%     xlabel('Time (s)');
%     ylabel('Amplitude (\muV)');
% 
%     title(sprintf('%s | NaN: %d/%d', ...
%         data_percept_epoched.label{ch}, ...
%         sum(nan_trial_mask(:, ch)), nTrials), ...
%         'Interpreter', 'none');
% 
%     xlim([time_ax(1) time_ax(end)]);
%     grid on; box off;
% 
% end
% 
% %% --- Legend --------------------------------------------------------------
% ax_leg = axes('Position', [0.92 0.3 0.01 0.3], 'Visible', 'off');
% 
% for i = 1:length(trial_types)
%     line(ax_leg, NaN, NaN, ...
%         'Color', colors(i,:), ...
%         'LineWidth', 2, ...
%         'DisplayName', sprintf('Code %d', trial_types(i)));
% end
% 
% line(ax_leg, NaN, NaN, ...
%     'Color', [0.85 0.15 0.15], ...
%     'LineWidth', 1.5, ...
%     'LineStyle', '--', ...
%     'DisplayName', 'NaN trial');
% 
% legend(ax_leg, 'show', 'Location', 'west');
% 
% sgtitle(sprintf('All epochs — %d trials | %d NaN (%.1f%%)', ...
%     nTrials, n_nan, pct_nan), ...
%     'FontSize', 13, 'FontWeight', 'bold');
% 
% %% -------------------------------------------------------------------------
% % BLOCK 22 — Grand average across all trials and digit positions
% % -------------------------------------------------------------------------
% 
% % --- 22a: Exclude NaN trials ----------------------------------------------
% clean_trials = find(~nan_any_ch);
% 
% cfg_clean        = [];
% cfg_clean.trials = clean_trials;
% 
% data_clean = ft_selectdata(cfg_clean, data_percept_epoched);
% 
% % --- 22b: Compute grand average (keep trials to compute SEM) -------------
% cfg_ga            = [];
% cfg_ga.trials     = 'all';
% cfg_ga.keeptrials = 'yes';
% 
% grandavg = ft_timelockanalysis(cfg_ga, data_clean);
% 
% % --- 22c: Baseline correction (−0.5 to 0 s) ------------------------------
% cfg_bl          = [];
% cfg_bl.baseline = [-0.5 0];
% 
% grandavg_bl = ft_timelockbaseline(cfg_bl, grandavg);
% 
% % Compute mean and SEM manually from the kept trials (nTrials x nChan x nTime)
% avg_signal = squeeze(mean(grandavg_bl.trial, 1));          % nChan x nTime
% sem_signal = squeeze(std(grandavg_bl.trial, 0, 1)) ...
%              ./ sqrt(length(clean_trials));                 % nChan x nTime
% 
% % --- 22d: Plot ------------------------------------------------------------
% time_ax = grandavg_bl.time;
% 
% figure('Color','w', 'Position', [100 100 1400 900]);
% 
% for ch = 1:nChannels
% 
%     subplot(2, 3, ch); hold on;
% 
%     sem = sem_signal(ch, :);
%     avg = avg_signal(ch, :);
% 
%     fill([time_ax fliplr(time_ax)], ...
%          [avg + sem, fliplr(avg - sem)], ...
%          [0.4 0.6 0.9], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
% 
%     plot(time_ax, avg, 'Color', [0.15 0.35 0.75], 'LineWidth', 2);
% 
%     xlims = [time_ax(1) time_ax(end)];
%     yl    = ylim;
%     fill([-0.5 0 0 -0.5], [yl(1) yl(1) yl(2) yl(2)], ...
%          [0.85 0.85 0.85], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     xline(0, 'k--', 'LineWidth', 1.2);
%     xlabel('Time (s)');
%     ylabel('\Delta Amplitude (\muV)');
%     title(sprintf('%s  |  n = %d', grandavg_bl.label{ch}, length(clean_trials)), ...
%           'Interpreter','none');
%     xlim(xlims);
%     grid on; box off;
% 
% end
% 
% sgtitle(sprintf('Grand average (baselined −0.5–0 s)  |  %d clean trials (%.1f%%)', ...
%         length(clean_trials), 100*length(clean_trials)/nTrials), ...
%         'FontSize', 13, 'FontWeight', 'bold');
% 
% 
% 
% %% -------------------------------------------------------------------------
% % BLOCK 22 — Grand average (LFP)
% % -------------------------------------------------------------------------
% 
% %% --- 22a: Remove NaN trials ---------------------------------------------
% clean_trials = find(~nan_any_ch);
% 
% cfg = [];
% cfg.trials = clean_trials;
% 
% data_clean = ft_selectdata(cfg, data_percept_epoched);
% 
% fprintf('Using %d clean trials (%.1f%%)\n', ...
%     length(clean_trials), 100*length(clean_trials)/length(nan_any_ch));
% 
% %% --- 22b: Timelock (keep trials for SEM) --------------------------------
% cfg = [];
% cfg.keeptrials = 'yes';
% 
% grandavg = ft_timelockanalysis(cfg, data_clean);
% 
% %% --- 22c: Baseline correction -------------------------------------------
% cfg = [];
% cfg.baseline = [-1 -0.5];   % adjust if needed
% 
% grandavg_bl = ft_timelockbaseline(cfg, grandavg);
% 
% %% --- 22d: Compute mean + SEM --------------------------------------------
% avg_signal = squeeze(mean(grandavg_bl.trial, 1));   % [chan x time]
% sem_signal = squeeze(std(grandavg_bl.trial, 0, 1)) ...
%              ./ sqrt(size(grandavg_bl.trial,1));
% 
% time_ax   = grandavg_bl.time;
% nChannels = length(grandavg_bl.label);
% 
% %% --- 22e: Plot -----------------------------------------------------------
% figure('Color','w', 'Position', [100 100 1400 900]);
% 
% for ch = 1:nChannels
% 
%     subplot(ceil(nChannels/2), 2, ch); hold on;
% 
%     avg = avg_signal(ch,:);
%     sem = sem_signal(ch,:);
% 
%     % SEM shading
%     fill([time_ax fliplr(time_ax)], ...
%          [avg+sem fliplr(avg-sem)], ...
%          [0.4 0.6 0.9], ...
%          'FaceAlpha', 0.25, 'EdgeColor', 'none');
% 
%     % Mean signal
%     plot(time_ax, avg, 'Color', [0.15 0.35 0.75], 'LineWidth', 2);
% 
%     % Baseline window shading
%     yl = ylim;
%     fill([-0.5 0 0 -0.5], [yl(1) yl(1) yl(2) yl(2)], ...
%          [0.85 0.85 0.85], ...
%          'FaceAlpha', 0.3, 'EdgeColor', 'none');
% 
%     % Stim onset
%     xline(0, 'k--', 'LineWidth', 1.2);
% 
%     xlabel('Time (s)');
%     ylabel('\Delta Amplitude (\muV)');
% 
%     title(sprintf('%s | n = %d', ...
%         grandavg_bl.label{ch}, length(clean_trials)), ...
%         'Interpreter','none');
% 
%     xlim([time_ax(1) time_ax(end)]);
%     grid on; box off;
% 
% end
% 
% sgtitle(sprintf('Grand average (baseline -0.5–0 s) | %d trials', ...
%     length(clean_trials)), ...
%     'FontSize', 13, 'FontWeight', 'bold');
% 
% 



% % Keep only first 48 trials
% data_percept_epoched.trial = data_percept_epoched.trial(1:48);
% data_percept_epoched.time  = data_percept_epoched.time(1:48);
% 
% % Also trim trialinfo if present
% if isfield(data_percept_epoched, 'trialinfo')
%     data_percept_epoched.trialinfo = data_percept_epoched.trialinfo(1:4);
% end


save('lfp_secondhalf2.mat', 'data_percept_epoched');