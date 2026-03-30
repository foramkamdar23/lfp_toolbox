%% =========================================================================
%  EEG–PERCEPT SYNCHRONIZATION PIPELINE
%  TENS artifact alignment using EXG7
%  =========================================================================

%% -------------------------------------------------------------------------
% BLOCK 1 — Initialize MATLAB environment
% -------------------------------------------------------------------------
clear; close all; clc;

% Add FieldTrip
addpath(genpath('C:\Users\saosorio\Toolboxes\fieldtrip-20250106'));
ft_defaults;

%% -------------------------------------------------------------------------
% BLOCK 2 — Load EEG recording
% -------------------------------------------------------------------------
dataDir = 'C:\Users\saosorio\Projects\WorkingMemory_CPEEG';
eegFile = 'cpeeg01_b02_offstimoffmed_wm.bdf';
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
% BLOCK 9 — Load Percept intracranial data
% -------------------------------------------------------------------------
load("cpeeg01_block_02_wm_offmed_offstim_percept.mat");

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
% -------------------------------------------------------------------------

% --- 20a: Read and correct BioSemi STATUS events --------------------------
% BioSemi adds 512 (bit 9) to all STATUS values. Masking to 8 bits isolates
% the trigger byte sent by Psychtoolbox.
events_raw    = ft_read_event(eegFullPath);
status_events = events_raw(strcmp({events_raw.type}, 'STATUS'));

% Handle both numeric and string value fields (FieldTrip is inconsistent)
if ischar(status_events(1).value) || isstring(status_events(1).value)
    raw_values = str2double({status_events.value});
else
    raw_values = double([status_events.value]);
end

corrected_values = bitand(raw_values - 512, 255);

% Print corrected code summary for inspection
fprintf('\nCorrected BioSemi event codes:\n');
fprintf('%-8s  %-8s  %-12s\n', 'Code', 'Count', 'First sample');
fprintf('%s\n', repmat('-', 1, 35));
unique_vals = unique(corrected_values);
for i = 1:length(unique_vals)
    v     = unique_vals(i);
    n     = sum(corrected_values == v);
    first = status_events(find(corrected_values == v, 1)).sample;
    fprintf('%-8d  %-8d  %-12d\n', v, n, first);
end

% --- 20b: Set epoching parameters -----------------------------------------
PRE_STIM_S  = 1.0;
POST_STIM_S = 2.0;
EPOCH_CODES = 141:145;      % per-digit ON codes, position 1..5 (numeric)

% --- 20c: Build trl matrix from corrected events --------------------------
% Col 1: begsample (Percept space)
% Col 2: endsample (Percept space)
% Col 3: offset (pre-stim samples, relative — no lag shift)
% Col 4: digit position (1..5), custom column preserved by ft_redefinetrial

% --- 20c: Build trl matrix from corrected events --------------------------
pre_samples  = round(PRE_STIM_S  * Fs_percept);
post_samples = round(POST_STIM_S * Fs_percept);

% Fs_bdf is the original BDF sampling rate before downsampling
Fs_bdf = data.fsample;      % pull from the original data structure loaded in Block 2

trl_percept = [];

for pos = 1:length(EPOCH_CODES)
    code          = EPOCH_CODES(pos);
    epoch_events  = status_events(corrected_values == code);
    epoch_samples = [epoch_events.sample];

    % Rescale from BDF sample space (e.g. 2048 Hz) to Percept space (250 Hz)
    epoch_samples_ds = round(epoch_samples * (Fs_percept / Fs_bdf));

    for k = 1:length(epoch_samples_ds)
        beg = epoch_samples_ds(k) - pre_samples  + lag_samples;
        en  = epoch_samples_ds(k) + post_samples + lag_samples;
        trl_percept(end+1, :) = [beg, en, -pre_samples, pos]; %#ok<AGROW>
    end
end

% Sort chronologically
[~, sort_idx] = sort(trl_percept(:,1));
trl_percept   = trl_percept(sort_idx, :);

% --- 20d: Remove trials outside Percept data range ------------------------
valid = trl_percept(:,1) >= 1 & trl_percept(:,2) <= nSamples;
if any(~valid)
    warning('%d / %d trial(s) fall outside Percept range and will be removed.', ...
            sum(~valid), size(trl_percept, 1));
end
trl_percept = trl_percept(valid, :);

% --- 20e: Epoch -----------------------------------------------------------
cfg_epoch     = [];
cfg_epoch.trl = trl_percept;

data_percept_epoched = ft_redefinetrial(cfg_epoch, data_percept_ft);

% Store digit position as a named trialinfo field for easy access later
data_percept_epoched.trialinfo = trl_percept(:, 4);   % Nx1 digit position

fprintf('\nEpoched: %d trials  |  [%.1f, %.1f] s  |  %d LFP channels\n', ...
        length(data_percept_epoched.trial), -PRE_STIM_S, POST_STIM_S, nChannels);

%% -------------------------------------------------------------------------
% BLOCK 21 — Visual sanity check: all epochs per channel
% -------------------------------------------------------------------------

time_ax = data_percept_epoched.time{1};
nTrials = length(data_percept_epoched.trial);
colors  = lines(5);

% --- 21a: Identify NaN-contaminated trials --------------------------------
nan_trial_mask = false(nTrials, nChannels);

for tr = 1:nTrials
    for ch = 1:nChannels
        nan_trial_mask(tr, ch) = any(isnan(data_percept_epoched.trial{tr}(ch, :)));
    end
end

% A trial is rejected if ANY channel contains a NaN
nan_any_ch   = any(nan_trial_mask, 2);          % nTrials x 1 logical
n_nan        = sum(nan_any_ch);
pct_nan      = 100 * n_nan / nTrials;

% Per-channel breakdown
fprintf('\n--- NaN trial report ---\n');
fprintf('Total trials:         %d\n', nTrials);
fprintf('Trials with NaN:      %d (%.1f%% of total)\n', n_nan, pct_nan);
fprintf('Clean trials:         %d (%.1f%% of total)\n', nTrials - n_nan, 100 - pct_nan);
fprintf('\nPer-channel breakdown:\n');
fprintf('%-20s  %-15s  %s\n', 'Channel', 'NaN trials', 'Percentage');
fprintf('%s\n', repmat('-', 1, 50));
for ch = 1:nChannels
    n_ch   = sum(nan_trial_mask(:, ch));
    pct_ch = 100 * n_ch / nTrials;
    fprintf('%-20s  %-15d  %.1f%%\n', data_percept_epoched.label{ch}, n_ch, pct_ch);
end

% Per digit position breakdown
fprintf('\nPer digit position breakdown:\n');
fprintf('%-10s  %-15s  %s\n', 'Position', 'NaN trials', 'Percentage');
fprintf('%s\n', repmat('-', 1, 40));
for pos = 1:5
    pos_mask  = data_percept_epoched.trialinfo == pos;
    n_pos     = sum(nan_any_ch & pos_mask);
    tot_pos   = sum(pos_mask);
    fprintf('Digit %-4d  %-15d  %.1f%%\n', pos, n_pos, 100 * n_pos / tot_pos);
end

% --- 21b: Plot ------------------------------------------------------------
figure('Color','w', 'Position', [100 100 1400 900]);

for ch = 1:nChannels

    subplot(2, 3, ch); hold on;

    for tr = 1:nTrials
        pos    = data_percept_epoched.trialinfo(tr);
        signal = data_percept_epoched.trial{tr}(ch, :);

        if nan_trial_mask(tr, ch)
            % NaN trials: thin red line
            plot(time_ax, signal, 'Color', [0.85 0.15 0.15 0.3], ...
                 'LineWidth', 0.6, 'LineStyle', '--');
        else
            plot(time_ax, signal, 'Color', [colors(pos,:), 0.4], ...
                 'LineWidth', 0.8);
        end
    end

    xline(0, 'k--', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(sprintf('%s  |  NaN: %d/%d', ...
          data_percept_epoched.label{ch}, sum(nan_trial_mask(:,ch)), nTrials), ...
          'Interpreter', 'none');
    xlim([time_ax(1) time_ax(end)]);
    grid on; box off;

end

% Legend
ax_leg = axes('Position', [0.92 0.3 0.01 0.3], 'Visible', 'off');
for pos = 1:5
    line(ax_leg, NaN, NaN, 'Color', colors(pos,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('Digit %d', pos));
end
line(ax_leg, NaN, NaN, 'Color', [0.85 0.15 0.15], 'LineWidth', 1.5, ...
     'LineStyle', '--', 'DisplayName', 'NaN trial');
legend(ax_leg, 'show', 'Location', 'west');

sgtitle(sprintf('All epochs — %d trials | %d NaN (%.1f%% rejected)', ...
        nTrials, n_nan, pct_nan), 'FontSize', 13, 'FontWeight', 'bold');

%% -------------------------------------------------------------------------
% BLOCK 22 — Grand average across all trials and digit positions
% -------------------------------------------------------------------------

% --- 22a: Exclude NaN trials ----------------------------------------------
clean_trials = find(~nan_any_ch);

cfg_clean        = [];
cfg_clean.trials = clean_trials;

data_clean = ft_selectdata(cfg_clean, data_percept_epoched);

% --- 22b: Compute grand average (keep trials to compute SEM) -------------
cfg_ga            = [];
cfg_ga.trials     = 'all';
cfg_ga.keeptrials = 'yes';

grandavg = ft_timelockanalysis(cfg_ga, data_clean);

% --- 22c: Baseline correction (−0.5 to 0 s) ------------------------------
cfg_bl          = [];
cfg_bl.baseline = [-0.5 0];

grandavg_bl = ft_timelockbaseline(cfg_bl, grandavg);

% Compute mean and SEM manually from the kept trials (nTrials x nChan x nTime)
avg_signal = squeeze(mean(grandavg_bl.trial, 1));          % nChan x nTime
sem_signal = squeeze(std(grandavg_bl.trial, 0, 1)) ...
             ./ sqrt(length(clean_trials));                 % nChan x nTime

% --- 22d: Plot ------------------------------------------------------------
time_ax = grandavg_bl.time;

figure('Color','w', 'Position', [100 100 1400 900]);

for ch = 1:nChannels

    subplot(2, 3, ch); hold on;

    sem = sem_signal(ch, :);
    avg = avg_signal(ch, :);

    fill([time_ax fliplr(time_ax)], ...
         [avg + sem, fliplr(avg - sem)], ...
         [0.4 0.6 0.9], 'FaceAlpha', 0.25, 'EdgeColor', 'none');

    plot(time_ax, avg, 'Color', [0.15 0.35 0.75], 'LineWidth', 2);

    xlims = [time_ax(1) time_ax(end)];
    yl    = ylim;
    fill([-0.5 0 0 -0.5], [yl(1) yl(1) yl(2) yl(2)], ...
         [0.85 0.85 0.85], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    xline(0, 'k--', 'LineWidth', 1.2);
    xlabel('Time (s)');
    ylabel('\Delta Amplitude (\muV)');
    title(sprintf('%s  |  n = %d', grandavg_bl.label{ch}, length(clean_trials)), ...
          'Interpreter','none');
    xlim(xlims);
    grid on; box off;

end

sgtitle(sprintf('Grand average (baselined −0.5–0 s)  |  %d clean trials (%.1f%%)', ...
        length(clean_trials), 100*length(clean_trials)/nTrials), ...
        'FontSize', 13, 'FontWeight', 'bold');