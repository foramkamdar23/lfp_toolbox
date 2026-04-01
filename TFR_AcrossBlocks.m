%% =========================================================================
%  EEG–PERCEPT SYNCHRONIZATION PIPELINE (MULTI-BLOCK)
% =========================================================================

clear; clc;

addpath(genpath('C:\Users\saosorio\Toolboxes\fieldtrip-20250106'));
addpath('C:\Users\saosorio\Projects\WorkingMemory_CPEEG\helper_functions')
ft_defaults;

subject = 'cpeeg02';
task    = 'rest';
blocks  = {'b01'};
events2plot = 934:938;

dataDir = fullfile('C:\Users\saosorio\Projects\WorkingMemory_CPEEG',subject);

all_data_epoched = cell(1,numel(blocks));   % ✅ NEW

%% =========================================================================
% LOOP OVER BLOCKS
% =========================================================================
for bi = 1:numel(blocks)

    block = blocks{bi};
    fprintf('\n================ Processing %s ================\n', block);

    %% -------------------------------------------------------------------------
    % BLOCK 2 — Load EEG recording
    % -------------------------------------------------------------------------
    eegFile = sprintf('%s_%s_offstimoffmed_%s.bdf', subject, block, task);
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
    % BLOCK 4 — Downsample EEG
    % -------------------------------------------------------------------------
    cfg = [];
    cfg.resamplefs = 250;
    cfg.detrend = 'no';
    data_ds = ft_resampledata(cfg, data_clean);
    Fs_eeg = data_ds.fsample;

    %% -------------------------------------------------------------------------
    % BLOCK 5–14 (UNCHANGED)
    % -------------------------------------------------------------------------
    cfg = []; cfg.channel = 'EXG7';
    data_exg = ft_selectdata(cfg, data_ds);

    cfg = [];
    cfg.bpfilter='yes'; cfg.bpfreq=[75 85]; cfg.bpfiltord=4;
    cfg.hilbert='abs';
    data_exg_filt = ft_preprocessing(cfg, data_exg);

    signal_eeg_env = zscore(data_exg_filt.trial{1}(1,:));

    load(fullfile(dataDir, sprintf("%s_%s_offmedoffstim_%s_percept.mat", subject, block, task)));
    signal_orig = processed_data(1).signal;
    Fs_percept  = processed_data(1).fs;

    signal_sync = signal_orig;
    nanIdx = isnan(signal_sync);
    signal_sync(nanIdx) = interp1(find(~nanIdx),signal_sync(~nanIdx),find(nanIdx),'linear','extrap');

    [b,a] = butter(4,[75 85]/(Fs_percept/2),'bandpass');
    signal_percept_env = zscore(abs(filtfilt(b,a,signal_sync)));

    [xc,lags] = xcorr(signal_percept_env,signal_eeg_env);
    [~,idx] = max(xc);
    lag_samples = lags(idx);

    %% -------------------------------------------------------------------------
    % BLOCK 16 — Align
    % -------------------------------------------------------------------------
    lag_abs = abs(lag_samples);

    if lag_samples > 0
        signal_percept_aligned = [signal_orig; nan(lag_abs,1)];
    elseif lag_samples < 0
        signal_percept_aligned = [nan(lag_abs,1); signal_orig];
    else
        signal_percept_aligned = signal_orig;
    end

    t_eeg = (0:length(signal_eeg_env)-1)/Fs_eeg;

    %% -------------------------------------------------------------------------
    % BLOCK 19 — Build FT structure
    % -------------------------------------------------------------------------
    nChannels = length(processed_data);
    nSamples  = length(t_eeg);

    percept_aligned_all = nan(nChannels,nSamples);

    for ch = 1:nChannels

        sig = processed_data(ch).signal;

        % --- Apply lag ---
        if lag_samples > 0
            sig = [sig; nan(lag_abs,1)];
        elseif lag_samples < 0
            sig = [nan(lag_abs,1); sig];
        end

        % --- FORCE ROW VECTOR ---
        sig = sig(:)';

        % --- Trim / pad ---
        len = length(sig);

        if len >= nSamples
            percept_aligned_all(ch,:) = sig(1:nSamples);
        else
            percept_aligned_all(ch,1:len) = sig;
        end

    end

    data_percept_ft = [];
    data_percept_ft.label = {processed_data.channel}';
    data_percept_ft.fsample = Fs_percept;
    data_percept_ft.trial = {percept_aligned_all};
    data_percept_ft.time  = {t_eeg};
    data_percept_ft.sampleinfo = [1 nSamples];

    data_percept_ft = ft_checkdata(data_percept_ft,'datatype','raw');

    %% -------------------------------------------------------------------------
    % BLOCK 20 — Epoch
    % -------------------------------------------------------------------------

    PRE_STIM_S  = 2.5;
    POST_STIM_S = 2.5;

    pre_samples  = round(PRE_STIM_S  * Fs_percept);
    post_samples = round(POST_STIM_S * Fs_percept);
    epoch_len    = pre_samples + post_samples;

    trl = [];

    if strcmp(task,'rest')

        % -------------------------------------------------
        % DETECT SYNC ARTIFACT PERIODS
        % -------------------------------------------------
        artifact_env = data_exg_filt.trial{1}(1,:);
        % -------------------------------------------------
        % FIND CLEAN WINDOW BOUNDARIES
        % -------------------------------------------------
        artifact_thresh = prctile(artifact_env, 91);
        artifact_idx = artifact_env > artifact_thresh;

        % Find all samples where artifact is active
        artifact_samples = find(artifact_idx);

        diffs = diff(artifact_samples);
        [~, max_gap_idx] = max(diffs);

        % Clean window is strictly between the two clusters
        clean_start = artifact_samples(max_gap_idx) + 1;
        clean_end   = artifact_samples(max_gap_idx + 1) - 1;

        fprintf('Clean window: sample %d to %d (%.1f s to %.1f s)\n', ...
            clean_start, clean_end, clean_start/Fs_percept, clean_end/Fs_percept);
        % -------------------------------------------------
        % CREATE DUMMY EPOCHS WITHIN CLEAN WINDOW ONLY
        % -------------------------------------------------
        step     = round(epoch_len / 2);
        start_idx = clean_start;
        trial_id  = 1;

        while (start_idx + epoch_len - 1) <= clean_end

            beg = start_idx;
            en  = start_idx + epoch_len - 1;

            offset = -pre_samples;
            trl(end+1,:) = [beg en offset trial_id];

            start_idx = start_idx + step;
            trial_id  = trial_id + 1;
        end

        fprintf('Created %d dummy epochs within clean window.\n', size(trl,1));

    else

        % -------------------------------------------------
        % ORIGINAL EVENT-BASED EPOCHING
        % -------------------------------------------------

        events_raw    = ft_read_event(eegFullPath);
        status_events = events_raw(strcmp({events_raw.type}, 'STATUS'));

        corrected_values = check_and_correct_biosemi(status_events, events2plot(1));

        Fs_bdf = data.fsample;

        for pos = 1:length(events2plot)
            code = events2plot(pos);
            ev   = status_events(corrected_values == code);
            samp = round([ev.sample] * (Fs_percept/Fs_bdf));

            for k = 1:length(samp)
                beg = samp(k)-pre_samples+lag_samples;
                en  = samp(k)+post_samples+lag_samples;
                trl(end+1,:) = [beg en -pre_samples pos];
            end
        end

        valid = trl(:,1)>=1 & trl(:,2)<=nSamples;
        trl = trl(valid,:);

    end

    % --- Apply to FieldTrip ---
    cfg = [];
    cfg.trl = trl;

    data_epoched = ft_redefinetrial(cfg, data_percept_ft);
    data_epoched.trialinfo = trl(:,4);

end % ===== END LOOP =====

%% =========================================================================
% CONCATENATE BLOCKS
% =========================================================================
cfg = [];

if isscalar(data_epoched)
    % No concatenation needed
    data_all = data_epoched;
else
    data_all = ft_appenddata(cfg, data_epoched{:});
end

fprintf('\nTotal trials after concatenation: %d\n', length(data_all.trial));

%% =========================================================================
% NaN rejection (GLOBAL)
% =========================================================================
nTrials = length(data_all.trial);
nChannels = length(data_all.label);

nan_mask = false(nTrials,nChannels);

for tr = 1:nTrials
    for ch = 1:nChannels
        nan_mask(tr,ch) = any(isnan(data_all.trial{tr}(ch,:)));
    end
end

nan_any = any(nan_mask,2);
clean_trials = find(~nan_any);

cfg = [];
cfg.trials = clean_trials;
data_clean = ft_selectdata(cfg, data_all);

%% =========================================================================
% GRAND AVERAGE
% =========================================================================
cfg = [];
cfg.keeptrials = 'yes';

grandavg = ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.baseline = [-0.5 0];
grandavg_bl = ft_timelockbaseline(cfg, grandavg);

avg_signal = squeeze(mean(grandavg_bl.trial,1));
sem_signal = squeeze(std(grandavg_bl.trial,0,1))/sqrt(length(clean_trials));

%% =========================================================================
% POWER SPECTRUM (FFT)
% =========================================================================

cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'pow';
cfg.taper     = 'dpss';
cfg.foilim    = [1 100];       % frequency range
cfg.tapsmofrq =  2;             % smoothing in Hz (adjust for smoother lines)
cfg.keeptrials = 'no';         % average across trials

freq_data = ft_freqanalysis(cfg, data_clean);

figure('Color','w', 'Position',[100 100 1200 800]);

nChannels = length(freq_data.label);

for ch = 1:nChannels

    subplot(2,3,ch); hold on;

    pow = squeeze(freq_data.powspctrm(ch,:));

    plot(freq_data.freq, pow, 'LineWidth', 1.5);

    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(freq_data.label{ch}, 'Interpreter','none');

    xlim([1 100]);

    % Optional: log scale (VERY common for LFP/EEG)
    set(gca, 'YScale', 'log', 'FontSize', 16);

end

sgtitle(sprintf('Power Spectra per contact for %s in %s', subject, task), ...
     'FontSize', 18, 'FontWeight', 'normal');
%% =========================================================================
% TFR (COMBINED DATA)
% =========================================================================
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.foi = 4:2:100;
cfg.toi = -2:0.01:2;
cfg.width = 7;
cfg.keeptrials = 'yes';

tfr_data = ft_freqanalysis(cfg, data_all);

%%
figure('Color','w', 'Position',[100 100 1400 900]);

% Time and frequency vectors
time_ax = tfr_data.time;        % 1 x nTime
freq_ax = tfr_data.freq;        % 1 x nFreq

baseline_win = [-0.5 0];        % seconds for baseline

for ch = 1:nChannels
    subplot(2,3,ch); hold on;

    % Extract power: nFreq x nTime x nTrials
    pow_trials = squeeze(tfr_data.powspctrm(:, ch, :, :));
    pow_mean = squeeze(nanmean(pow_trials, 1));

    % Ensure baseline indices are valid
    baseline_idx = find(time_ax >= baseline_win(1) & time_ax <= baseline_win(2));

    % Clip indices to the available time points
    baseline_idx = baseline_idx(baseline_idx <= size(pow_mean,2));

    % Fallback if baseline window does not overlap
    if isempty(baseline_idx)
        warning('Baseline window [%g %g]s does not overlap TFR time points. Using first time point as baseline.', ...
            baseline_win(1), baseline_win(2));
        baseline_idx = 1;
    end

    % Compute baseline mean (per frequency)
    baseline_mean = mean(pow_mean(:, baseline_idx), 2);  % nFreq x 1

    % Relative change from baseline
    pow_rel = (pow_mean - baseline_mean) ./ baseline_mean;

    % Plot TFR
    imagesc(time_ax, freq_ax, pow_rel, [-.7 .7]);
    axis xy;                % flip y-axis so low freq is at bottom
    colormap jet;
    colorbar;
    % caxis([-1 1]);         % optional: normalize scale

    % Plot formatting
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(tfr_data.label{ch}, 'Interpreter','none', 'FontWeight','normal');
    xlim([-0.5 1.5]);
    ylim([4 100]);

    % Overlay vertical line at stimulus onset
    plot([0 0], ylim, 'k--', 'LineWidth', 1.2);
    set(gca,'FontSize',16)

end

if strcmp(task,'rest')
    task = task;
else
    task = 'Working Memory';
end

sgtitle(sprintf('TFR per contact for %s in %s', subject, task), ...
     'FontSize', 18, 'FontWeight', 'normal');