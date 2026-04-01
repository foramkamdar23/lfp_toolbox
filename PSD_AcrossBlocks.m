%% =========================================================================
%  POWER SPECTRUM — REST vs WM OVERLAY
% =========================================================================

clear; clc;

addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
addpath('C:\Users\fkamdar\Desktop\repos\WorkingMemory_CPEEG\helper_functions')
ft_defaults;

subject     = 'cpeeg02';
task_blocks = struct('rest', {{'b01'}}, 'wm', {{'b01','b02'}});
tasks       = {'rest', 'wm'};
events2plot = 934:938;

dataDir = fullfile('C:\Users\fkamdar\Desktop\repos\WorkingMemory_CPEEG\', subject);

% colors = struct('rest', [0.2 0.4 0.8], 'wm', [0.8 0.2 0.2]);

freq_results = struct();

%% =========================================================================
% LOOP OVER TASKS
% =========================================================================
for ti = 1:numel(tasks)

    task   = tasks{ti};
    blocks = task_blocks.(task);

    all_data_epoched = cell(1, numel(blocks));

    %% -----------------------------------------------------------------------
    % LOOP OVER BLOCKS
    % -----------------------------------------------------------------------
    for bi = 1:numel(blocks)

        block = blocks{bi};
        fprintf('\n==== Processing task: %s | block: %s ====\n', task, block);

        % --- Load EEG ---
        eegFile     = sprintf('%s_%s_offstimoffmed_%s.bdf', subject, block, task);
        eegFullPath = fullfile(dataDir, eegFile);

        cfg            = [];
        cfg.dataset    = eegFullPath;
        cfg.continuous = 'yes';
        data           = ft_preprocessing(cfg);

        % --- Select channels ---
        cfg         = [];
        cfg.channel = ft_channelselection({'EEG*','EXG*'}, data.label);
        data_clean  = ft_selectdata(cfg, data);

        % --- Downsample ---
        cfg           = [];
        cfg.resamplefs = 250;
        cfg.detrend   = 'no';
        data_ds       = ft_resampledata(cfg, data_clean);
        Fs_eeg        = data_ds.fsample;

        % --- EXG7 envelope ---
        cfg         = [];
        cfg.channel = 'EXG7';
        data_exg    = ft_selectdata(cfg, data_ds);

        cfg            = [];
        cfg.bpfilter   = 'yes'; cfg.bpfreq = [75 85]; cfg.bpfiltord = 4;
        cfg.hilbert    = 'abs';
        data_exg_filt  = ft_preprocessing(cfg, data_exg);

        signal_eeg_env = zscore(data_exg_filt.trial{1}(1,:));

        % --- Load Percept ---
        load(fullfile(dataDir, sprintf('%s_%s_offmedoffstim_%s_percept.mat', subject, block, task)));
        signal_orig  = processed_data(1).signal;
        Fs_percept   = processed_data(1).fs;

        % --- NaN interpolation ---
        signal_sync        = signal_orig;
        nanIdx             = isnan(signal_sync);
        signal_sync(nanIdx)= interp1(find(~nanIdx), signal_sync(~nanIdx), find(nanIdx), 'linear', 'extrap');

        % --- Compute lag ---
        [b,a]               = butter(4, [75 85]/(Fs_percept/2), 'bandpass');
        signal_percept_env  = zscore(abs(filtfilt(b, a, signal_sync)));
        [xc, lags]          = xcorr(signal_percept_env, signal_eeg_env);
        [~, idx]            = max(xc);
        lag_samples         = lags(idx);
        lag_abs             = abs(lag_samples);

        % --- Build FT structure ---
        t_eeg      = (0:length(signal_eeg_env)-1) / Fs_eeg;
        nChannels  = length(processed_data);
        nSamples   = length(t_eeg);

        percept_aligned_all = nan(nChannels, nSamples);

        for ch = 1:nChannels
            sig = processed_data(ch).signal;
            if lag_samples > 0
                sig = [sig; nan(lag_abs,1)];
            elseif lag_samples < 0
                sig = [nan(lag_abs,1); sig];
            end
            sig = sig(:)';
            len = length(sig);
            if len >= nSamples
                percept_aligned_all(ch,:) = sig(1:nSamples);
            else
                percept_aligned_all(ch,1:len) = sig;
            end
        end

        data_percept_ft            = [];
        data_percept_ft.label      = {processed_data.channel}';
        data_percept_ft.fsample    = Fs_percept;
        data_percept_ft.trial      = {percept_aligned_all};
        data_percept_ft.time       = {t_eeg};
        data_percept_ft.sampleinfo = [1 nSamples];
        data_percept_ft            = ft_checkdata(data_percept_ft, 'datatype', 'raw');

        % --- Epoching ---
        PRE_STIM_S   = .5;
        POST_STIM_S  = 1;
        pre_samples  = round(PRE_STIM_S  * Fs_percept);
        post_samples = round(POST_STIM_S * Fs_percept);
        epoch_len    = pre_samples + post_samples;
        trl          = [];

        if strcmp(task, 'rest')

            artifact_env    = data_exg_filt.trial{1}(1,:);
            artifact_thresh = prctile(artifact_env, 91);
            artifact_idx    = artifact_env > artifact_thresh;
            artifact_samples= find(artifact_idx);

            diffs           = diff(artifact_samples);
            [~, max_gap_idx]= max(diffs);
            clean_start     = artifact_samples(max_gap_idx) + 1;
            clean_end       = artifact_samples(max_gap_idx + 1) - 1;

            fprintf('Clean window: %.1f s to %.1f s\n', clean_start/Fs_percept, clean_end/Fs_percept);

            step      = round(epoch_len / 2);
            start_idx = clean_start;
            trial_id  = 1;

            while (start_idx + epoch_len - 1) <= clean_end
                beg = start_idx;
                en  = start_idx + epoch_len - 1;
                trl(end+1,:) = [beg en -pre_samples trial_id];
                start_idx = start_idx + step;
                trial_id  = trial_id + 1;
            end

            fprintf('Created %d dummy epochs.\n', size(trl,1));
            
            % get diagnostic plot to make sure epochs are artifact free
            plot_artifact_epochs(artifact_env, artifact_thresh, trl, ...
                clean_start, clean_end, subject, block, Fs_percept);
        else

            events_raw       = ft_read_event(eegFullPath);
            status_events    = events_raw(strcmp({events_raw.type}, 'STATUS'));
            corrected_values = check_and_correct_biosemi(status_events, events2plot(1));
            Fs_bdf           = data.fsample;

            for pos = 1:length(events2plot)
                code = events2plot(pos);
                ev   = status_events(corrected_values == code);
                samp = round([ev.sample] * (Fs_percept/Fs_bdf));
                for k = 1:length(samp)
                    beg = samp(k) - pre_samples + lag_samples;
                    en  = samp(k) + post_samples + lag_samples;
                    trl(end+1,:) = [beg en -pre_samples pos];
                end
            end

            valid = trl(:,1) >= 1 & trl(:,2) <= nSamples;
            trl   = trl(valid,:);

        end

        cfg     = [];
        cfg.trl = trl;
        data_epoched          = ft_redefinetrial(cfg, data_percept_ft);
        data_epoched.trialinfo= trl(:,4);

        all_data_epoched{bi} = data_epoched;

    end % block loop

    %% -----------------------------------------------------------------------
    % CONCATENATE BLOCKS
    % -----------------------------------------------------------------------
    cfg = [];
    if numel(all_data_epoched) == 1
        data_all = all_data_epoched{1};
    else
        data_all = ft_appenddata(cfg, all_data_epoched{:});
    end

    % --- NaN rejection ---
    nTrials   = length(data_all.trial);
    nan_mask  = false(nTrials, nChannels);
    for tr = 1:nTrials
        for ch = 1:nChannels
            nan_mask(tr,ch) = any(isnan(data_all.trial{tr}(ch,:)));
        end
    end
    clean_trials = find(~any(nan_mask, 2));

    cfg        = [];
    cfg.trials = clean_trials;
    data_clean_task = ft_selectdata(cfg, data_all);

    % --- Power spectrum ---
    cfg            = [];
    cfg.method     = 'mtmfft';
    cfg.output     = 'pow';
    cfg.taper      = 'dpss';
    cfg.foilim     = [1 100];
    cfg.tapsmofrq  = 1.5;
    cfg.keeptrials = 'no';

    freq_results.(task) = ft_freqanalysis(cfg, data_clean_task);

    fprintf('\nDone with task: %s\n', task);

end % task loop

%% =========================================================================
% PLOT — OVERLAY REST vs WM PER CONTACT
% =========================================================================
chan_labels = freq_results.rest.label;
nChannels   = numel(chan_labels);

figure('Color','w', 'Position',[100 100 1200 800]);
colors = struct('rest', [0.00 0.62 0.65], 'wm', [0.84 0.37 0.00]);

for ch = 1:nChannels

    subplot(2, 3, ch); hold on;

    for ti = 1:numel(tasks)
        task = tasks{ti};
        pow  = squeeze(freq_results.(task).powspctrm(ch,:));
        plot(freq_results.(task).freq, pow, ...
            'Color', colors.(task), 'LineWidth', 1.8);
    end

    set(gca, 'YScale', 'log', 'FontSize', 14);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(chan_labels{ch}, 'Interpreter', 'none', 'FontWeight', 'normal');
    xlim([1 100]);

    legend({'Rest', 'Working Memory'}, 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');

end

sgtitle(sprintf('Power Spectra per contact — %s', subject), ...
    'FontSize', 18, 'FontWeight', 'normal');