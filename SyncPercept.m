
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



%% -------------------------------------------------------------------------
% BLOCK 22 — Grand average (LFP)
% -------------------------------------------------------------------------

%% --- 22a: Remove NaN trials ---------------------------------------------
clean_trials = find(~nan_any_ch);

cfg = [];
cfg.trials = clean_trials;

data_clean = ft_selectdata(cfg, data_percept_epoched);

fprintf('Using %d clean trials (%.1f%%)\n', ...
    length(clean_trials), 100*length(clean_trials)/length(nan_any_ch));

%% --- 22b: Timelock (keep trials for SEM) --------------------------------
cfg = [];
cfg.keeptrials = 'yes';

grandavg = ft_timelockanalysis(cfg, data_clean);

%% --- 22c: Baseline correction -------------------------------------------
cfg = [];
cfg.baseline = [-0.5 0];   % adjust if needed

grandavg_bl = ft_timelockbaseline(cfg, grandavg);

%% --- 22d: Compute mean + SEM --------------------------------------------
avg_signal = squeeze(mean(grandavg_bl.trial, 1));   % [chan x time]
sem_signal = squeeze(std(grandavg_bl.trial, 0, 1)) ...
             ./ sqrt(size(grandavg_bl.trial,1));

time_ax   = grandavg_bl.time;
nChannels = length(grandavg_bl.label);

%% --- 22e: Plot -----------------------------------------------------------
figure('Color','w', 'Position', [100 100 1400 900]);

for ch = 1:nChannels

    subplot(ceil(nChannels/2), 2, ch); hold on;

    avg = avg_signal(ch,:);
    sem = sem_signal(ch,:);

    % SEM shading
    fill([time_ax fliplr(time_ax)], ...
         [avg+sem fliplr(avg-sem)], ...
         [0.4 0.6 0.9], ...
         'FaceAlpha', 0.25, 'EdgeColor', 'none');

    % Mean signal
    plot(time_ax, avg, 'Color', [0.15 0.35 0.75], 'LineWidth', 2);

    % Baseline window shading
    yl = ylim;
    fill([-0.5 0 0 -0.5], [yl(1) yl(1) yl(2) yl(2)], ...
         [0.85 0.85 0.85], ...
         'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Stim onset
    xline(0, 'k--', 'LineWidth', 1.2);

    xlabel('Time (s)');
    ylabel('\Delta Amplitude (\muV)');

    title(sprintf('%s | n = %d', ...
        grandavg_bl.label{ch}, length(clean_trials)), ...
        'Interpreter','none');

    xlim([time_ax(1) time_ax(end)]);
    grid on; box off;

end

sgtitle(sprintf('Grand average (baseline -0.5–0 s) | %d trials', ...
    length(clean_trials)), ...
    'FontSize', 13, 'FontWeight', 'bold');