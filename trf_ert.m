%% =========================================================================
% TFR — FEEL vs TONE vs REST (FINAL)
% =========================================================================

% --- PARAMETERS ---
cfg = [];
cfg.method       = 'wavelet';
cfg.output       = 'pow';
cfg.foi          = 4:2:100;      % frequencies
cfg.toi          = -1:0.02:5;    % time window (match your epochs)
cfg.width        = 7;
cfg.keeptrials   = 'yes';

% --- COMPUTE TFR ---
fprintf('Computing TFR...\n');

tfr_feel = ft_freqanalysis(cfg, data_feel);
tfr_tone = ft_freqanalysis(cfg, data_tone);
tfr_rest = ft_freqanalysis(cfg, data_rest);

fprintf('TFR done\n');

%% =========================================================================
% BASELINE + PLOTTING
% =========================================================================

baseline_win = [-0.5 0];

nChannels = length(tfr_feel.label);

figure('Color','w','Position',[100 100 1400 900]);

for ch = 1:nChannels

    % =========================
    % FEEL
    % =========================
    subplot(3, nChannels, ch); hold on;

    pow = squeeze(nanmean(tfr_feel.powspctrm(:, ch, :, :),1)); % freq x time

    baseline_idx = find(tfr_feel.time >= baseline_win(1) & ...
                        tfr_feel.time <= baseline_win(2));

    if isempty(baseline_idx)
        baseline_idx = 1;
    end

    baseline_mean = mean(pow(:,baseline_idx),2);
    pow_rel = (pow - baseline_mean) ./ baseline_mean;

    imagesc(tfr_feel.time, tfr_feel.freq, pow_rel, [-0.7 0.7]);
    axis xy; colormap jet; colorbar;

    title(['FEEL - ' tfr_feel.label{ch}], 'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Freq (Hz)');
    xlim([-0.5 2]);
    ylim([4 100]);
    xline(0,'k--','LineWidth',1.2);

    % =========================
    % TONE
    % =========================
    subplot(3, nChannels, ch + nChannels); hold on;

    pow = squeeze(nanmean(tfr_tone.powspctrm(:, ch, :, :),1));

    baseline_mean = mean(pow(:,baseline_idx),2);
    pow_rel = (pow - baseline_mean) ./ baseline_mean;

    imagesc(tfr_tone.time, tfr_tone.freq, pow_rel, [-0.7 0.7]);
    axis xy; colormap jet; colorbar;

    title(['TONE - ' tfr_tone.label{ch}], 'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Freq (Hz)');
    xlim([-0.5 2]);
    ylim([4 100]);
    xline(0,'k--','LineWidth',1.2);

    % =========================
    % REST
    % =========================
    subplot(3, nChannels, ch + 2*nChannels); hold on;

    pow = squeeze(nanmean(tfr_rest.powspctrm(:, ch, :, :),1));

    baseline_idx_rest = find(tfr_rest.time >= baseline_win(1) & ...
                             tfr_rest.time <= baseline_win(2));

    if isempty(baseline_idx_rest)
        baseline_idx_rest = 1;
    end

    baseline_mean = mean(pow(:,baseline_idx_rest),2);
    pow_rel = (pow - baseline_mean) ./ baseline_mean;

    imagesc(tfr_rest.time, tfr_rest.freq, pow_rel, [-0.7 0.7]);
    axis xy; colormap jet; colorbar;

    title(['REST - ' tfr_rest.label{ch}], 'Interpreter','none');
    xlabel('Time (s)');
    ylabel('Freq (Hz)');
    xlim([-0.5 2]);
    ylim([4 100]);
    xline(0,'k--','LineWidth',1.2);

end

sgtitle('Time-Frequency Representation (Baseline corrected)', ...
    'FontSize', 14, 'FontWeight', 'bold');