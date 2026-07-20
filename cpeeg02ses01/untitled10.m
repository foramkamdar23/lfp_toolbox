%% Sanity check: use plot_tfr_grid for power, then add stats as a 3rd panel

chan_name = plot_layout{1};

% --- Averaged TFR per condition, full time range including baseline ---
cfg = [];
cfg.trials = feel_neg_idx;
cfg.avgoverrpt = 'yes';
freq_avg_feelneg_full = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = feel_neu_idx;
cfg.avgoverrpt = 'yes';
freq_avg_feelneu_full = ft_selectdata(cfg, freq_clean_bl);

pow_clim = [-2 2];

% Power panels - reuse plot_tfr_grid (1x1 layout, single channel each)
plot_tfr_grid(freq_avg_feelneg_full, {chan_name}, pow_clim, ...
              sprintf('%s - Emotion Perception Condition: Negative', label_map(chan_name)), label_map);

plot_tfr_grid(freq_avg_feelneu_full, {chan_name}, pow_clim, ...
              sprintf('%s - Emotion Perception Condition: Neutral', label_map(chan_name)), label_map)

%% Stats panel - separate figure for now
figure('Name', sprintf('%s: Negative vs Neutral t-stat', label_map(chan_name)));
imagesc(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
        stat_freq_feelneg_vs_feelneu_ZTR.freq, ...
        squeeze(stat_freq_feelneg_vs_feelneu_ZTR.stat));
axis xy;
xline(0, 'k--');
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim_val = max(abs(squeeze(stat_freq_feelneg_vs_feelneu_ZTR.stat)), [], 'all');
clim([-clim_val clim_val]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Passive viewing: Negative vs. Neutral (t-statistic)');

hold on;
sig_alpha = 0.025;
if ~isempty(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
    sig_k = find([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelneg_vs_feelneu_ZTR.posclusterslabelmat == k);
        contour(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
                log10(stat_freq_feelneg_vs_feelneu_ZTR.freq), ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
if ~isempty(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)
    sig_k = find([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelneg_vs_feelneu_ZTR.negclusterslabelmat == k);
        contour(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
                log10(stat_freq_feelneg_vs_feelneu_ZTR.freq), ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
hold off;