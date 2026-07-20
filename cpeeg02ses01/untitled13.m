%% plot neg

%chan_name = plot_layout{1};

% --- Averaged TFR per condition, full time range including baseline ---
cfg = [];
cfg.trials = neg_idx;
cfg.avgoverrpt = 'yes';
freq_avg_neg_full = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = neu_idx;
cfg.avgoverrpt = 'yes';
freq_avg_neu_full = ft_selectdata(cfg, freq_clean_bl);

pow_clim = [-2 2];

% Power panels - reuse plot_tfr_grid (1x1 layout, single channel each)
plot_tfr_grid(freq_avg_neg_full, {chan_name}, pow_clim, ...
              sprintf('%s Negative', label_map(chan_name)), label_map);

plot_tfr_grid(freq_avg_neu_full, {chan_name}, pow_clim, ...
              sprintf('%s Neutral', label_map(chan_name)), label_map)

% Stats panel - separate figure for now
figure('Name', sprintf('%s: Negative vs Neutral t-stat', label_map(chan_name)));
imagesc(stat_freq_neg_vs_neu.time, stat_freq_neg_vs_neu.freq,  squeeze(stat_freq_neg_vs_neu.stat));
axis xy;
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim(pow_clim);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Negative vs. Neutral (t-statistic)');
yticks([2 4 5 6 7 8 16 28]);

%
hold on;

% --- Largest positive cluster (regardless of significance) ---
if ~isempty(stat_freq_neg_vs_neu.posclusters)
    pos_sizes = zeros(1, numel(stat_freq_neg_vs_neu.posclusters));
    for k = 1:numel(stat_freq_neg_vs_neu.posclusters)
        mask = squeeze(stat_freq_neg_vs_neu.posclusterslabelmat == k);
        pos_sizes(k) = sum(mask(:));
    end
    [~, biggest_pos_k] = max(pos_sizes);
    mask = squeeze(stat_freq_neg_vs_neu.posclusterslabelmat == biggest_pos_k);
    contour(stat_freq_neg_vs_neu.time, ...
            log10(stat_freq_neg_vs_neu.freq), ...
            double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
end

% --- Largest negative cluster (regardless of significance) ---
if ~isempty(stat_freq_neg_vs_neu.negclusters)
    neg_sizes = zeros(1, numel(stat_freq_neg_vs_neu.negclusters));
    for k = 1:numel(stat_freq_neg_vs_neu.negclusters)
        mask = squeeze(stat_freq_neg_vs_neu.negclusterslabelmat == k);
        neg_sizes(k) = sum(mask(:));
    end
    [~, biggest_neg_k] = max(neg_sizes);
    mask = squeeze(stat_freq_neg_vs_neu.negclusterslabelmat == biggest_neg_k);
    contour(stat_freq_neg_vs_neu.time, ...
           log10(stat_freq_neg_vs_neu.freq), ...
            double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
end

hold off;

%% plot pos

%chan_name = plot_layout{1};

% --- Averaged TFR per condition, full time range including baseline ---
cfg = [];
cfg.trials = pos_idx;
cfg.avgoverrpt = 'yes';
freq_avg_pos_full = ft_selectdata(cfg, freq_clean_bl);

cfg = [];
cfg.trials = neu_idx;
cfg.avgoverrpt = 'yes';
freq_avg_neu_full = ft_selectdata(cfg, freq_clean_bl);

pow_clim = [-2 2];

% Power panels - reuse plot_tfr_grid (1x1 layout, single channel each)
plot_tfr_grid(freq_avg_pos_full, {chan_name}, pow_clim, ...
              sprintf('%s Positive', label_map(chan_name)), label_map);

plot_tfr_grid(freq_avg_neu_full, {chan_name}, pow_clim, ...
              sprintf('%s Neutral', label_map(chan_name)), label_map)

% Stats panel - separate figure for now
figure('Name', sprintf('%s: Positive vs Neutral t-stat 2', label_map(chan_name)));
imagesc(stat_freq_pos_vs_neu.time, stat_freq_pos_vs_neu.freq,  squeeze(stat_freq_pos_vs_neu.stat));
axis xy;
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim(pow_clim);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Positive vs. Neutral (t-statistic)');
yticks([2 4 5 6 7 8 16 28]);


%
hold on;

% --- Largest positive cluster (regardless of significance) ---
if ~isempty(stat_freq_pos_vs_neu.posclusters)
    pos_sizes = zeros(1, numel(stat_freq_pos_vs_neu.posclusters));
    for k = 1:numel(stat_freq_pos_vs_neu.posclusters)
        mask = squeeze(stat_freq_pos_vs_neu.posclusterslabelmat == k);
        pos_sizes(k) = sum(mask(:));
    end
    [~, biggest_pos_k] = max(pos_sizes);
    mask = squeeze(stat_freq_pos_vs_neu.posclusterslabelmat == biggest_pos_k);
    contour(stat_freq_pos_vs_neu.time, ...
            stat_freq_pos_vs_neu.freq, ...
            double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
end

% --- Largest negative cluster (regardless of significance) ---
if ~isempty(stat_freq_pos_vs_neu.negclusters)
    neg_sizes = zeros(1, numel(stat_freq_pos_vs_neu.negclusters));
    for k = 1:numel(stat_freq_pos_vs_neu.negclusters)
        mask = squeeze(stat_freq_pos_vs_neu.negclusterslabelmat == k);
        neg_sizes(k) = sum(mask(:));
    end
    [~, biggest_neg_k] = max(neg_sizes);
    mask = squeeze(stat_freq_pos_vs_neu.negclusterslabelmat == biggest_neg_k);
    contour(stat_freq_pos_vs_neu.time, ...
            stat_freq_pos_vs_neu.freq, ...
            double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
end

hold off;

%% USING PCOLOR
% Stats panel - separate figure for now
figure('Name', sprintf('%s: Positive vs Neutral t-stat', label_map(chan_name)));

t = stat_freq_pos_vs_neu.time;
f = stat_freq_pos_vs_neu.freq;
tstat = squeeze(stat_freq_pos_vs_neu.stat);

pcolor(t, f, tstat);
shading interp;
set(gca, 'YScale', 'log');
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim(pow_clim);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Positive vs. Neutral (t-statistic)');
yticks([2 4 5 6 7 8 16 28]);

hold on;
line([0.808 1.316], [5 5], 'Color', 'w', 'LineWidth', 4);
line([0.604 1.696], [6 6], 'Color', 'w', 'LineWidth', 4);
hold off;

%% USING PCOLOR
% Stats panel - separate figure for now
figure('Name', sprintf('%s: Positive vs Neutral t-stat', label_map(chan_name)));

t = stat_freq_neg_vs_neu.time;
f = stat_freq_neg_vs_neu.freq;
tstat = squeeze(stat_freq_neg_vs_neu.stat);

pcolor(t, f, tstat);
shading interp;
set(gca, 'YScale', 'log');
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim(pow_clim);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Negative vs. Neutral (t-statistic)');
yticks([2 4 5 6 7 8 16 28]);

hold on;

% --- Largest positive cluster (regardless of significance) ---
if ~isempty(stat_freq_neg_vs_neu.posclusters)
    pos_sizes = zeros(1, numel(stat_freq_neg_vs_neu.posclusters));
    for k = 1:numel(stat_freq_neg_vs_neu.posclusters)
        mask = squeeze(stat_freq_neg_vs_neu.posclusterslabelmat == k);
        pos_sizes(k) = sum(mask(:));
    end
    [~, biggest_pos_k] = max(pos_sizes);
    mask = squeeze(stat_freq_neg_vs_neu.posclusterslabelmat == biggest_pos_k);
    contour(stat_freq_neg_vs_neu.time, ...
            stat_freq_neg_vs_neu.freq, ...
            double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
end

% --- Largest negative cluster (regardless of significance) ---
if ~isempty(stat_freq_neg_vs_neu.negclusters)
    neg_sizes = zeros(1, numel(stat_freq_neg_vs_neu.negclusters));
    for k = 1:numel(stat_freq_neg_vs_neu.negclusters)
        mask = squeeze(stat_freq_neg_vs_neu.negclusterslabelmat == k);
        neg_sizes(k) = sum(mask(:));
    end
    [~, biggest_neg_k] = max(neg_sizes);
    mask = squeeze(stat_freq_neg_vs_neu.negclusterslabelmat == biggest_neg_k);
    contour(stat_freq_neg_vs_neu.time, ...
           stat_freq_neg_vs_neu.freq, ...
            double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
end

hold off;
