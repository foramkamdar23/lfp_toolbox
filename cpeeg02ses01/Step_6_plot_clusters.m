%% Plot stats as TFR: t-values with cluster outline (Pos vs Neu, Neg vs Neu)

sig_alpha = 0.025;

figure('Name', sprintf('%s: Cluster stats TFR', label_map(chan_name)));

% --- Pos vs Neu ---
subplot(1,2,1)
imagesc(stat_freq_feelpos_vs_feelneu_ZTR.time, ...
        stat_freq_feelpos_vs_feelneu_ZTR.freq, ...
        squeeze(stat_freq_feelpos_vs_feelneu_ZTR.stat));
axis xy;
xline(0, 'k--');
cb = colorbar; cb.Label.String = 't-value';
colormap(jet);
clim_val = max(abs(squeeze(stat_freq_feelpos_vs_feelneu_ZTR.stat)), [], 'all');
clim([-clim_val clim_val]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('FEEL Pos vs Neu (t-values)');

% Outline significant clusters (pos)
hold on;
if ~isempty(stat_freq_feelpos_vs_feelneu_ZTR.posclusters)
    sig_k = find([stat_freq_feelpos_vs_feelneu_ZTR.posclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelpos_vs_feelneu_ZTR.posclusterslabelmat == k);
        contour(stat_freq_feelpos_vs_feelneu_ZTR.time, ...
                stat_freq_feelpos_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
% Outline significant clusters (neg)
if ~isempty(stat_freq_feelpos_vs_feelneu_ZTR.negclusters)
    sig_k = find([stat_freq_feelpos_vs_feelneu_ZTR.negclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelpos_vs_feelneu_ZTR.negclusterslabelmat == k);
        contour(stat_freq_feelpos_vs_feelneu_ZTR.time, ...
                stat_freq_feelpos_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
hold off;

% --- Neg vs Neu ---
subplot(1,2,2)
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
title('FEEL Neg vs Neu (t-values)');

hold on;
if ~isempty(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
    sig_k = find([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelneg_vs_feelneu_ZTR.posclusterslabelmat == k);
        contour(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
                stat_freq_feelneg_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
if ~isempty(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)
    sig_k = find([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob] < sig_alpha);
    for k = sig_k
        mask = squeeze(stat_freq_feelneg_vs_feelneu_ZTR.negclusterslabelmat == k);
        contour(stat_freq_feelneg_vs_feelneu_ZTR.time, ...
                stat_freq_feelneg_vs_feelneu_ZTR.freq, ...
                double(mask), [0.5 0.5], 'k', 'LineWidth', 2);
    end
end
hold off;

%sgtitle(sprintf('%s: t-statistic TFR)', label_map(chan_name), sig_alpha), 'Interpreter', 'none');
sgtitle(sprintf('%s: t-statistic TFR', label_map(chan_name)), 'Interpreter', 'none');


fprintf('Positive clusters:\n'); disp(stat_freq_feelpos_vs_feelneu_ZTR.posclusters)
fprintf('Negative clusters:\n'); disp(stat_freq_feelpos_vs_feelneu_ZTR.negclusters)
disp([stat_freq_feelpos_vs_feelneu_ZTR.negclusters.prob])

fprintf('Positive clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.posclusters)
fprintf('Negative clusters:\n'); disp(stat_freq_feelneg_vs_feelneu_ZTR.negclusters)

disp([stat_freq_feelneg_vs_feelneu_ZTR.negclusters.prob])
disp([stat_freq_feelneg_vs_feelneu_ZTR.posclusters.prob])