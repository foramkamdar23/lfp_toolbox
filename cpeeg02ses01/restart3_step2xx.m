psd_twin = [0 3];     % s, time window to average over
f = freq_clean_bl.freq;
t_idx = freq_clean_bl.time >= psd_twin(1) & freq_clean_bl.time <= psd_twin(2);

figure;

% --- Subplot 1: ONE_THREE_LEFT ---
lab = 'ONE_THREE_LEFT';
ch  = find(strcmp(freq_clean_bl.label, lab), 1);
d      = squeeze(freq_clean_bl.powspctrm(:, ch, :, t_idx));
d_tavg = mean(d, 3);
psd_mean    = mean(d_tavg, 1);
psd_sem     = std(d_tavg, 0, 1) / sqrt(size(d_tavg, 1));
upper_bound = psd_mean + psd_sem;
lower_bound = psd_mean - psd_sem;
subplot(2, 2, 1); hold on;
fill([f(:)' fliplr(f(:)')], [upper_bound(:)' fliplr(lower_bound(:)')], ...
    cond_colors_light{k}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(f, psd_mean, 'b', 'LineWidth', 2);
hold off; grid on;
xlabel('Frequency (Hz)'); ylabel('Power (relchange)');
title(lab, 'Interpreter', 'none');
set(gca, 'FontSize', 12);
