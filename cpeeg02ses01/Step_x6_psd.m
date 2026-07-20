%% Figure: means (top) + raw effects (bottom), ZERO_TWO_RIGHT

figure;

% --- Top: all 3 means with SEM ---
subplot(2,1,1); hold on;

fill([freq_axis, fliplr(freq_axis)], [mean_pos+sem_pos, fliplr(mean_pos-sem_pos)], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([freq_axis, fliplr(freq_axis)], [mean_neu+sem_neu, fliplr(mean_neu-sem_neu)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([freq_axis, fliplr(freq_axis)], [mean_neg+sem_neg, fliplr(mean_neg-sem_neg)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

plot(freq_axis, mean_pos, 'g', 'LineWidth', 2);
plot(freq_axis, mean_neu, 'b', 'LineWidth', 2);
plot(freq_axis, mean_neg, 'r', 'LineWidth', 2);

xlabel('Frequency (Hz)');
ylabel('Power (relchange)');
title(sprintf('%s: FEEL Pos vs Neu vs Neg', label_map(chan_name)), 'Interpreter', 'none');
legend('', '', '', 'Pos', 'Neu', 'Neg', 'Location', 'best');
xlim([1 40]);
hold off;

% --- Bottom: raw effects ---
subplot(2,1,2); hold on;

plot(freq_axis, raweffect_neg_neu, 'r', 'LineWidth', 2);
plot(freq_axis, raweffect_pos_neu, 'g', 'LineWidth', 2);
yline(0, 'k--');

xlabel('Frequency (Hz)');
ylabel('Power diff');
title('Raw effects: Neg - Neu and Pos - Neu');
legend('Neg - Neu', 'Pos - Neu', 'Location', 'best');
xlim([1 40]);
hold off;