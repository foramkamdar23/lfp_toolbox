function plot_artifact_epochs(artifact_env, artifact_thresh, trl, clean_start, clean_end, subject, block, Fs)
% PLOT_ARTIFACT_EPOCHS  Diagnostic plot of artifact envelope with epoch patches
%
% INPUTS:
%   artifact_env    - 1 x nSamples envelope signal
%   artifact_thresh - scalar threshold value
%   trl             - nTrials x 4 trial matrix [beg end offset id]
%   clean_start     - first sample of clean window
%   clean_end       - last sample of clean window
%   subject         - string, subject ID
%   block           - string, block ID
%   Fs              - sampling rate (for x-axis in seconds)

figure('Color','w', 'Position',[100 100 1400 500]);
hold on;

% --- Shade clean window in light yellow ---
y_limits = [0, max(artifact_env) * 1.05];
patch([clean_start clean_start clean_end clean_end], ...
      [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
      [1 1 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
      'DisplayName', 'Clean window');

% --- Shade each epoch as green patch ---
for tr = 1:size(trl, 1)
    beg = trl(tr, 1);
    en  = trl(tr, 2);
    patch([beg beg en en], ...
          [y_limits(1) y_limits(2) y_limits(2) y_limits(1)], ...
          'g', 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
          'HandleVisibility', 'off');
end

% --- Plot envelope ---
plot(artifact_env, 'b', 'LineWidth', 1, 'DisplayName', 'Artifact env');

% --- Threshold line ---
yline(artifact_thresh, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');

% --- Clean window boundaries ---
xline(clean_start, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Clean start');
xline(clean_end,   'k-', 'LineWidth', 1.5, 'DisplayName', 'Clean end');

% Add a single dummy patch for the legend
patch(nan, nan, 'g', 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'Epochs');

ylim(y_limits);
xlabel('Sample index');
ylabel('Envelope amplitude');
title(sprintf('Artifact envelope with epochs — %s %s', subject, block));
legend('Location', 'northeast', 'FontSize', 11);

end