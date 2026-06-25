function plot_lfp_with_gaps(t, x)
    % Plot a 1-channel signal and shade NaN regions in gray.
    figure('Name','Merged LFP with NaN gaps'); hold on;
    plot(t, x, 'b');
    nan_mask = isnan(x);
    d = diff([0 nan_mask 0]);
    s = find(d == 1);
    e = find(d == -1) - 1;
    yl = ylim;
    for i = 1:numel(s)
        patch([t(s(i)) t(e(i)) t(e(i)) t(s(i))], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.7 0.7 0.7], ...
              'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title('Merged LFP with NaN gaps'); xlabel('Time (s)'); ylabel('LFP');
end


