function plot_lfp_trials(data, layout, ttl)
    % Per-channel grid (LEFT col = left hemisphere, RIGHT col = right):
    % all trials in gray + black mean. layout is an {nRows x 2} cell of labels.
    [nr, nc] = size(layout);
    figure('Name', ttl);
    for r = 1:nr
        for c = 1:nc
            lab = layout{r,c};
            if isempty(lab), continue; end
            ch = find(strcmp(data.label, lab), 1);
            if isempty(ch), continue; end
            subplot(nr, nc, (r-1)*nc + c); hold on;
            n_tr = numel(data.trial);
            M = zeros(n_tr, numel(data.time{1}));
            for i = 1:n_tr
                M(i,:) = data.trial{i}(ch,:);
            end
            plot(data.time{1}, M', 'Color', [0.8 0.8 0.8]);
            plot(data.time{1}, mean(M, 1, 'omitnan'), 'k', 'LineWidth', 2);
            xline(0, 'r--');
            title(lab, 'Interpreter', 'none'); xlabel('Time (s)'); ylabel('LFP');
        end
    end
    sgtitle(ttl);
end
