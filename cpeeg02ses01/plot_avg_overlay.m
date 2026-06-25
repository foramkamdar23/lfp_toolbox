function plot_avg_overlay(avgs, colors, names, layout, ttl)
    % Overlay timelock averages on a LEFT/RIGHT grid laid out by `layout`.
    [nr, nc] = size(layout);
    figure('Name', ttl);
    for r = 1:nr
        for c = 1:nc
            lab = layout{r,c};
            if isempty(lab), continue; end
            subplot(nr, nc, (r-1)*nc + c); hold on;
            for k = 1:numel(avgs)
                ch = find(strcmp(avgs{k}.label, lab), 1);
                if isempty(ch), continue; end
                plot(avgs{k}.time, avgs{k}.avg(ch,:), 'Color', colors{k}, 'LineWidth', 2);
            end
            xline(0, 'k--');
            title(lab, 'Interpreter', 'none'); xlabel('Time (s)'); ylabel('LFP');
        end
    end
    legend(names);
    sgtitle(ttl);
end
