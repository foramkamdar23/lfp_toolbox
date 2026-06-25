function plot_tfr_grid(freq, layout, clim, ttl, label_map)
    % Per-channel TFR grid (LEFT col = left hemisphere, RIGHT col = right).
    % freq.powspctrm is [channels x freq x time] after averaging over trials.
    % label_map: containers.Map from channel label -> descriptive title (optional)
    [nr, nc] = size(layout);
    figure('Name', ttl);
    for r = 1:nr
        for c = 1:nc
            lab = layout{r,c};
            ch = find(strcmp(freq.label, lab), 1);
            subplot(nr, nc, (r-1)*nc + c);
            imagesc(freq.time, freq.freq, squeeze(freq.powspctrm(ch,:,:)), clim);
            axis xy;
            xline(0, 'k--');                 % image onset

            if nargin >= 5 && isKey(label_map, lab)
                plot_title = label_map(lab);
            else
                plot_title = lab;
            end
            title(plot_title, 'Interpreter', 'none');

            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            colorbar;
        end
    end
    colormap jet
    sgtitle(ttl);
end