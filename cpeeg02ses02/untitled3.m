%% TFR - No Baseline (Raw Power)
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -2:0.05:4; % Updated time range
cfg.keeptrials = 'yes';
freq_raw = ft_freqanalysis(cfg, data_clean);

% Average across trials for plotting
cfg = [];
cfg.avgoverrpt = 'yes';
freq_avg_raw = ft_selectdata(cfg, freq_raw);

n_ch = length(freq_avg_raw.label);

% Plotting
figure('Name', 'Raw TFR (No Baseline)');
for ch = 1:n_ch
    subplot(ceil(n_ch/2), 2, ch);
    
    % Note: Raw power is always positive, so we don't use [-2 2]
    % imagesc will auto-scale, or you can set a specific range
    imagesc(freq_avg_raw.time, freq_avg_raw.freq, ...
        squeeze(freq_avg_raw.powspctrm(ch,:,:))); 
    
    axis xy;
    title(['Raw Power: ' freq_avg_raw.label{ch}]);
    
    xline(0, 'w--', 'LineWidth', 1.5); % Mark stimulus onset
    colorbar;
    colormap jet;
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');