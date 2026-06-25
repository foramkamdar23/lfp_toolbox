%% TFR
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 1:1:100;
%cfg.toi = -0.5:0.05:3;
cfg.toi = -1:0.05:3.5; % Updated time range
cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%avg
cfg = [];
cfg.avgoverrpt = 'yes';

freq_avg = ft_selectdata(cfg, freq_clean);

n_ch = length(freq_avg.label);

%plots
figure;

for ch = 1:n_ch
    
    subplot(ceil(n_ch/2), 2, ch);
    
    imagesc(freq_avg.time, freq_avg.freq, ...
        squeeze(freq_avg.powspctrm(ch,:,:)), [1 400]); 
    
    axis xy;
    title(freq_avg.label{ch});
    
    xline(0,'w--');
    colorbar;

    colormap jet
end

xlabel('Time (s)');
ylabel('Frequency (Hz)');
sgtitle('Raw TFR')