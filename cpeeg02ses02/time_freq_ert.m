%% Time frequency power 

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';

cfg.foi = 2:1:100;

cfg.t_ftimwin = 0.5 * ones(size(cfg.foi));   % 500 ms window
cfg.toi = -0.5:0.05:3.5;                    % sliding windows

cfg.keeptrials = 'yes';

freq = ft_freqanalysis(cfg, data_bl);

cfg = [];
cfg.avgoverrpt = 'yes';
freq_avg = ft_selectdata(cfg, freq);

ch = 1;

figure;
imagesc(freq_avg.time, freq_avg.freq, squeeze(freq_avg.powspctrm(ch,:,:)));
axis xy;

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['MTMCONVOL (pwelch-like) | ' freq_avg.label{ch}]);

colorbar;
xline(0,'w--','LineWidth',2);

% perception and regulation
%
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.taper = 'hanning';
cfg.foi = 2:2:100;
cfg.t_ftimwin = 0.5 * ones(size(cfg.foi));
cfg.toi = -0.5:0.05:3.5;
cfg.keeptrials = 'yes';

freq_feel = ft_freqanalysis(cfg, data_feel);
freq_tone = ft_freqanalysis(cfg, data_tone);

cfg = [];
cfg.avgoverrpt = 'yes';

freq_feel_avg = ft_selectdata(cfg, freq_feel);
freq_tone_avg = ft_selectdata(cfg, freq_tone);

diff_pow = freq_feel_avg.powspctrm - freq_tone_avg.powspctrm;

ch = 1;

figure;
imagesc(freq_feel_avg.time, freq_feel_avg.freq, squeeze(diff_pow(ch,:,:)));
axis xy;

xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['FEEL - TONE | ' freq_feel_avg.label{ch}]);

colorbar;
xline(0,'w--','LineWidth',2);