%% TFR
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 1:1:40;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

%avg
cfg = [];
cfg.avgoverrpt = 'yes';

freq_avg = ft_selectdata(cfg, freq_clean_bl);

n_ch = length(freq_avg.label);

%plot
figure;
for ch = 1:n_ch
    
    subplot(ceil(n_ch/2), 2, ch);
    
    imagesc(freq_avg.time, freq_avg.freq, ...
        squeeze(freq_avg.powspctrm(ch,:,:)), [-2 2]); 
    
    axis xy;
    title(freq_avg.label{ch});
    
    xline(0,'w--');
    colorbar;

    colormap jet
end

xlabel('Time (s)');
ylabel('Frequency (Hz)');
%%
%% Power Spectral Density (PSD) for All Trials
% We use 'freq_clean' (pre-baseline) to see absolute power 
% or 'freq_clean_bl' to see change from baseline.

% 1. Average over the time dimension to get power per frequency per trial
cfg = [];
cfg.avgovertime = 'yes';
psd_data = ft_selectdata(cfg, freq_clean_bl); % Use freq_clean for raw power


% 2. Plotting

n_ch = 2;
for ch = 1:n_ch
    %subplot(ceil(n_ch/2), 2, ch);
    
    % Extract power: [trials x freq]
    % squeeze removes the channel and time singleton dimensions
    trial_pow = squeeze(psd_data.powspctrm(:, ch, :, :));
    
    % Calculate Mean and SEM
    mean_pow = mean(trial_pow, 1);
    sem_pow  = std(trial_pow, 0, 1) / sqrt(size(trial_pow, 1));
    x = psd_data.freq;
    upper = mean_pow + sem_pow;
    lower = mean_pow - sem_pow;

    
    % Plot shaded area
    fill([x, fliplr(x)], [upper, fliplr(lower)], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
    
    % Plot the Mean line
    plot(x, mean_pow, 'b', 'LineWidth', 1.5);
    
    title(['PSD: ' psd_data.label{ch}]);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Power');
    xticks([1 2 3 4 5 10 20 30 40])
end

n_ch = 3;
for ch = 1:n_ch
    % 1. Extract and Calculate
    trial_pow = squeeze(psd_data.powspctrm(:, ch, :, :));
    mean_pow = mean(trial_pow, 1);
    sem_pow  = std(trial_pow, 0, 1) / sqrt(size(trial_pow, 1));
    x = psd_data.freq;
    
    % 2. Create Figure
    figure; 
    hold on; % <--- THIS IS THE KEY: Apply to the current figure immediately
    
    % 3. Plot SEM first
    fill([x, fliplr(x)], [mean_pow + sem_pow, fliplr(mean_pow - sem_pow)], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
    
    % 4. Plot Mean on top
    plot(x, mean_pow, 'b', 'LineWidth', 1.5);
    
    % 5. Formatting
    title(['PSD: ' psd_data.label{ch}]);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Power Change (%)'); % Updated label for your percentage change
    xticks([1 2 3 4 5 10 20 30 40 50 60 70 80 90 100]); % Expanded for your foi
end


n_ch = 3;
figure; % Open ONE figure
hold on; % Keep EVERYTHING we draw on this figure

for ch = 1:n_ch
    trial_pow = squeeze(psd_data.powspctrm(:, ch, :, :));
    
    mean_pow = mean(trial_pow, 1);
    sem_pow  = std(trial_pow, 0, 1) / sqrt(size(trial_pow, 1));
    x = psd_data.freq;
    
    % Plot shaded area
    fill([x, fliplr(x)], [mean_pow + sem_pow, fliplr(mean_pow - sem_pow)], ...
         [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.2); % Lower Alpha helps when overlapping
    
    % Plot the Mean line
    plot(x, mean_pow, 'LineWidth', 1.5, 'DisplayName', psd_data.label{ch});
end

grid on;
xlabel('Frequency (Hz)');
ylabel('Power Change'); % Since you converted to %
legend('show'); % This will show which line is which channel
xticks([1 2 3 4 5 10 20 30 40]) 
%%
data = freq_avg.powspctrm;

neg_vals = data(data < 0);
length(neg_vals)
pos_vals = data(data > 0);
length(pos_vals)


%% TFR for single trials

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

%plotting
n_trials = size(freq_clean_bl.powspctrm, 1);
block_size = 23;
n_blocks = floor(n_trials / block_size);
n_ch = 6;
for b = 1:n_blocks
    
    idx_start = (b-1)*block_size + 1;
    idx_end   = b*block_size;

    for ch = 1:n_ch
        
        figure;
        
        for t = 1:block_size
            
            trial_idx = idx_start + t - 1;
            
            subplot(5,5,t); % 5x5 grid
            
            data_trial = squeeze(freq_clean_bl.powspctrm(trial_idx,ch,:,:));
            
            imagesc(freq_clean_bl.time, freq_clean_bl.freq, data_trial, [-2 2]);
            
            axis xy;
            title(['T' num2str(trial_idx)]);
            
            xline(0,'w--');
            colormap jet
        end
        
        sgtitle(['Block ' num2str(b) ' - ' freq_clean_bl.label{ch}]);
    end
end

%% PSD for all trails

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

%plotting
n_trials = size(freq_clean_bl.powspctrm, 1);
block_size = 23;
n_blocks = floor(n_trials / block_size);
n_ch = 6;

for ch = 1:n_ch
        
    figure;
   
        
    trial_idx = idx_start + t - 1;
    
    subplot(5,5,t); % 5x5 grid
    
    data_trial = squeeze(freq_clean_bl.powspctrm(trial_idx,ch,:,:));
    
    imagesc(freq_clean_bl.time, freq_clean_bl.freq, data_trial, [-2 2]);
    
    axis xy;
    title(['T' num2str(trial_idx)]);
    
    xline(0,'w--');
    colormap jet
end
 

%% TFR for single trials

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);
  % create ONE figure window

t_idx = freq_clean_bl.time >= 0 & freq_clean_bl.time <= 2;
f = freq_clean_bl.freq;


for ch = 1:n_ch
    % FEEL + POS
    idx = feel_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    psd_feel_pos = squeeze(mean(mean(data,4),1));
    
    % FEEL + NEU
    idx = feel_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    psd_feel_neu = squeeze(mean(mean(data,4),1));
    
    % FEEL + NEG
    idx = feel_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    psd_feel_neg = squeeze(mean(mean(data,4),1));
    
    % TONE + POS
    idx = tone_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    psd_tone_pos = squeeze(mean(mean(data,4),1));
    
    % TONE + NEU
    idx = tone_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    psd_tone_neu = squeeze(mean(mean(data,4),1));
    
    % TONE + NEG
    idx = tone_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    psd_tone_neg = squeeze(mean(mean(data,4),1));


    
    figure(); hold on;
    
    plot(f, psd_feel_pos, 'r',  'LineWidth', 2);
    plot(f, psd_feel_neu, 'g',  'LineWidth', 2);
    plot(f, psd_feel_neg, 'b',  'LineWidth', 2);
    
    plot(f, psd_tone_pos, '--r', 'LineWidth', 2);
    plot(f, psd_tone_neu, '--g', 'LineWidth', 2);
    plot(f, psd_tone_neg, '--b', 'LineWidth', 2);
    
    legend({'FEEL-POS','FEEL-NEU','FEEL-NEG', ...
            'TONE-POS','TONE-NEU','TONE-NEG'});
    
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['PSD (0–2s) - Channel ' freq_clean_bl.label{ch}]);
    grid on;

end



%% psd for single trials with SEM

cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);
  % create ONE figure window

t_idx = freq_clean_bl.time >= 0 & freq_clean_bl.time <= 2;
f = freq_clean_bl.freq;
n_ch = 3;

for ch = 1:n_ch
    % FEEL + POS
    idx = feel_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));   % average over time (per trial) [trials x freq]
    psd_feel_pos = mean(data_tavg, 1);   % mean across trials
    sem_feel_pos = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1)); % SEM across trials
    
    % FEEL + NEU
    idx = feel_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_feel_neu  = mean(data_tavg, 1);
    sem_feel_neu = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % FEEL + NEG
    idx = feel_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_feel_neg  = mean(data_tavg, 1);
    sem_feel_neg = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % TONE + POS
    idx = tone_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_pos  = mean(data_tavg, 1);
    sem_tone_pos = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
   
    % TONE + NEU
    idx = tone_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_neu = mean(data_tavg, 1);
    sem_tone_neu = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % TONE + NEG
    idx = tone_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_neg = mean(data_tavg, 1);
    sem_tone_neg = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));

    figure(); hold on;
    
    % shaded
    fill([f fliplr(f)], ...
     [psd_feel_pos+sem_feel_pos fliplr(psd_feel_pos-sem_feel_pos)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, psd_feel_pos, 'r',  'LineWidth', 2);


    fill([f fliplr(f)], ...
     [psd_feel_neu+sem_feel_neu fliplr(psd_feel_neu-sem_feel_neu)], ...
     'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, psd_feel_neu, 'g',  'LineWidth', 2);

    fill([f fliplr(f)], ...
     [psd_feel_neg+sem_feel_neg fliplr(psd_feel_neg-sem_feel_neg)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, psd_feel_neg, 'b',  'LineWidth', 2);
    
    fill([f fliplr(f)], ...
     [psd_tone_pos+sem_tone_pos fliplr(psd_tone_pos-sem_tone_pos)], ...
     'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, psd_tone_pos, '--r', 'LineWidth', 2);

    fill([f fliplr(f)], ...
     [psd_tone_neu+sem_tone_neu fliplr(psd_tone_neu-sem_tone_neu)], ...
     'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, psd_tone_neu, '--g', 'LineWidth', 2);

    fill([f fliplr(f)], ...
     [psd_tone_neg+sem_tone_neg fliplr(psd_tone_neg-sem_tone_neg)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, psd_tone_neg, '--b', 'LineWidth', 2);
    
    legend({'FEEL-POS','FEEL-NEU','FEEL-NEG', ...
            'TONE-POS','TONE-NEU','TONE-NEG'});
    
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['PSD (0–2s) - Channel ' freq_clean_bl.label{ch}]);
    grid on;

end

%% psd for single trials with SEM proper colors


cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 2:1:100;
cfg.toi = -0.5:0.05:3;

cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);
  % create ONE figure window

t_idx = freq_clean_bl.time >= 0 & freq_clean_bl.time <= 1;
f = freq_clean_bl.freq;
n_ch = 6;

for ch = 1:n_ch
    % FEEL + POS
    idx = feel_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));   % average over time (per trial) [trials x freq]
    psd_feel_pos = mean(data_tavg, 1);   % mean across trials
    sem_feel_pos = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1)); % SEM across trials
    
    % FEEL + NEU
    idx = feel_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_feel_neu  = mean(data_tavg, 1);
    sem_feel_neu = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % FEEL + NEG
    idx = feel_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_feel_neg  = mean(data_tavg, 1);
    sem_feel_neg = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % TONE + POS
    idx = tone_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_pos  = mean(data_tavg, 1);
    sem_tone_pos = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
   
    % TONE + NEU
    idx = tone_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_neu = mean(data_tavg, 1);
    sem_tone_neu = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % TONE + NEG
    idx = tone_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_neg = mean(data_tavg, 1);
    sem_tone_neg = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));

    figure(); hold on;
    h = gobjects(3,1); % store line handles

    % define 6 distinct colors
    c1 = [0.90 0.10 0.10]; % red
    c2 = [0.10 0.70 0.10]; % green
    c3 = [0.10 0.30 0.90]; % blue
    c4 = [0.90 0.60 0.10]; % orange
    c5 = [0.60 0.20 0.80]; % purple
    c6 = [0.10 0.70 0.70]; % cyan
    
    % FEEL-POS
    fill([f fliplr(f)], ...
         [psd_feel_pos+sem_feel_pos fliplr(psd_feel_pos-sem_feel_pos)], ...
         c1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
   h(1) = plot(f, psd_feel_pos, 'Color', c1, 'LineWidth', 2);
    
    % FEEL-NEU
    fill([f fliplr(f)], ...
         [psd_feel_neu+sem_feel_neu fliplr(psd_feel_neu-sem_feel_neu)], ...
         c2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h(2) = plot(f, psd_feel_neu, 'Color', c2, 'LineWidth', 2);
    
    % FEEL-NEG
    fill([f fliplr(f)], ...
         [psd_feel_neg+sem_feel_neg fliplr(psd_feel_neg-sem_feel_neg)], ...
         c3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h(3) = plot(f, psd_feel_neg, 'Color', c3, 'LineWidth', 2);
    
    % % TONE-POS
    % fill([f fliplr(f)], ...
    %      [psd_tone_pos+sem_tone_pos fliplr(psd_tone_pos-sem_tone_pos)], ...
    %      c4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % h(4) = plot(f, psd_tone_pos, '--', 'Color', c4, 'LineWidth', 2);
    % 
    % % TONE-NEU
    % fill([f fliplr(f)], ...
    %      [psd_tone_neu+sem_tone_neu fliplr(psd_tone_neu-sem_tone_neu)], ...
    %      c5, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % h(5) = plot(f, psd_tone_neu, '--', 'Color', c5, 'LineWidth', 2);
    % 
    % % TONE-NEG
    % fill([f fliplr(f)], ...
    %      [psd_tone_neg+sem_tone_neg fliplr(psd_tone_neg-sem_tone_neg)], ...
    %      c6, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    %  h(6) =plot(f, psd_tone_neg, '--', 'Color', c6, 'LineWidth', 2);
    % 

    legend(h, {'FEEL-POS','FEEL-NEU','FEEL-NEG'});%, ...
        %   'TONE-POS','TONE-NEU','TONE-NEG'});
    
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['PSD (0–1s) - Channel ' freq_clean_bl.label{ch}]);
    axis tight
    xticks()
    grid on;

end





for ch = 1:n_ch
    % FEEL + POS
    idx = feel_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));   % average over time (per trial) [trials x freq]
    psd_feel_pos = mean(data_tavg, 1);   % mean across trials
    sem_feel_pos = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1)); % SEM across trials
    
    % FEEL + NEU
    idx = feel_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_feel_neu  = mean(data_tavg, 1);
    sem_feel_neu = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % FEEL + NEG
    idx = feel_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_feel_neg  = mean(data_tavg, 1);
    sem_feel_neg = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % TONE + POS
    idx = tone_idx & pos_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_pos  = mean(data_tavg, 1);
    sem_tone_pos = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
   
    % TONE + NEU
    idx = tone_idx & neu_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_neu = mean(data_tavg, 1);
    sem_tone_neu = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));
    
    % TONE + NEG
    idx = tone_idx & neg_idx;
    data = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
    data_tavg = squeeze(mean(data, 4));
    psd_tone_neg = mean(data_tavg, 1);
    sem_tone_neg = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1));

    figure(); hold on;
    h = gobjects(3,1); % store line handles

    % define 6 distinct colors
    c1 = [0.90 0.10 0.10]; % red
    c2 = [0.10 0.70 0.10]; % green
    c3 = [0.10 0.30 0.90]; % blue
    c4 = [0.90 0.60 0.10]; % orange
    c5 = [0.60 0.20 0.80]; % purple
    c6 = [0.10 0.70 0.70]; % cyan
    
   %  % FEEL-POS
   %  fill([f fliplr(f)], ...
   %       [psd_feel_pos+sem_feel_pos fliplr(psd_feel_pos-sem_feel_pos)], ...
   %       c1, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
   % h(1) = plot(f, psd_feel_pos, 'Color', c1, 'LineWidth', 2);
   % 
   %  % FEEL-NEU
   %  fill([f fliplr(f)], ...
   %       [psd_feel_neu+sem_feel_neu fliplr(psd_feel_neu-sem_feel_neu)], ...
   %       c2, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
   %  h(2) = plot(f, psd_feel_neu, 'Color', c2, 'LineWidth', 2);
   % 
   %  % FEEL-NEG
   %  fill([f fliplr(f)], ...
   %       [psd_feel_neg+sem_feel_neg fliplr(psd_feel_neg-sem_feel_neg)], ...
   %       c3, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
   %  h(3) = plot(f, psd_feel_neg, 'Color', c3, 'LineWidth', 2);
   % 
    % TONE-POS
    fill([f fliplr(f)], ...
         [psd_tone_pos+sem_tone_pos fliplr(psd_tone_pos-sem_tone_pos)], ...
         c4, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h(1) = plot(f, psd_tone_pos, '--', 'Color', c4, 'LineWidth', 2);

    % TONE-NEU
    fill([f fliplr(f)], ...
         [psd_tone_neu+sem_tone_neu fliplr(psd_tone_neu-sem_tone_neu)], ...
         c5, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h(2) = plot(f, psd_tone_neu, '--', 'Color', c5, 'LineWidth', 2);

    % TONE-NEG
    fill([f fliplr(f)], ...
         [psd_tone_neg+sem_tone_neg fliplr(psd_tone_neg-sem_tone_neg)], ...
         c6, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
     h(3) =plot(f, psd_tone_neg, '--', 'Color', c6, 'LineWidth', 2);


    legend(h, {'TONE-POS','TONE-NEU','TONE-NEG'});
    

    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(['PSD (0–1s) - Channel ' freq_clean_bl.label{ch}]);
    axis tight
    %xticks()
    grid on;

end
%% PSD of all trials
cfg = [];
cfg.method = 'wavelet';
cfg.output = 'pow';
cfg.width  = 6;
cfg.foi = 1:1:40;
cfg.toi = -0.5:0.05:3;
cfg.keeptrials = 'yes';

freq_clean = ft_freqanalysis(cfg, data_clean);

%baseline
cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'relchange';   % or 'db'
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

t_idx = freq_clean_bl.time >= 0 & freq_clean_bl.time <= 3;
f = freq_clean_bl.freq;
ch = 1;

data = freq_clean_bl.powspctrm(:, ch, :, t_idx);
% psd_all1= squeeze(mean(mean(data,4),1));

data_tavg = squeeze(mean(data, 4));   % average over time (per trial) [trials x freq]
psd_all = mean(data_tavg, 1);   % mean across trials
sem_all = std(data_tavg, 0, 1) / sqrt(size(data_tavg,1)); 


% --- Plotting ---
figure('Color', 'w');
hold on;

% 1. Define the Upper and Lower bounds for the shade
upper_bound = psd_all + sem_all;
lower_bound = psd_all - sem_all;

% 2. Create the filled area (Shaded SEM)
% We concatenate the x-axis (f) and the bounds to create a closed loop
x_patch = [f, fliplr(f)];
y_patch = [upper_bound, fliplr(lower_bound)];

fill(x_patch, y_patch, [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

% 3. Plot the Mean Line on top
plot(f, psd_all, 'Color', [0 0.4470 0.7410], 'LineWidth', 2);

% --- Formatting ---
grid on;
xlabel('Frequency (Hz)');
ylabel('Power (uV^2/Hz)');
title(['PSD with SEM Shading - Channel ', num2str(ch)]);
set(gca, 'FontSize', 12);