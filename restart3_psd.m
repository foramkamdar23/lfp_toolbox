% Plot layout: left column = LEFT hemisphere, right column = RIGHT.
clear all;
clc;

load("C:\Users\fkamdar\Desktop\repos\lfp_toolbox\derivatives\sub-CPEEG01_ses-02_data_clean.mat")

% Plot layout: left column = LEFT hemisphere, right column = RIGHT.
plot_layout = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT' ...
                 'ZERO_TWO_LEFT', 'ZERO_TWO_RIGHT'};
 
fprintf('Loaded %d clean trials, channels: %s\n', ...
        numel(data_clean.trial), strjoin(data_clean.label, ', '));

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
tfr_baseline = [-0.5 -0.2];     % s, pre-image baseline window
tfr_basetype = 'relchange';     % 'relchange' (fractional) or 'db' % change

cfg = [];
cfg.baseline = tfr_baseline;
cfg.baselinetype = tfr_basetype; 
freq_clean_bl = ft_freqbaseline(cfg, freq_clean);

%avg
cfg = [];
cfg.avgoverrpt = 'yes';
freq_avg = ft_selectdata(cfg, freq_clean_bl);

% PLOT TFR
tfr_clim = [-2 2];
plot_tfr_grid(freq_avg, plot_layout, tfr_clim, ...
              sprintf('TFR (%s baseline %g to %g s)', tfr_basetype, tfr_baseline));
 




%% PSD for all trials of 1 channel:
% cfg = [];
% cfg.method = 'wavelet';
% cfg.output = 'pow';
% cfg.width  = 6;
% cfg.foi = 1:1:40;
% cfg.toi = -0.5:0.05:3;
% cfg.keeptrials = 'yes';
% 
% freq_clean = ft_freqanalysis(cfg, data_clean);
% 
% %baseline
% tfr_baseline = [-0.5 -0.2];     % s, pre-image baseline window
% tfr_basetype = 'relchange';     % 'relchange' (fractional) or 'db' % change
% 
% cfg = [];
% cfg.baseline = tfr_baseline;
% cfg.baselinetype = tfr_basetype; 
% freq_clean_bl = ft_freqbaseline(cfg, freq_clean);
%%
%% PSD (time-averaged, with SEM)
psd_twin = [0 3];     % s, time window to average over
f = freq_clean_bl.freq;
t_idx = freq_clean_bl.time >= psd_twin(1) & freq_clean_bl.time <= psd_twin(2);

plot_layout = {'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT' };
[nr, nc] = size(plot_layout);
cond_colors_light = {[0.7 0.8 1.0], [1.0 0.7 0.7]}; 
figure;
for r = 1:nr
    for c = 1:nc
        lab = plot_layout{r,c};
        ch  = find(strcmp(freq_clean_bl.label, lab), 1);

        % [trials x freq]: average over time window per trial
        d      = squeeze(freq_clean_bl.powspctrm(:, ch, :, t_idx));  % [trials x freq x t_win]
        d_tavg = mean(d, 3);                                          % [trials x freq]

        psd_mean    = mean(d_tavg, 1);
        psd_sem     = std(d_tavg, 0, 1) / sqrt(size(d_tavg, 1));
        upper_bound = psd_mean + psd_sem;
        lower_bound = psd_mean - psd_sem;

        subplot(nr, nc, (r-1)*nc + c); hold on;
        
        fill([f(:)' fliplr(f(:)')], ...
            [upper_bound(:)' fliplr(lower_bound(:)')], ...
            cond_colors_light{c}, 'EdgeColor', 'none', 'FaceAlpha', 0.3 );  

        % plot(f, upper_bound, 'Color', [0.9 0.9 0.9], 'LineWidth', 1, 'LineStyle', '-');
        % plot(f, lower_bound, 'Color', [0.9 0.9 0.9], 'LineWidth', 1, 'LineStyle', '-');
        plot(f, psd_mean, 'b', 'LineWidth', 2);
        hold off
        grid on;
        xlabel('Frequency (Hz)'); ylabel('Power (relchange)');
        title(lab, 'Interpreter', 'none');
        set(gca, 'FontSize', 12);
    end
end
sgtitle(sprintf('PSD %g-%g s, mean +/- SEM', psd_twin));


%% PSD per categories 

pos_idx = strcmp(trial_info_clean.Valence, 'Pos');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');


% cfg = [];
% cfg.trials = find(feel_idx);
% data_feel = ft_selectdata(cfg, data_clean);

%% PSD per condition (FEEL vs TONE)
psd_twin = [0 3];
f     = freq_clean_bl.freq;
t_idx = freq_clean_bl.time >= psd_twin(1) & freq_clean_bl.time <= psd_twin(2);

feel_idx = strcmp(trial_info_clean.Condition, 'FEEL');
tone_idx = strcmp(trial_info_clean.Condition, 'TONE');

% plot_layout = { 'ONE_THREE_LEFT', 'ONE_THREE_RIGHT'; ...
%                 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT' };

plot_layout = {'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT' };

conditions  = {feel_idx, tone_idx};
cond_names  = {'FEEL', 'TONE'};
cond_colors = {[0.2 0.4 1.0], [1.0 0.2 0.2]};       % blue, red
cond_colors_light = {[0.7 0.8 1.0], [1.0 0.7 0.7]};  % light blue, light red

[nr, nc] = size(plot_layout);
figure;
for r = 1:nr
    for c = 1:nc
        lab = plot_layout{r,c};
        ch  = find(strcmp(freq_clean_bl.label, lab), 1);

        subplot(nr, nc, (r-1)*nc + c); hold on;

        h = gobjects(2,1);   % store main line handles for legend

        for k = 1:2
            idx = conditions{k};

            % d      = squeeze(freq_clean_bl.powspctrm(idx, ch, :, t_idx));  % [trials x freq x t_win]
            % d_tavg = mean(d, 3);                                            % [trials x freq]
            d = freq_clean_bl.powspctrm(idx, ch, :, t_idx);   % [trials x 1 x freq x t_win]
            d_tavg = mean(d, 4);                               % average over time -> [trials x 1 x freq]
            d_tavg = reshape(d_tavg, size(d_tavg,1), size(d_tavg,3));  % [trials x freq]

            psd_mean    = mean(d_tavg, 1);
            psd_sem     = std(d_tavg, 0, 1) / sqrt(size(d_tavg, 1));
            upper_bound = psd_mean + psd_sem;
            lower_bound = psd_mean - psd_sem;
            
            fill([f(:)' fliplr(f(:)')], ...
                [upper_bound(:)' fliplr(lower_bound(:)')], ...
                cond_colors_light{k}, 'EdgeColor', 'none', 'FaceAlpha', 0.3 );
            %plot(f, upper_bound, 'Color', cond_colors_light{k}, 'LineWidth', 1);
            %plot(f, lower_bound, 'Color', cond_colors_light{k}, 'LineWidth', 1);
            h(k) = plot(f(:)', psd_mean(:)', 'Color', cond_colors{k}, 'LineWidth', 2);
        end
        
        grid on;
        xlabel('Frequency (Hz)'); ylabel('Power (relchange)');
        title(lab, 'Interpreter', 'none');
        set(gca, 'FontSize', 12);
    end
end

legend(h, cond_names, 'Location', 'best');
sgtitle(sprintf('PSD %g-%g s: FEEL vs TONE', psd_twin));


%% PSD per valence (Neg vs Neu vs Pos)
psd_twin = [0 3];
f     = freq_clean_bl.freq;
t_idx = freq_clean_bl.time >= psd_twin(1) & freq_clean_bl.time <= psd_twin(2);

neg_idx = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx = strcmp(trial_info_clean.Valence, 'Pos');

plot_layout = {'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT' };


conditions  = {neg_idx, neu_idx, pos_idx};
cond_names  = {'Neg', 'Neu', 'Pos'};
cond_colors = {[0.8 0.1 0.1], [0.4 0.4 0.4], [0.1 0.6 0.2]};            % red, gray, green
cond_colors_light = {[1.0 0.7 0.7], [0.75 0.75 0.75], [0.7 0.9 0.75]};  % light red, light gray, light green

[nr, nc] = size(plot_layout);
figure;
for r = 1:nr
    for c = 1:nc
        lab = plot_layout{r,c};
        ch  = find(strcmp(freq_clean_bl.label, lab), 1);
        subplot(nr, nc, (r-1)*nc + c); hold on;
        h = gobjects(3,1);

        for k = 1:3
            idx = conditions{k};

            d      = freq_clean_bl.powspctrm(idx, ch, :, t_idx);   % [trials x 1 x freq x t_win]
            d_tavg = mean(d, 4);                                     % [trials x 1 x freq]
            d_tavg = reshape(d_tavg, size(d_tavg,1), size(d_tavg,3)); % [trials x freq]

            psd_mean    = mean(d_tavg, 1);
            psd_sem     = std(d_tavg, 0, 1) / sqrt(size(d_tavg, 1));
            upper_bound = psd_mean + psd_sem;
            lower_bound = psd_mean - psd_sem;

            fill([f(:)' fliplr(f(:)')], ...
                 [upper_bound(:)' fliplr(lower_bound(:)')], ...
                 cond_colors_light{k}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

            h(k) = plot(f(:)', psd_mean(:)', 'Color', cond_colors{k}, 'LineWidth', 2);
        end

        grid on;
        xlabel('Frequency (Hz)'); ylabel('Power (relchange)');
        title(lab, 'Interpreter', 'none');
        set(gca, 'FontSize', 12);
    end
end
legend(h, cond_names, 'Location', 'best');
sgtitle(sprintf('PSD %g-%g s: Neg vs Neu vs Pos', psd_twin));

%% PSD per condition x valence (FEEL and TONE, split by Neg/Neu/Pos)
psd_twin = [0 3];
f     = freq_clean_bl.freq;
t_idx = freq_clean_bl.time >= psd_twin(1) & freq_clean_bl.time <= psd_twin(2);

feel_idx = strcmp(trial_info_clean.Condition, 'FEEL');
tone_idx = strcmp(trial_info_clean.Condition, 'TONE');
neg_idx  = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx  = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx  = strcmp(trial_info_clean.Valence, 'Pos');

plot_layout = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT' };

cond_names  = {'Neg', 'Neu', 'Pos'};
cond_colors       = {[0.8 0.1 0.1], [0.4 0.4 0.4], [0.1 0.6 0.2]};
cond_colors_light = {[1.0 0.7 0.7], [0.75 0.75 0.75], [0.7 0.9 0.75]};

main_conditions = {feel_idx, tone_idx};
main_names      = {'FEEL', 'TONE'};

[nr, nc] = size(plot_layout);

for m = 1:2
    main_idx = main_conditions{m};

    conditions = { main_idx & neg_idx, ...
                   main_idx & neu_idx, ...
                   main_idx & pos_idx };

    figure;
    for r = 1:nr
        for c = 1:nc
            lab = plot_layout{r,c};
            ch  = find(strcmp(freq_clean_bl.label, lab), 1);
            subplot(nr, nc, (r-1)*nc + c); hold on;
            h = gobjects(3,1);

            for k = 1:3
                idx = conditions{k};

                d      = freq_clean_bl.powspctrm(idx, ch, :, t_idx);    % [trials x 1 x freq x t_win]
                d_tavg = mean(d, 4);                                      % [trials x 1 x freq]
                d_tavg = reshape(d_tavg, size(d_tavg,1), size(d_tavg,3)); % [trials x freq]

                psd_mean    = mean(d_tavg, 1);
                psd_sem     = std(d_tavg, 0, 1) / sqrt(size(d_tavg, 1));
                upper_bound = psd_mean + psd_sem;
                lower_bound = psd_mean - psd_sem;

                fill([f(:)' fliplr(f(:)')], ...
                     [upper_bound(:)' fliplr(lower_bound(:)')], ...
                     cond_colors_light{k}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

                h(k) = plot(f(:)', psd_mean(:)', 'Color', cond_colors{k}, 'LineWidth', 2);
            end

            grid on;
            xlabel('Frequency (Hz)'); ylabel('Power (relchange)');
            title(lab, 'Interpreter', 'none');
            set(gca, 'FontSize', 12);
        end
    end
    legend(h, cond_names, 'Location', 'best');
    sgtitle(sprintf('PSD %g-%g s: %s — Neg vs Neu vs Pos', psd_twin(1), psd_twin(2), main_names{m}));
end


%% PSD per condition x valence — one figure per channel, FEEL | TONE side by side
psd_twin = [0 3];
f     = freq_clean_bl.freq;
t_idx = freq_clean_bl.time >= psd_twin(1) & freq_clean_bl.time <= psd_twin(2);

feel_idx = strcmp(trial_info_clean.Condition, 'FEEL');
tone_idx = strcmp(trial_info_clean.Condition, 'TONE');
neg_idx  = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx  = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx  = strcmp(trial_info_clean.Valence, 'Pos');

channels = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT'  };

main_conditions = {feel_idx, tone_idx};
main_names      = {'FEEL', 'TONE'};

valence_idx       = {neg_idx, neu_idx, pos_idx};
cond_names        = {'Neg', 'Neu', 'Pos'};
cond_colors       = {[0.8 0.1 0.1], [0.4 0.4 0.4], [0.1 0.6 0.2]};
cond_colors_light = {[1.0 0.7 0.7], [0.75 0.75 0.75], [0.7 0.9 0.75]};

for r = 1:size(channels, 1)
    for c = 1:size(channels, 2)
        lab = channels{r,c};
        ch  = find(strcmp(freq_clean_bl.label, lab), 1);

        figure;
        all_ax = gobjects(2,1);

        for m = 1:2
            all_ax(m) = subplot(1, 2, m); hold on;
            h = gobjects(3,1);

            for k = 1:3
                idx = main_conditions{m} & valence_idx{k};

                d      = freq_clean_bl.powspctrm(idx, ch, :, t_idx);
                d_tavg = mean(d, 4);
                d_tavg = reshape(d_tavg, size(d_tavg,1), size(d_tavg,3));

                psd_mean    = mean(d_tavg, 1);
                psd_sem     = std(d_tavg, 0, 1) / sqrt(size(d_tavg, 1));
                upper_bound = psd_mean + psd_sem;
                lower_bound = psd_mean - psd_sem;

                fill([f(:)' fliplr(f(:)')], ...
                     [upper_bound(:)' fliplr(lower_bound(:)')], ...
                     cond_colors_light{k}, 'EdgeColor', 'none', 'FaceAlpha', 0.3);

                h(k) = plot(f(:)', psd_mean(:)', 'Color', cond_colors{k}, 'LineWidth', 2);
            end

            grid on;
            xlabel('Frequency (Hz)'); ylabel('Power (relchange)');
            title(main_names{m}, 'Interpreter', 'none');
            legend(h, cond_names, 'Location', 'best');
            set(gca, 'FontSize', 12);
        end

        % link y axes so both subplots share the same scale
        linkaxes(all_ax, 'y');

        sgtitle(sprintf('PSD %g-%g s: %s', psd_twin(1), psd_twin(2), lab), 'Interpreter', 'none');
    end
end