fieldtrip_path = "C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507";

addpath(fieldtrip_path);
ft_defaults;

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
cfg.foi = 1:1:99;
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
 


% function plot_tfr_grid(freq, layout, clim, ttl)
%     % Per-channel TFR grid (LEFT col = left hemisphere, RIGHT col = right).
%     % freq.powspctrm is [channels x freq x time] after averaging over trials.
%     [nr, nc] = size(layout);
%     figure('Name', ttl);
%     for r = 1:nr
%         for c = 1:nc
%             lab = layout{r,c};
%             ch = find(strcmp(freq.label, lab), 1);
%             subplot(nr, nc, (r-1)*nc + c);
%             imagesc(freq.time, freq.freq, squeeze(freq.powspctrm(ch,:,:)), clim);
%             axis xy;
%             xline(0, 'k--');                 % image onset
%             title(lab, 'Interpreter', 'none');
%             xlabel('Time (s)'); ylabel('Frequency (Hz)');
%             colorbar;
%         end
%     end
%     colormap jet
%     sgtitle(ttl);
% end
% 


%% TFR by valence group (Neg / Neu / Pos)
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx = strcmp(trial_info_clean.Valence, 'Pos');

valence_idx   = {neg_idx,  neu_idx,  pos_idx};
valence_names = {'Neg',    'Neu',    'Pos'};

for v = 1:3
    % select trials for this valence
    cfg = [];
    cfg.trials = find(valence_idx{v});
    freq_val = ft_selectdata(cfg, freq_clean_bl);

    % average over trials
    cfg = [];
    cfg.avgoverrpt = 'yes';
    freq_val_avg = ft_selectdata(cfg, freq_val);

    plot_tfr_grid(freq_val_avg, plot_layout, tfr_clim, ...
        sprintf('TFR (%s baseline %g to %g s) — %s', ...
                tfr_basetype, tfr_baseline(1), tfr_baseline(2), valence_names{v}));
end

%% TFR by valence group — one figure per channel, 3 valence subplots side by side
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx = strcmp(trial_info_clean.Valence, 'Pos');

valence_idx   = {neg_idx,  neu_idx,  pos_idx};
valence_names = {'Neg',    'Neu',    'Pos'};

channels = {'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT'  };

% pre-compute per-valence averages once
freq_val_avg = cell(3,1);
for v = 1:3
    cfg        = [];
    cfg.trials = find(valence_idx{v});
    freq_val   = ft_selectdata(cfg, freq_clean_bl);

    cfg            = [];
    cfg.avgoverrpt = 'yes';
    freq_val_avg{v} = ft_selectdata(cfg, freq_val);
end

for r = 1:size(channels, 1)
    for c = 1:size(channels, 2)
        lab = channels{r,c};
        ch  = find(strcmp(freq_clean_bl.label, lab), 1);

        figure('Name', lab);
        all_ax = gobjects(3,1);

        for v = 1:3
            all_ax(v) = subplot(1, 3, v);
            imagesc(freq_val_avg{v}.time, freq_val_avg{v}.freq, ...
                    squeeze(freq_val_avg{v}.powspctrm(ch,:,:)), tfr_clim);
            axis xy;
            xline(0, 'k--');
            title(valence_names{v});
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            colorbar;
        end

        colormap jet;
        linkaxes(all_ax, 'xy');
        sgtitle(sprintf('TFR (%s) — %s', tfr_basetype, lab), 'Interpreter', 'none');
    end
end

%% TFR by valence group — one figure per channel, 3 valence subplots + 1 empty (2x2)
neg_idx = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx = strcmp(trial_info_clean.Valence, 'Pos');

valence_idx   = {neg_idx,  neu_idx,  pos_idx};
valence_names = {'Neg',    'Neu',    'Pos'};

channels = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT'  };

% pre-compute per-valence averages once
freq_val_avg = cell(3,1);
for v = 1:3
    cfg        = [];
    cfg.trials = find(valence_idx{v});
    freq_val   = ft_selectdata(cfg, freq_clean_bl);

    cfg             = [];
    cfg.avgoverrpt  = 'yes';
    freq_val_avg{v} = ft_selectdata(cfg, freq_val);
end

for r = 1:size(channels, 1)
    for c = 1:size(channels, 2)
        lab = channels{r,c};
        ch  = find(strcmp(freq_clean_bl.label, lab), 1);

        figure('Name', lab);
        all_ax = gobjects(3,1);

        for v = 1:3
            all_ax(v) = subplot(2, 2, v);
            imagesc(freq_val_avg{v}.time, freq_val_avg{v}.freq, ...
                    squeeze(freq_val_avg{v}.powspctrm(ch,:,:)), tfr_clim);
            axis xy;
            xline(0, 'k--');
            title(valence_names{v});
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            colorbar;
            ylim([1 99]);
            xlim([-0.5 3]);
        end

        % 4th subplot: empty
        subplot(2, 2, 4);
        axis off;

        colormap jet;
        linkaxes(all_ax, 'xy');
        sgtitle(sprintf('TFR (%s) — %s', tfr_basetype, lab), 'Interpreter', 'none');
    end
end

%% TFR by condition x valence — one figure per channel per condition (2x2: Neg/Neu/Pos/empty)
feel_idx = strcmp(trial_info_clean.Condition, 'FEEL');
tone_idx = strcmp(trial_info_clean.Condition, 'TONE');
neg_idx  = strcmp(trial_info_clean.Valence, 'Neg');
neu_idx  = strcmp(trial_info_clean.Valence, 'Neu');
pos_idx  = strcmp(trial_info_clean.Valence, 'Pos');

valence_idx   = {neg_idx,  neu_idx,  pos_idx};
valence_names = {'Neg',    'Neu',    'Pos'};

main_conditions = {feel_idx, tone_idx};
main_names      = {'FEEL',   'TONE'};

channels = { 'ZERO_TWO_LEFT',  'ZERO_TWO_RIGHT'  };

% pre-compute per-condition x valence averages once
% freq_val_avg{m,v}: condition m, valence v
freq_val_avg = cell(2, 3);
for m = 1:2
    for v = 1:3
        cfg        = [];
        cfg.trials = find(main_conditions{m} & valence_idx{v});
        freq_val   = ft_selectdata(cfg, freq_clean_bl);

        cfg             = [];
        cfg.avgoverrpt  = 'yes';
        freq_val_avg{m,v} = ft_selectdata(cfg, freq_val);
    end
end

for m = 1:2
    for r = 1:size(channels, 1)
        for c = 1:size(channels, 2)
            lab = channels{r,c};
            ch  = find(strcmp(freq_clean_bl.label, lab), 1);

            figure('Name', sprintf('%s — %s', main_names{m}, lab));
            all_ax = gobjects(3,1);

            for v = 1:3
                all_ax(v) = subplot(2, 2, v);
                imagesc(freq_val_avg{m,v}.time, freq_val_avg{m,v}.freq, ...
                        squeeze(freq_val_avg{m,v}.powspctrm(ch,:,:)), tfr_clim);
                axis xy;
                xline(0, 'k--');
                title(valence_names{v});
                xlabel('Time (s)'); ylabel('Frequency (Hz)');
                colorbar;
                ylim([1 99]);
                xlim([-0.5 3]);
            end

            % 4th subplot: empty
            subplot(2, 2, 4);
            axis off;

            colormap jet;
            linkaxes(all_ax, 'xy');
            sgtitle(sprintf('TFR (%s) | %s — %s', tfr_basetype, main_names{m}, lab), ...
                    'Interpreter', 'none');
        end
    end
end