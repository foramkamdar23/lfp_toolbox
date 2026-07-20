%%  - LFP loading & preprocessing
% Pipeline: load Percept LFP (two halves) -> align to EEG via TENS 80 Hz sync
% -> merge halves on a common timeline -> select channels of interest
% -> decode BioSemi triggers -> epoch around image onset
% -> attach labels from the curated CSV -> baseline -> split by condition/valence.
%
% Run cell by cell (Run Section). Each section keeps its own settings at the top.
% Helper functions are at the bottom of the file.

clearvars; close all; clc;

%% ===================== SETUP FIELDTRIP =====================
fieldtrip_path = "C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507";

addpath(fieldtrip_path);
ft_defaults;

%% ===================== LOAD PERCEPT LFP =====================
lfp_first_file  = "C:\Users\fkamdar\Desktop\repos\data\ok_new_cp\cpeeg02_ses01_ert1_1437.mat";
lfp_second_file = "C:\Users\fkamdar\Desktop\repos\data\ok_new_cp\cpeeg02_ses01_ert2_1453.mat";
n_lfp_channels  = 6;                 % Percept IndefiniteStreaming channels per file
align_channel   = 'ZERO_THREE_LEFT';  % channel used to estimate the alignment lag

data1 = load(lfp_first_file);
data2 = load(lfp_second_file);

lfp1_all = data1.data.IndefiniteStreaming;
lfp2_all = data2.data.IndefiniteStreaming;  

all_chan_names = strings(n_lfp_channels, 1);
for c = 1:n_lfp_channels
    all_chan_names(c) = string(lfp1_all(c).Channel);
end

align_idx = find(all_chan_names == align_channel, 1);
assert(~isempty(align_idx), 'Alignment channel %s not found in LFP file.', align_channel);

fs_lfp = lfp1_all(align_idx).SampleRateInHz;

sig1 = lfp1_all(align_idx).TimeDomainData(:)';   % alignment channel only, for QC + lag
sig2 = lfp2_all(align_idx).TimeDomainData(:)';
t_lfp1 = (0:numel(sig1)-1) / fs_lfp;
t_lfp2 = (0:numel(sig2)-1) / fs_lfp;

% QC: raw alignment channel, both halves
figure('Name','Raw LFP (alignment channel)');
subplot(2,1,1); plot(t_lfp1, sig1); title(['LFP1: ' char(align_channel)]); xlabel('Time (s)');
subplot(2,1,2); plot(t_lfp2, sig2); title(['LFP2: ' char(align_channel)]); xlabel('Time (s)');

%% ===================== LOAD EEG SYNC CHANNEL =====================
%bdf_file1 = "C:\Users\fkamdar\Desktop\repos\data\ok_new_cp\cpeeg01_ses02_ert_b02.bdf";
bdf_file =  "C:\Users\fkamdar\Desktop\repos\data\ok_new_cp\cpeeg02_ses01_ert_b01.bdf";
sync_channel = "EXG7";               % EEG channel carrying the TENS sync burst

cfg = [];
cfg.dataset = char(bdf_file);
cfg.channel = cellstr(sync_channel);
exg7_data = ft_preprocessing(cfg);

exg7_eeg = exg7_data.trial{1};
t_eeg    = exg7_data.time{1};
fs_eeg   = exg7_data.fsample;

% Downsample EEG sync to the LFP rate so both live on the same sample grid
exg7_ds   = resample(exg7_eeg, fs_lfp, fs_eeg);
t_exg7_ds = (0:numel(exg7_ds)-1) / fs_lfp;

figure('Name','EEG sync downsampling');
subplot(2,1,1); plot(t_eeg, exg7_eeg);    title('Original EEG sync');
subplot(2,1,2); plot(t_exg7_ds, exg7_ds); title(sprintf('Downsampled to %g Hz', fs_lfp));

%% ===================== ESTIMATE ALIGNMENT LAG =====================
% The TENS burst (~80 Hz) appears in both EEG and LFP. Cross-correlating the
% bandpass envelopes gives the time offset between each LFP half and the EEG.
tens_band     = [78 82];   % Hz, bandpass around the 80 Hz TENS burst
tens_filt_ord = 4;

[b_tens, a_tens] = butter(tens_filt_ord, tens_band / (fs_lfp/2), 'bandpass');

exg7_env = bandpass_envelope(exg7_ds, b_tens, a_tens);
lfp1_env = bandpass_envelope(sig1,    b_tens, a_tens);
lfp2_env = bandpass_envelope(sig2,    b_tens, a_tens);

lag_sec1 = estimate_lag(exg7_env, lfp1_env, fs_lfp);
lag_sec2 = estimate_lag(exg7_env, lfp2_env, fs_lfp);

t_lfp1_aligned = t_lfp1 + lag_sec1;
t_lfp2_aligned = t_lfp2 + lag_sec2;
figure('Name','Envelope alignment');
plot(t_exg7_ds, exg7_env); hold on;
plot(t_lfp1_aligned, lfp1_env);
plot(t_lfp2_aligned, lfp2_env);
legend('EEG env','LFP1 env (aligned)','LFP2 env (aligned)');
title('TENS envelope alignment'); xlabel('Time (s)');
fprintf('Alignment done. lag1 = %.3f s, lag2 = %.3f s\n', lag_sec1, lag_sec2);

%% ===================== MERGE LFP HALVES (ALL CHANNELS) =====================
% Build one continuous LFP matrix on the EEG (master) timeline. Gaps between
% the two halves stay as NaN and are visualised below.
t_master = t_exg7_ds;
n_master = numel(t_master);
lfp_merged_all = nan(n_lfp_channels, n_master);

for ch = 1:n_lfp_channels
    s1 = lfp1_all(ch).TimeDomainData(:)';
    s2 = lfp2_all(ch).TimeDomainData(:)';

    idx1 = round(((0:numel(s1)-1)/fs_lfp + lag_sec1) * fs_lfp) + 1;
    idx2 = round(((0:numel(s2)-1)/fs_lfp + lag_sec2) * fs_lfp) + 1;

    v1 = idx1 > 0 & idx1 <= n_master;
    v2 = idx2 > 0 & idx2 <= n_master;

    lfp_merged_all(ch, idx1(v1)) = s1(v1);
    lfp_merged_all(ch, idx2(v2)) = s2(v2);
end

plot_lfp_with_gaps(t_master, lfp_merged_all(align_idx,:));

%% ===================== DECODE BIOSEMI TRIGGERS =====================
% bdf_file is set in the "LOAD EEG SYNC CHANNEL" section above.
drop_first_events = 3;       % leading system events to discard
status_offset     = 512;     % offset removed before masking
trig_map_codes  = [1 2 11 20 12 30 13 40 99 200 201 202];
trig_map_labels = {'intro1','intro2','fix1','cue','fix2','image', ...
                   'fix3','valence_onset','end_screen','break_start', ...
                   'resync_start','resync_end'};

events = ft_read_event(char(bdf_file));
events(1:drop_first_events) = [];

trig = double(bitand(int32([events.value]' - status_offset), 255));
trig_time = [events.sample]' / exg7_data.hdr.Fs;

trig_table = table(trig, trig_time, 'VariableNames', {'TriggerCode','TriggerTime'});
trig_table = trig_table(trig_table.TriggerCode ~= 0, :);

trig_map = containers.Map(trig_map_codes, trig_map_labels);
event_labels = strings(height(trig_table),1);
for i = 1:height(trig_table)
    code = trig_table.TriggerCode(i);
    if isKey(trig_map, code)
        event_labels(i) = trig_map(code);
    elseif code >= 51 && code <= 59          % valence responses 1-9
        event_labels(i) = "valence_resp_" + string(code - 50);
    else
        event_labels(i) = "unknown";
    end
end
trig_table.Event = event_labels;

figure('Name','LFP + fix1 triggers');
plot(t_master, lfp_merged_all(align_idx,:)); hold on;
fix1_times = trig_table.TriggerTime(trig_table.Event == "fix1");
for i = 1:numel(fix1_times)
    xline(fix1_times(i), 'r');
end
title('Merged LFP + fix1 triggers'); ylabel('LFP'); xlabel('Time (s)');

%% ===================== BUILD FIELDTRIP RAW & SELECT CHANNELS =====================
% Channels of interest and the plot layout (left column = LEFT hemisphere,
% right column = RIGHT). plot_layout is reused by the plotting sections below.
channels_of_interest = {'ZERO_TWO_LEFT','ONE_THREE_RIGHT','ZERO_TWO_RIGHT'};
% labels = {'ZERO_TWO_LEFT' : 'ventral channel left hemisphere',
% 'ONE_THREE_RIGHT': 'dorsal channel right hemisphere' , 'ZERO_TWO_RIGHT':
% 'ventral channel right hemisphere'

plot_layout = { 'ZERO_TWO_LEFT','ONE_THREE_RIGHT','ZERO_TWO_RIGHT'};

channel_names = cellstr(all_chan_names(:));

data_ft = [];
data_ft.label      = channel_names;
data_ft.fsample    = fs_lfp;
data_ft.trial      = {lfp_merged_all};
data_ft.time       = {t_master};
data_ft.sampleinfo = [1 n_master];
data_ft.hdr.label    = channel_names;
data_ft.hdr.Fs       = fs_lfp;
data_ft.hdr.nChans   = n_lfp_channels;
data_ft.hdr.nSamples = n_master;
data_ft.hdr.nTrials  = 1;

data_ft = ft_checkdata(data_ft, 'datatype','raw','hassampleinfo','yes');

cfg = [];
cfg.channel = channels_of_interest;
data_ft = ft_selectdata(cfg, data_ft);
fprintf('FieldTrip raw structure created with channels: %s\n', strjoin(data_ft.label, ', '));

%% ===================== EPOCH AROUND IMAGE ONSET =====================
% 4th trl column = original image-trial number (1..N), stored in trialinfo.
n_trials_exp = 96;     % image trials presented (48 + 48 across the break)
epoch_pre_s  = 4;      % seconds before image onset (image presentation = 3 s)
epoch_post_s = 5;      % seconds after image onset

image_samples = round(trig_table.TriggerTime(trig_table.Event == "image") * fs_lfp);
fprintf('Found %d image events (expected %d)\n', numel(image_samples), n_trials_exp);

pre  = round(epoch_pre_s  * fs_lfp);
post = round(epoch_post_s * fs_lfp);

trl = [];
for i = 1:numel(image_samples)
    s   = image_samples(i);
    beg = s - pre;
    fin = s + post;
    if beg > 0 && fin <= n_master
        trl = [trl; beg fin -pre i];     %#ok<AGROW>  4th col = original trial #
    else
        fprintf('  trial %d skipped (window outside recording)\n', i);
    end
end

cfg = [];
cfg.trl = trl;
data_epoched = ft_redefinetrial(cfg, data_ft);
fprintf('Epoching complete: %d epochs\n', numel(data_epoched.trial));

%% ===================== QC FLAGS (TENS + NaN, advisory only) =====================
% These no longer drop trials. The curated CSV is the source of truth for which
% trials are kept (next section). QC flags are reported for cross-check.
% b_tens / a_tens come from the "ESTIMATE ALIGNMENT LAG" section.
tens_reject_sd = 2;    % flag trials with 80 Hz power > mean + N*SD

n_epochs = numel(data_epoched.trial);
orig_all = data_epoched.trialinfo(:,1);
power_80hz = zeros(n_epochs,1);
nan_bad    = false(n_epochs,1);

for i = 1:n_epochs
    td  = data_epoched.trial{i};                          % [channels x time]
    env = abs(hilbert(filtfilt(b_tens, a_tens, td')'));   % per channel
    power_80hz(i) = mean(env, 'all');
    nan_bad(i)    = any(isnan(td), 'all');
end
tens_bad = power_80hz > (mean(power_80hz) + tens_reject_sd*std(power_80hz));
qc_flag  = tens_bad | nan_bad;

qc_info = table(orig_all, tens_bad, nan_bad, qc_flag, ...
                'VariableNames', {'OrigTrial','TENS_Flag','NaN_Flag','QC_Flag'});
fprintf('QC flags (advisory): %d TENS, %d NaN, %d total flagged\n', ...
        sum(tens_bad), sum(nan_bad), sum(qc_flag));

%% ===================== ATTACH LABELS FROM CURATED CSV =====================
%labels_csv      = "C:\Users\fkamdar\Desktop\repos\lfp_toolbox\cpeeg02_block_96.csv";
labels_csv      = "C:\Users\fkamdar\Desktop\repos\lfp_toolbox\cpeeg02ses01\cpeeg02_block_96_neu4-6.csv"
csv_trialid_col = "trial";          % column holding the original 1..96 trial number
csv_cond_col    = "condition";
csv_val_col     = "Valence_group";

csv_data = readtable(char(labels_csv));
csv_vars = string(csv_data.Properties.VariableNames);

% resolve the trial-id column (auto-detect if not set)
id_candidates = ["trial","trial_num","trial_number","trial_idx","trial_id","OrigTrial","index"];
if strlength(csv_trialid_col) > 0
    idcol = csv_trialid_col;
else
    hit = id_candidates(ismember(id_candidates, csv_vars));
    idcol = ""; if ~isempty(hit), idcol = hit(1); end
end

if strlength(idcol) > 0
    % ---- safe path: match epochs to CSV by original trial number ----
    csv_ids = csv_data.(char(idcol));
    keep    = ismember(orig_all, csv_ids);
    sel     = find(keep);

    cfg = []; cfg.trials = sel;
    data_clean = ft_selectdata(cfg, data_epoched);

    [~, loc] = ismember(orig_all(keep), csv_ids);   % CSV row for each kept epoch
    trial_info_clean = table(orig_all(keep), ...
        string(csv_data.(char(csv_cond_col))(loc)), ...
        string(csv_data.(char(csv_val_col))(loc)), ...
        'VariableNames', {'OrigTrial','Instruction','Valence'});

    % ---- reconciliation report ----
    csv_no_epoch     = setdiff(csv_ids, orig_all);    % in CSV, no matching epoch
    kept_but_flagged = orig_all(keep & qc_flag);      % kept despite QC flag
    fprintf('Labeling by trial id "%s": kept %d epochs.\n', idcol, numel(sel));
    if ~isempty(csv_no_epoch)
        fprintf('  WARNING: %d CSV trials have no matching epoch: %s\n', ...
                numel(csv_no_epoch), mat2str(csv_no_epoch(:)'));
    end
    if ~isempty(kept_but_flagged)
        fprintf('  NOTE: %d kept trials were QC-flagged: %s\n', ...
                numel(kept_but_flagged), mat2str(kept_but_flagged(:)'));
    end
else
    % ---- fallback: positional match, only if counts agree ----
    sel = find(~qc_flag);
    assert(height(csv_data) == numel(sel), ...
        ['No trial-id column found in CSV and counts disagree ' ...
         '(CSV=%d, pipeline-kept=%d). Add a trial-number column to the CSV ' ...
         'or set csv_trialid_col.'], height(csv_data), numel(sel));
    warning(['No trial-id column found. Matching CSV rows to kept epochs by ' ...
             'POSITION. Only valid if manual deletions match the QC-flagged trials.']);
    cfg = []; cfg.trials = sel;
    data_clean = ft_selectdata(cfg, data_epoched);
    trial_info_clean = table(orig_all(sel), ...
        string(csv_data.(char(csv_cond_col))), ...
        string(csv_data.(char(csv_val_col))), ...
        'VariableNames', {'OrigTrial','Instruction','Valence'});
end

disp(trial_info_clean);

%% ===================== BASELINE CORRECTION =====================
baseline_win = [-0.5 0];    % seconds, pre-image baseline

cfg = [];
cfg.demean         = 'yes';
cfg.baselinewindow = baseline_win;
data_bl = ft_preprocessing(cfg, data_clean);

% plot_layout is defined in "BUILD FIELDTRIP RAW & SELECT CHANNELS".
plot_lfp_trials(data_clean, plot_layout, 'Clean trials (trials + mean)');
plot_lfp_trials(data_bl,    plot_layout, 'Baseline-corrected trials');

%% ===================== SPLIT BY CONDITION & VALENCE =====================
data_feel = ft_selectdata(struct('trials', find(trial_info_clean.Instruction == "FEEL")), data_clean);
data_tone = ft_selectdata(struct('trials', find(trial_info_clean.Instruction == "TONE")), data_clean);

data_pos = ft_selectdata(struct('trials', find(trial_info_clean.Valence == "Pos")), data_clean);
data_neu = ft_selectdata(struct('trials', find(trial_info_clean.Valence == "Neu")), data_clean);
data_neg = ft_selectdata(struct('trials', find(trial_info_clean.Valence == "Neg")), data_clean);

fprintf('FEEL=%d TONE=%d | Pos=%d Neu=%d Neg=%d\n', ...
    sum(trial_info_clean.Instruction=="FEEL"), sum(trial_info_clean.Instruction=="TONE"), ...
    sum(trial_info_clean.Valence=="Pos"), sum(trial_info_clean.Valence=="Neu"), ...
    sum(trial_info_clean.Valence=="Neg"));

%% ===================== PLOT CONDITION & VALENCE AVERAGES =====================
avg_feel = ft_timelockanalysis([], data_feel);
avg_tone = ft_timelockanalysis([], data_tone);
plot_avg_overlay({avg_feel, avg_tone}, {'b','r'}, ...
                 {'Feel (perception)','Tone (regulation)'}, plot_layout, 'Condition averages');

avg_neg = ft_timelockanalysis([], data_neg);
avg_neu = ft_timelockanalysis([], data_neu);
avg_pos = ft_timelockanalysis([], data_pos);
plot_avg_overlay({avg_neg, avg_neu, avg_pos}, {'r','g','b'}, ...
                 {'Neg','Neu','Pos'}, plot_layout, 'Valence averages');



subject = 'CPEEG02';  session = 'ses-01';
deriv_root = 'C:\Users\fkamdar\Desktop\repos\lfp_toolbox\derivatives';
if ~exist(deriv_root, 'dir'); mkdir(deriv_root); end
save(fullfile(deriv_root, sprintf('sub-%s_%s_data_clean_neu4-6.mat', subject, session)), ...
     'data_clean', 'trial_info_clean', '-v7.3');