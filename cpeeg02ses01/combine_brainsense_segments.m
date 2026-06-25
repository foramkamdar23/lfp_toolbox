function combined = combine_brainsense_segments(raw_struct)
% Groups entries of a Percept BrainSenseTimeDomain/IndefiniteStreaming struct
% array by Channel name and concatenates multiple segments (e.g. from a
% paused/resumed recording) into one continuous trace per channel. Any gap
% between segments is filled with NaN, sized from FirstPacketDateTime.
%
% If there is already exactly one entry per channel, each channel just
% passes through unchanged.

chan_names = strings(numel(raw_struct), 1);
for i = 1:numel(raw_struct)
    chan_names(i) = string(raw_struct(i).Channel);
end

unique_chans = unique(chan_names, 'stable');
combined = struct('Channel', {}, 'SampleRateInHz', {}, 'TimeDomainData', {}, 'FirstPacketDateTime', {});

for c = 1:numel(unique_chans)
    this_chan = unique_chans(c);
    seg_idx = find(chan_names == this_chan);

    % sort this channel's segments chronologically
    seg_times = datetime(string({raw_struct(seg_idx).FirstPacketDateTime}), ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSX', 'TimeZone', 'UTC');
    [seg_times, order] = sort(seg_times);
    seg_idx = seg_idx(order);

    fs = raw_struct(seg_idx(1)).SampleRateInHz;
    full_trace = raw_struct(seg_idx(1)).TimeDomainData(:)';

    for s = 2:numel(seg_idx)
        prev_end_time = seg_times(s-1) + seconds(numel(raw_struct(seg_idx(s-1)).TimeDomainData) / fs);
        gap_sec = seconds(seg_times(s) - prev_end_time);
        gap_samples = max(round(gap_sec * fs), 0);

        if gap_sec < 0
            warning('Channel %s: segment %d starts before segment %d ends (overlap %.3f s). Check timestamps.', ...
                this_chan, s, s-1, -gap_sec);
        end

        full_trace = [full_trace, nan(1, gap_samples), raw_struct(seg_idx(s)).TimeDomainData(:)']; %#ok<AGROW>
        fprintf('  %s: stitched segment %d, gap = %.3f s (%d samples)\n', this_chan, s, gap_sec, gap_samples);
    end

    combined(end+1).Channel = char(this_chan); %#ok<AGROW>
    combined(end).SampleRateInHz = fs;
    combined(end).TimeDomainData = full_trace;
    combined(end).FirstPacketDateTime = raw_struct(seg_idx(1)).FirstPacketDateTime;
end
end