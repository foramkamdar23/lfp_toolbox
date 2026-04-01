function [signal_with_nan, info] = detect_packet_loss_percept(signal, GlobalSequences, GlobalPacketSizes, fs)
%DETECT_PACKET_LOSS_PERCEPT Insert NaNs for missing packets and compute gap info
%
% Inputs:
%   signal             - raw time-domain samples
%   GlobalSequences    - string of packet sequence numbers
%   GlobalPacketSizes  - string of packet sample counts
%   fs                 - sampling rate in Hz (required for gap durations)
%
% Outputs:
%   signal_with_nan    - reconstructed signal with NaNs for missing packets
%   info               - structure with packet loss info and gap info

% Convert strings to numeric arrays
seq  = sscanf(GlobalSequences,'%f,')';
pkts = sscanf(GlobalPacketSizes,'%f,')';

seq  = seq(:);
pkts = pkts(:);

% Detect packet jumps
d = diff(seq);

% Valid steps: +1 or wrap 255->0
expected = (d == 1) | (d == -255);
bad_idx = find(~expected);

% Preallocate output
signal_with_nan = [];
cursor = 1;

missing_packets_total = 0;
missing_samples_total = 0;

for p = 1:length(pkts)

    % Append current packet samples
    packet_samples = pkts(p);
    signal_with_nan = [signal_with_nan; signal(cursor:cursor+packet_samples-1)];
    cursor = cursor + packet_samples;

    % If packet loss occurs after this packet
    if ismember(p, bad_idx)

        jump = seq(p+1) - seq(p);

        if jump < 0
            jump = jump + 256; % handle wrap
        end

        missing_packets = jump - 1;

        if missing_packets > 0
            missing_packets_total = missing_packets_total + missing_packets;

            % Estimate missing samples using packet size
            missing_samples = missing_packets * pkts(p);
            missing_samples_total = missing_samples_total + missing_samples;

            signal_with_nan = [signal_with_nan; nan(missing_samples,1)];
        end
    end
end

% Detect NaN gaps
nan_mask = isnan(signal_with_nan);
d_nan = diff([0; nan_mask; 0]);
start_idx = find(d_nan == 1);
end_idx   = find(d_nan == -1) - 1;

gap_samples = end_idx - start_idx + 1;
gap_ms = (gap_samples / fs) * 1000;  % convert to milliseconds

% Store gap info
info.total_packets      = length(seq);
info.missing_packets    = missing_packets_total;
info.packet_loss_rate   = missing_packets_total / length(seq);
info.missing_samples    = missing_samples_total;
info.bad_packet_indices = bad_idx;
info.gap_info           = table(start_idx, end_idx, gap_samples, gap_ms , ...
                                'VariableNames', {'StartIdx','EndIdx','Samples','Milliseconds'});

end