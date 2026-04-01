function corrected_values = check_and_correct_biosemi(status_events, target_code)
% -------------------------------------------------------------------------
% CHECK_AND_CORRECT_BIOSEMI
%
% Ensures that a target trigger/event code exists in a BioSemi STATUS channel.
% If necessary, corrects raw BioSemi STATUS values (subtract 512, mask lower 8 bits)
%
% INPUTS:
%   status_events - structure array from ft_read_event(...), STATUS type
%   target_code   - numeric code to check for (e.g., 141)
%
% OUTPUT:
%   corrected_values - uint16 array of corrected STATUS values
% -------------------------------------------------------------------------

% Convert raw values to numeric (handles string/numeric inconsistencies)
if ischar(status_events(1).value) || isstring(status_events(1).value)
    raw_values = str2double({status_events.value});
else
    raw_values = double([status_events.value]);
end

% Count occurrences in raw data
n_occurrences_raw = sum(raw_values == target_code);

if n_occurrences_raw ~= 8
    fprintf('>>> Target code %d appears %d times in raw events. Attempting correction...\n', ...
            target_code, n_occurrences_raw);

    % Correct BioSemi STATUS values (subtract 512, mask lower 8 bits)
    corrected_values = bitand(uint16(raw_values) - 512, uint16(255));

    n_occurrences_corrected = sum(corrected_values == target_code);

    if n_occurrences_corrected == 8
        fprintf('>>> Target code %d appears %d times after correction ✅\n', ...
                target_code, n_occurrences_corrected);
    else
        error('>>> Target code %d not found correctly after correction ❌', target_code);
    end
else
    fprintf('>>> Target code %d appears %d times in raw events ✅\n', ...
            target_code, n_occurrences_raw);
    corrected_values = uint16(raw_values);  % keep as uint16 for downstream processing
end
end
