folder = "C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG02\ses-01-selected";
files = dir(fullfile(folder, "*.mat"));

for i = 1:length(files)
    
    filePath = fullfile(folder, files(i).name);
    S = load(filePath);
    
    % Check data exists
    if ~isfield(S, 'data')
        fprintf("Skipped (no data): %s\n", files(i).name);
        continue;
    end
    
    data = S.data;
    converted = false;
    
    for j = 1:numel(data)
        
        % Skip if missing
        if ~isfield(data(j), 'IndefiniteStreaming')
            continue;
        end
        
        IS = data(j).IndefiniteStreaming;
        
        % Skip if empty
        if isempty(IS) || ~isfield(IS, 'FirstPacketDateTime')
            continue;
        end
        
        % ✅ Convert once (all channels same timestamp)
        utc_str = IS(1).FirstPacketDateTime;
        
        t_utc = datetime(utc_str, ...
            'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', ...
            'TimeZone','UTC');
        
        t_pacific = datetime(t_utc, 'TimeZone','America/Los_Angeles');
        
        % Assign to all elements
        for k = 1:numel(IS)
            data(j).IndefiniteStreaming(k).FirstPacketDateTime_PDT = t_pacific;
            data(j).IndefiniteStreaming(k).FirstPacketDateTime_PDT_str = ...
                datestr(t_pacific, 'yyyy-mm-dd HH:MM:SS');
        end
        
        converted = true;
    end
    
    % Save only if something was converted
    if converted
        save(filePath, 'data');
        fprintf("Processed: %s\n", files(i).name);
    else
        fprintf("Skipped (no valid timestamps): %s\n", files(i).name);
    end
end