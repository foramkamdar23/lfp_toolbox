clear;
clc;

%% ========== CONFIG ==========
rootDir = 'C:\Users\fkamdar\Desktop\repos\data\percept';
skipIfExists = true;
%% ============================

% Find all JSON files recursively
jsonFiles = dir(fullfile(rootDir, '**', '*.json'));

for i = 1:length(jsonFiles)

    try
        % Full path
        jsonPath = fullfile(jsonFiles(i).folder, jsonFiles(i).name);

        % Output MAT path (same folder)
        [path, name, ext] = fileparts(jsonFiles(i).name);
        matPath = fullfile(jsonFiles(i).folder, [name '.mat']);

        % Skip if already exists
        if skipIfExists && exist(matPath, 'file')
            fprintf('[SKIP] %s already exists\n', matPath);
            continue;
        end

        %% Load JSON
        jsonText = fileread(jsonPath);
        data = jsondecode(jsonText);

        %% Save MAT
        save(matPath, 'data', '-v7.3');

        fprintf('  -> Saved: %s\n\n', matPath);

    catch ME
        fprintf('[ERROR] %s\n', jsonFiles(i).name);
        fprintf('  %s\n\n', ME.message);
    end

end