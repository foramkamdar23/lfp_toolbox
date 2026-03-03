clear all;
clc;
close all;

filePath = 'C:\Users\fkamdar\Desktop\percept\testing_data\Report_Json_Session_Report_20250903T170810.json';

% Read file as text
jsonText = fileread(filePath);

% Convert JSON text into MATLAB struct
data = jsondecode(jsonText);

% Inspect structure
disp(data)

% Save clean MAT file for faster reuse
save('Report_Json_Session_Report_20250903T170810.mat','data','-v7.3')
