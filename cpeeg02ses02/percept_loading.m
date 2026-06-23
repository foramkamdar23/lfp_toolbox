clear all;
clc;
close all;

filePath = 'C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG01\session01\Report_Json_Session_Report_20260312T145645.json';

% Read file as text
jsonText = fileread(filePath);

% Convert JSON text into MATLAB struct
data = jsondecode(jsonText);

% Save clean MAT file for faster reuse
save('Report_Json_Session_Report_20260312T145645.mat','data','-v7.3')


