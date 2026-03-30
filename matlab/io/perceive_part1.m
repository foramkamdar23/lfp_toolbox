addpath(genpath('C:\Users\fkamdar\Documents\MATLAB\perceive-master'));
filepath = 'C:\Users\fkamdar\Desktop\repos\data\percept\cp_nonmotor\CPEEG01\session01\Report_Json_Session_Report_20260312T145645.json';


%% Load data
% INPUTs:
medState = ''; %  Choose from allowed list: "","MedOff01","MedOn01","MedOff02","MedOn02","MedOff03","MedOn03","MedOffOn01"
extended = ''; % ["","yes"]  gives an extensive output of chronic, calibration, lastsignalcheck, diagnostic, impedance and snapshot data
gui = ''; %  ["","yes"]  gives option to skip gui by default settings
% Load data
% perceive(filename, subjectID, medState, extended, gui);

% perceive(filepath, 'cpeeg01') 

path_of_loaded_data = 'C:\Users\fkamdar\Documents\MATLAB\perceive-master\cpeeg01\ses-2026031205361790\ieeg\cpeeg01_ses-2026031205361790_task-Rest_acq-StimOff_mod-ISRing_run-1.mat';
load(path_of_loaded_data)
