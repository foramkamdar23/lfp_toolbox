
% Lets look into the triggers now
addpath('C:\Users\fkamdar\Documents\MATLAB\fieldtrip-20210507');
ft_defaults;


events = ft_read_event('R:\Trainee_folders\foram\Testdata.bdf');
events(1:3) = [];
trig = [events.value]'; 
trig = trig-512;
trig = double(bitand(int32(trig),255)); % we have triggers now
