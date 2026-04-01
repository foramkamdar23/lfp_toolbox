function t_local = perceptUTC2local(utc_time_string, UTC_offset)

t = datetime(utc_time_string, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSS''Z''', ...
    'TimeZone','UTC');

sign_char = UTC_offset(1);
offset_h  = str2double(UTC_offset(2:3));
offset_m  = str2double(UTC_offset(5:6));

offset_hours = offset_h + offset_m/60;

if sign_char == '-'
    offset_hours = -offset_hours;
end

t_local = t + hours(offset_hours);

end