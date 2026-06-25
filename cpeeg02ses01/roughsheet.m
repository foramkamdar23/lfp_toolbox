%Plot power spectrum
figure(2);
[pxx,f] = pwelch(lfp, fs*2, [], [], fs);
figure;
plot(f,10*log10(pxx))
xlim([0 100])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title(['BrainSense LFP - ' bs.Channel])


% Convert to absolute time to align with eeg?) 
startTime = datetime(bs.FirstPacketDateTime, ...
    'InputFormat','yyyy-MM-dd''T''HH:mm:ss''Z''', ...
    'TimeZone','UTC');
absoluteTime = startTime + seconds(t);
