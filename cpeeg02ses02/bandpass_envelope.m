function z = bandpass_envelope(sig, b, a)
    % Bandpass -> Hilbert envelope -> zscore. Returns a row vector.
    sig = sig(:)';
    z = zscore(abs(hilbert(filtfilt(b, a, sig))));
end

