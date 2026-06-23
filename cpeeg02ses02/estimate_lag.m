function lag_sec = estimate_lag(ref_env, sig_env, fs)
    % Lag (seconds) of sig_env relative to ref_env from the cross-correlation peak.
    [xc, lags] = xcorr(ref_env, sig_env);
    [~, k] = max(xc);
    lag_sec = lags(k) / fs;
end
