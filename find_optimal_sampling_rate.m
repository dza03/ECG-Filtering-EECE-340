function [optimal_fs, error_curve] = find_optimal_sampling_rate(signal, t, max_fs, target_error)
    % Find optimal sampling rate for a given signal
    % Inputs:
    %   signal - original continuous signal
    %   t - time vector
    %   max_fs - maximum sampling frequency to consider
    % Outputs:
    %   optimal_fs - recommended sampling frequency
    %   error_curve - error vs. sampling rate data
    
    if nargin < 4
        target_error = 0.01;  % 1% error default
    end
    
    % Calculate original sampling rate
    fs_original = 1/(t(2)-t(1));
    
    % Estimate bandwidth more robustly. When running the code
    % ECG signals can have unusual spectral distributions that don't reach the 95% energy threshold in a clean way.
    %The fallback mechanism sets the bandwidth to another frequency.
    % This highlights the importance of combining signal processing theory
    % with domain expertise when designing systems for biomedical signals.
    NFFT = length(signal);
    f = (-NFFT/2:NFFT/2-1)*(fs_original/NFFT);
    X_f = abs(fftshift(fft(signal, NFFT)));
    
    % Only consider positive frequencies for bandwidth estimation
    pos_freq_idx = f >= 0;
    f_pos = f(pos_freq_idx);
    X_f_pos = X_f(pos_freq_idx);
    
    % Calculate cumulative energy (as percentage)
    if sum(X_f_pos.^2) > 0  % Make sure we have energy
        X_f_sq = X_f_pos.^2;
        cum_energy = cumsum(X_f_sq) / sum(X_f_sq) * 100;
        
        % Find 95% energy bandwidth (with safeguard)
        bw_idx = find(cum_energy >= 95, 1, 'first');
        
        % If we couldn't find 95%, use 90% or whatever we can find
        if isempty(bw_idx)
            bw_idx = find(cum_energy >= 90, 1, 'first');
            if isempty(bw_idx)
                bw_idx = length(f_pos);  % Use max frequency as fallback
            end
        end
        
        signal_bw = f_pos(bw_idx);
    else
        % Fallback if spectrum calculation fails
        signal_bw = fs_original / 4;
    end
    
    % Ensure non-zero bandwidth with fallback
    if signal_bw <= 0 || isnan(signal_bw)
        fprintf('Warning: Bandwidth estimation failed, using heart rate-based estimate.\n');
        
        %ECG-specific information suggested by ChatGPT
        % Estimate heart rate from signal
        [peaks, locs] = findpeaks(signal, 'MinPeakHeight', max(signal)*0.6, 'MinPeakDistance', round(fs_original*0.3));
        if length(locs) >= 2
            % Calculate heart rate from peak intervals
            mean_interval = mean(diff(locs))/fs_original;  % in seconds
            heart_rate = 60/mean_interval;  % in bpm
            fprintf('Estimated heart rate: %.1f bpm\n', heart_rate);
            
            % Adjust bandwidth based on heart rate
            % Higher heart rates have higher frequency content
            signal_bw = max(40, heart_rate/3 + 25);  % Formula: heart rate effect on bandwidth
        else
            % Fallback if peak detection fails
            signal_bw = 45;  % Default for normal ECG
        end
        
        fprintf('Using heart rate-adjusted bandwidth: %.1f Hz\n', signal_bw);
    end
    
    % Theoretical minimum (Nyquist)
    nyquist_fs = 2 * signal_bw;
    
    % Test range of sampling rates
    fs_range = linspace(max(nyquist_fs*0.5, 10), min(max_fs, nyquist_fs*4), 20);
    errors = zeros(size(fs_range));
    
    % Test each sampling rate
    for i = 1:length(fs_range)
        fs = fs_range(i);
        
        % Sample and reconstruct
        [t_sampled, x_sampled] = sample(t, signal, fs);
        x_recon = reconstruct(t, x_sampled, t_sampled, fs);
        
        % Calculate error
        errors(i) = norm(signal - x_recon) / norm(signal);
    end
    
    % Find optimal rate (lowest rate with error below target)
    error_idx = find(errors <= target_error, 1, 'first');
    if isempty(error_idx)
        optimal_fs = max(fs_range);
        warning('Target error not achieved. Consider higher sampling rates.');
    else
        optimal_fs = fs_range(error_idx);
    end
    
    error_curve = [fs_range; errors];
    
    % Plot error curve
    figure;
    plot(fs_range, errors, 'o-', 'LineWidth', 2);
    grid on;
    xlabel('Sampling Frequency (Hz)');
    ylabel('Relative Reconstruction Error');
    title('Error vs. Sampling Rate');
    
    % Mark theoretical Nyquist rate
    hold on;
    plot([nyquist_fs, nyquist_fs], [0, max(errors)], 'r--');
    text(nyquist_fs, max(errors)*0.8, 'Nyquist Rate', 'Color', 'r');
    
    % Mark optimal rate
    [~, closest_idx] = min(abs(fs_range - optimal_fs));
    plot(fs_range(closest_idx), errors(closest_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(fs_range(closest_idx), errors(closest_idx), sprintf('Optimal: %.2f Hz', optimal_fs));
    
    % Display results
    fprintf('Signal analysis results:\n');
    fprintf('Estimated bandwidth: %.2f Hz\n', signal_bw);
    fprintf('Theoretical Nyquist rate: %.2f Hz\n', nyquist_fs);
    fprintf('Recommended sampling rate: %.2f Hz (%.1fx Nyquist)\n', optimal_fs, optimal_fs/nyquist_fs);
    
    % Get the error at optimal rate (or closest available rate)
    fprintf('Expected reconstruction error: %.4f\n', errors(closest_idx));
end