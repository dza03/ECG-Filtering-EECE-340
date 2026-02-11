addpath('C:\Users\dinaa\OneDrive\Desktop\EECE 340_PROJ');

%   We need to load our synthetic ECG data
    load('ecg_data.mat');
    
    % This layout was taken from AI upon asking on how to implement the
    % comparison of the normal and tachycardic data
    ecg_type = input('Select ECG type (1: Normal, 2: Tachycardia): ');
    
    if ecg_type == 2
        xt = tachy_ecg;
        t = t_tachy;
        disp('Analyzing tachycardia ECG');
    else
        xt = normal_ecg;
        t = t_normal;
        disp('Analyzing normal ECG');
    end
    
    fs = 1/(t(2)-t(1));


% Extract a meaningful segment
segment_duration = 3;  % seconds
segment_samples = round(segment_duration * fs);
if length(xt) > segment_samples
    % Choose different starting points to capture representative segments
    if ecg_type == 2
        start_idx = 1000;  % Starting sample for tachycardia
    else
        start_idx = 800;   % Different starting sample for normal ECG
    end
    
    % Extract segment with proper time alignment
    xt = xt(start_idx:start_idx+segment_samples-1);
    t = t(start_idx:start_idx+segment_samples-1);  % Keep time aligned with signal
    
    % Plot the extracted segment to verify
    figure;
    plot(t, xt);
    if ecg_type == 2
        title('Extracted Tachycardia ECG Segment');
    else
        title('Extracted Normal ECG Segment');
    end
    xlabel('Time (s)'); ylabel('Amplitude');
    grid on;
end
% Apply lowpass filter
cutoff_freq = 40;  %40 Hz cutoff
filter_order = 101; % odd number for symmetry (Also recommended by AI)
h = design_lowpass_fir(cutoff_freq, fs, filter_order);

xt_filtered = conv(xt, h, 'same');

% Compare original and filtered signals in time domain
figure;
subplot(2,1,1);
plot(t, xt, 'b');
title('Original Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(t, xt_filtered, 'r');
title('Filtered Signal (Time Domain)');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

% Compare original and filtered signals in frequency domain
NFFT = length(xt);
f_axis = (-NFFT/2:NFFT/2-1)*(fs/NFFT);
X_f = abs(fftshift(fft(xt, NFFT)));
Xf_filtered = abs(fftshift(fft(xt_filtered, NFFT)));

figure;
subplot(2,1,1);
plot(f_axis, X_f);
title('Original Signal');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;
xlim([0 fs/2]);

subplot(2,1,2);
plot(f_axis, Xf_filtered, 'r');
title('Filtered Signal');
xlabel('Frequency (Hz)'); ylabel('Magnitude'); grid on;
xlim([0 fs/2]);

% Plot filter characteristics
figure;
subplot(2,1,1);
stem((0:length(h)-1) / fs, h);
title('FIR Filter Impulse Response');
xlabel('Time (s)'); ylabel('h[n]'); grid on;

H_f = abs(fftshift(fft(h, NFFT)));
f_h_axis = (-NFFT/2:NFFT/2-1)*(fs/NFFT);

subplot(2,1,2);
plot(f_h_axis, H_f);
title('FIR Filter Frequency Response');
xlabel('Frequency (Hz)'); ylabel('|H(f)|'); grid on;
xlim([0 fs/2]);

% Calculate Z-transform of the FIR filter
syms z;
H_z = 0;
for k = 1:length(h)
    H_z = H_z + h(k) * z^(-(k-1));
end
disp('System function H(z):');
pretty(H_z)

% Comment on stability and causality
disp('Stability analysis: Stable since it has a finite impulse response.');
if all(abs(roots(poly(h))) < 1)
    disp('All poles are within the unit circle, confirming stability.');
end

% Check causality
if mod(filter_order, 2) == 1
    disp('Causality: The filter is non-causal');

else
    disp('Causality: The filter may not have perfect linear phase.');
end


% Calculate optimal sampling rate for filtered signal from extra credit
% function
%This actually results in a warning on tachycardic data which is great, because it confirms our hypothesis that tachycardia signals
% require higher sampling rates for accurate reconstruction. 
[optimal_fs, ~] = find_optimal_sampling_rate(xt_filtered, t, fs/2, 0.01);

%Sample at various rates
nyquist_fs = 2 * cutoff_freq;  % Theoretical Nyquist rate based on filter cutoff
fs_under = nyquist_fs * 0.5;   % Undersampling (below Nyquist)
fs_nyquist = nyquist_fs;       % At Nyquist rate
fs_over = nyquist_fs * 2;      % Oversampling
fs_optimal = optimal_fs;       % Optimal rate from our function

sampling_rates = [fs_under, fs_nyquist, fs_over, fs_optimal];
rate_names = {'Undersampling', 'Nyquist Rate', 'Oversampling', 'Optimal Rate'};
errors = zeros(length(sampling_rates), 1);

for i = 1:length(sampling_rates)
    fs_sample = sampling_rates(i);
    
    % Sample filtered signal
    [t_sample, x_sample] = sample(t, xt_filtered, fs_sample);
    
    % Reconstruct
    x_recon = reconstruct(t, x_sample, t_sample, fs_sample);
    
    % Calculate error
    errors(i) = norm(xt_filtered - x_recon) / norm(xt_filtered);
    
    % Plot reconstruction
    subplot(5,1,i+1);
    plot(t, xt_filtered, 'k', t, x_recon, 'r--');
    legend('Original', 'Reconstructed');
    title(sprintf('%s (%.1f Hz): Error = %.4f', rate_names{i}, fs_sample, errors(i)));
    xlabel('Time (s)'); ylabel('Amplitude'); grid on;
end

% Plot Error results (Extra) 


%Noise robustness analysis
noise_levels_dB = [30, 20, 10, 0];  % SNR in dB (higher = less noise)
fs_sample = optimal_fs;  % Use the optimal sampling rate

figure;
subplot(length(noise_levels_dB)+1, 1, 1);
plot(t, xt_filtered);
title('Original Filtered Signal (No Noise)');
grid on;

errors_noise = zeros(length(noise_levels_dB), 1);

for i = 1:length(noise_levels_dB)
    snr_dB = noise_levels_dB(i);
    
    % Add noise to filtered signal (manual implementation of white gaussian
    % noise from AI
    signal_power = sum(xt_filtered.^2) / length(xt_filtered);
    noise_power = signal_power / (10^(snr_dB/10));
    noise = sqrt(noise_power) * randn(size(xt_filtered));
    xt_noisy = xt_filtered + noise;
        
    % Sample noisy signal
    [t_sample_noisy, x_sample_noisy] = sample(t, xt_noisy, fs_sample);
    
    % Reconstruct from noisy samples
    x_recon_noisy = reconstruct(t, x_sample_noisy, t_sample_noisy, fs_sample);
    
    % Calculate error
    errors_noise(i) = norm(xt_filtered - x_recon_noisy) / norm(xt_filtered);
    
    % Plot
    subplot(length(noise_levels_dB)+1, 1, i+1);
    plot(t, xt_filtered, 'k', t, x_recon_noisy, 'r--');
    title(sprintf('SNR = %d dB, Error = %.4f', snr_dB, errors_noise(i)));
    grid on;
end

% Plot error vs noise level
figure;
plot(noise_levels_dB, errors_noise, 'o-', 'LineWidth', 2);
set(gca, 'XDir', 'reverse');  % Higher SNR (less noise) on left
xlabel('Noise Level (SNR in dB)');
ylabel('Reconstruction Error');
title('Effect of Noise on Reconstruction Error');
grid on;

%The following was also recommended by AI 
fprintf('\nNoise robustness analysis:\n');
fprintf('As noise increases (SNR decreases), reconstruction error increases.\n');
fprintf('At SNR = %.0f dB, error = %.4f\n', noise_levels_dB(1), errors_noise(1));
fprintf('At SNR = %.0f dB, error = %.4f\n', noise_levels_dB(end), errors_noise(end));
fprintf('This demonstrates the importance of signal-to-noise ratio in sampling systems.\n');

% Summary of results
fprintf('\nSummary of filter and sampling analysis:\n');
fprintf('Filter cutoff frequency: %.1f Hz\n', cutoff_freq);
fprintf('Original sampling frequency: %.1f Hz\n', fs);
fprintf('Filter order: %d taps\n', filter_order);
fprintf('Theoretical Nyquist rate: %.1f Hz\n', nyquist_fs);
fprintf('Optimal sampling rate: %.1f Hz (%.1fx Nyquist)\n', optimal_fs, optimal_fs/nyquist_fs);
fprintf('This analysis demonstrates how proper filter design and sampling rate selection\n');
fprintf('affects reconstruction quality, particularly in the presence of noise.\n');