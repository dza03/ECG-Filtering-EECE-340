% Examine noise robustness
t_noise = 0:0.001:1;       % Time vector
f1 = 3;                    % First frequency (Hz)
f2 = 7;                    % Second frequency (Hz)
fmax = max(f1, f2);        % Maximum frequency

% Original signal
xt_noise = cos(2*pi*f1*t_noise) + 0.5*cos(2*pi*f2*t_noise);

% Try different sampling rates
fs_values = [1.5*fmax, 2*fmax, 4*fmax];
noise_levels = [0.01, 0.05, 0.1, 0.2];  % Different noise standard deviations

% Initialize error matrix
error_matrix = zeros(length(fs_values), length(noise_levels));

% Create figure for visualization
figure;
subplot_count = 1;

for i = 1:length(fs_values)
    fs = fs_values(i);
    
    % Sample the clean signal
    [t_sample_clean, x_sample_clean] = sample(t_noise, xt_noise, fs);
    
    for j = 1:length(noise_levels)
        noise_level = noise_levels(j);
        
        % Add Gaussian noise to the samples
        x_sample_noisy = x_sample_clean + noise_level * randn(size(x_sample_clean));
        
        % Reconstruct from noisy samples
        t_recon_noise = 0:0.0005:1;  % Fine grid for reconstruction
        xrcon_noisy = reconstruct(t_recon_noise, x_sample_noisy, t_sample_clean, fs);
        
        % Calculate reconstruction error
        xt_interp = interp1(t_noise, xt_noise, t_recon_noise, 'linear');
        error_matrix(i, j) = mean(abs(xt_interp - xrcon_noisy).^2);
        
        % Plot one example for each sampling rate and highest noise level
        if j == length(noise_levels)
            subplot(length(fs_values), 1, i);
            plot(t_noise, xt_noise, 'k-'); hold on;
            plot(t_sample_clean, x_sample_noisy, 'ro', 'MarkerSize', 6);
            plot(t_recon_noise, xrcon_noisy, 'b--', 'LineWidth', 1.5);
            title(['fs = ' num2str(fs) ' Hz, Noise Level = ' num2str(noise_level)]);
            xlabel('Time (s)'); ylabel('Amplitude'); grid on;
            legend('Original', 'Noisy Samples', 'Reconstructed');
        end
    end
end

% Plot error matrix as a heatmap
figure;
imagesc(noise_levels, fs_values/fmax, error_matrix);
colorbar;
title('Reconstruction Error vs. Sampling Rate and Noise Level');
xlabel('Noise Level (σ)');
ylabel('Sampling Rate (× fmax)');
set(gca, 'YDir', 'normal');
colormap('jet');% Connect Fourier Series and Transform
t = -1:0.001:1;       % Time vector
T = 2;                % Period
n = 50;               % Number of harmonics for FFS

% Define a test signal (using Gaussian)
xt_test = exp(-4*t.^2);  % Gaussian centered at 0

% Get Fourier series coefficients
[xhat_test, ck_test] = ffs(xt_test, t, n, T);

% Calculate the Fourier transform
[f_test, xf_test, W_test] = ftr(t, xt_test, T);

% Create frequency points corresponding to Fourier series harmonics
f_series = (-n:n) / T;  % Frequencies for Fourier series: k/T

% Extract CTFT values at those specific frequencies
xf_samples = zeros(size(f_series));
for i = 1:length(f_series)
    % Find closest frequency in the FT result
    [~, idx] = min(abs(f_test - f_series(i)));
    xf_samples(i) = xf_test(idx);
end

% Plot the relationship
figure;
subplot(2,1,1);
stem(-n:n, abs(ck_test), 'b-o'); 
title('Magnitude of Fourier Series Coefficients');
xlabel('Harmonic Index (k)'); ylabel('|c_k|'); grid on;

subplot(2,1,2);
stem(-n:n, abs(xf_samples), 'r-o');
title('Samples of the CTFT at f = k/T');
xlabel('Harmonic Index (k)'); ylabel('|X(k/T)|'); grid on;

% Calculate and display the relationship (should be approximately X(k/T) = T*c_k)
ratio = xf_samples ./ ck_test;
mean_ratio = mean(abs(ratio));
fprintf('Mean ratio |X(k/T)|/|c_k| = %f (should be close to T = %f)\n', mean_ratio, T);

% Display the theoretical relationship
disp('Theoretical relationship: X(k/T) = T·c_k');