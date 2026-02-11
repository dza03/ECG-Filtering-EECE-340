% Connect Fourier Series and Transform
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