% Test reconstruction function
% Reconstruct signals from their samples
% Time vector for reconstruction (finer than original for better visualization)
t_recon = 0:0.0005:1;

% Reconstruct from different sampling rates
xrcon_low = reconstruct(t_recon, x_low, t_low, fs_low);
xrcon_low2 = reconstruct(t_recon, x_low2, t_low2, fs_low2);
xrcon_nyquist = reconstruct(t_recon, x_nyquist, t_nyquist, fs_nyquist);


% Plot original and reconstructed signals
figure;
subplot(3,1,1);
plot(t, xt, 'k-'); hold on;
plot(t_low, x_low, 'ro', 'MarkerSize', 6);
plot(t_recon, xrcon_low, 'r--', 'LineWidth', 1.5);
title('Reconstruction from Undersampled Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples', 'Reconstructed');

subplot(3,1,3);
plot(t, xt, 'k-'); hold on;
plot(t_low2, x_low2, 'bo', 'MarkerSize', 6);
plot(t_recon, xrcon_low2, 'b--', 'LineWidth', 1.5);
title('Reconstruction from Second Undersampled Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples', 'Reconstructed');

subplot(3,1,2);
plot(t, xt, 'k-'); hold on;
plot(t_nyquist, x_nyquist, 'go', 'MarkerSize', 6);
plot(t_recon, xrcon_nyquist, 'g--', 'LineWidth', 1.5);
title('Reconstruction from Nyquist-Sampled Signal');
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples', 'Reconstructed');


% Calculate reconstruction errors
error_low = mean(abs(interp1(t, xt, t_recon) - xrcon_low).^2);
error_low2 = mean(abs(interp1(t, xt, t_recon) - xrcon_low2).^2);
error_nyquist = mean(abs(interp1(t, xt, t_recon) - xrcon_nyquist).^2);

fprintf('Mean square error (Undersampling): %e\n', error_low);
fprintf('Mean square error (Also Undersampling): %e\n', error_low2);
fprintf('Mean square error (Nyquist rate): %e\n', error_nyquist);
