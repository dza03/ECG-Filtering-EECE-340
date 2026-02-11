% Parameters
addpath('C:\Users\dinaa\OneDrive\Desktop\EECE 340_PROJ\functions');
T = 2;        % Period (since both signals are defined on -1 to 1)
n = 1000;       % Number of harmonics
t = -1:0.001:1;  % Time vector

xt_sin = sin(2*pi*t) .* (abs(t) < 1);  % Sine wave from -1 to 1

[xhat_sin, ck_sin] = ffs(xt_sin, t, n, T);

figure;
plot(t, xt_sin, 'k', 'LineWidth', 1.5); hold on;
plot(t, xhat_sin, 'r--');
legend('Original Sine', 'Fourier Approximation');
title(['Sine: Finite Fourier Series Approximation with n=', num2str(n)])
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

E_sin = trapz(t, abs(xt_sin - xhat_sin).^2);  % Error for sine
disp(['Sine Signal Error E = ', num2str(E_sin)]);


xt_gauss = exp(-t.^2);  % Gaussian pulse (time-limited naturally due to decay)

[xhat_gauss, ck_gauss] = ffs(xt_gauss, t, n, T);

figure;
plot(t, xt_gauss, 'k', 'LineWidth', 1.5); hold on;
plot(t, xhat_gauss, 'r--');
legend('Original Gaussian', 'Fourier Approximation');
title(['Gaussian: Finite Fourier Series Approximation with n=', num2str(n)])
xlabel('Time (s)'); ylabel('Amplitude'); grid on;

E_gauss = trapz(t, abs(xt_gauss - xhat_gauss).^2);  % Error for Gaussian
disp(['Gaussian Signal Error E = ', num2str(E_gauss)]);
