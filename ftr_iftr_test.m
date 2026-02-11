% Test Fourier Transform and Inverse
addpath('C:\Users\dinaa\OneDrive\Desktop\EECE 340_PROJ\functions');
t = -5:0.01:5;             % Time vector
T = 10;                    % Period to assume

% Test signal 1: Gaussian pulse
xt_gauss = exp(-t.^2);

% Apply CTFT
[f_gauss, xf_gauss, W_gauss] = ftr(t, xt_gauss, T);

% Apply inverse CTFT
[t_rec, xt_rec_gauss, T_rec] = iftr(f_gauss, xf_gauss, W_gauss);

% Plot original and reconstructed signals
figure;
subplot(2,2,1);
plot(t, xt_gauss); title('Original Gaussian'); xlabel('Time (s)'); grid on;

subplot(2,2,2);
plot(f_gauss, abs(xf_gauss)); title('Magnitude of CTFT'); xlabel('Frequency (Hz)'); grid on;

subplot(2,2,3);
plot(t_rec, xt_rec_gauss); title('Reconstructed Gaussian'); xlabel('Time (s)'); grid on;

subplot(2,2,4);
plot(t, xt_gauss - xt_rec_gauss); title('Reconstruction Error'); xlabel('Time (s)'); grid on;

% Test signal 2: Sinc function (known FT pair with rectangular pulse)
xt_sinc = sinc(2*t);     % sinc(x) = sin(πx)/(πx)

% Apply CTFT
[f_sinc, xf_sinc, W_sinc] = ftr(t, xt_sinc, T);

% Apply inverse CTFT
[t_rec, xt_rec_sinc, T_rec] = iftr(f_sinc, xf_sinc, W_sinc);

% Plot original and reconstructed signals
figure;
subplot(2,2,1);
plot(t, xt_sinc); title('Original Sinc'); xlabel('Time (s)'); grid on;

subplot(2,2,2);
plot(f_sinc, abs(xf_sinc)); title('Magnitude of CTFT'); xlabel('Frequency (Hz)'); grid on;

subplot(2,2,3);
plot(t_rec, xt_rec_sinc); title('Reconstructed Sinc'); xlabel('Time (s)'); grid on;

subplot(2,2,4);
plot(t, xt_sinc - xt_rec_sinc(1:length(t))); title('Reconstruction Error'); xlabel('Time (s)'); grid on;