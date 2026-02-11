% Test sampling function
addpath('C:\Users\dinaa\OneDrive\Desktop\EECE442_Project\functions');
t = 0:0.001:1;        % Fine time vector for "continuous" signal
f1 = 3;               % First frequency component (Hz)
f2 = 7;               % Second frequency component (Hz)
fmax = max(f1, f2);   % Maximum frequency in the signal

% Create test signal: x(t) = cos(2πf1t) + 0.5cos(2πf2t)
xt = cos(2*pi*f1*t) + 0.5*cos(2*pi*f2*t);

% Try different sampling rates
fs_low = 0.5 * fmax;   
fs_low2= fmax; % Below Nyquist rate (will cause aliasing)
fs_nyquist = fmax * 2;    % At Nyquist rate
     

% Sample the signal at different rates
[t_low, x_low] = sample(t, xt, fs_low);
[t_low2, x_low2] = sample(t, xt, fs_low2);
[t_nyquist, x_nyquist] = sample(t, xt, fs_nyquist);


% Plot original and sampled signals
figure;
subplot(3,1,1);
plot(t, xt, 'k-'); hold on;
plot(t_low, x_low, 'ro', 'MarkerSize', 8);
title(['Sampling at fs = 0.5·fmax = ' num2str(fs_low) ' Hz (Undersampling)']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples');

subplot(3,1,2);
plot(t, xt, 'k-'); hold on;
plot(t_low2, x_low2, 'go', 'MarkerSize', 8);
title(['Sampling at fs = fmax = ' num2str(fs_low2) ' Hz (Also Undersampling)']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples');

subplot(3,1,3);
plot(t, xt, 'k-'); hold on;
plot(t_nyquist, x_nyquist, 'bo', 'MarkerSize', 8);
title(['Sampling at fs = 2*fmax = ' num2str(fs_nyquist) ' Hz (Nyquist Rate)']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples');