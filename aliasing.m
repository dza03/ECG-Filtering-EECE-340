addpath('C:\Users\dinaa\OneDrive\Desktop\EECE 340_PROJ\functions')
t_demo = 0:0.001:1;           % Time vector
f_orig = 5;                   % Original frequency (Hz)
fs_demo = 8;                  % Sampling frequency (Hz)

% Generate the original 5 Hz cosine
x_orig = cos(2*pi*f_orig*t_demo);

% Sample the signal at 8 Hz
[t_sample_demo, x_sample_demo] = sample(t_demo, x_orig, fs_demo);
%cosines that could produce the same samples
f_alias1 = fs_demo - f_orig;  % First alias frequency (3 Hz)
f_alias2 = fs_demo + f_orig;  % Second alias frequency (13 Hz)
f_alias3 = 2*fs_demo - f_orig;  % Third alias frequency (11 Hz)

x_alias1 = cos(2*pi*f_alias1*t_demo);
x_alias2 = cos(2*pi*f_alias2*t_demo);
x_alias3 = cos(2*pi*f_alias3*t_demo);

% Plot to demonstrate aliasing
figure;
subplot(4,1,1);
plot(t_demo, x_orig, 'b-'); hold on;
plot(t_sample_demo, x_sample_demo, 'ro', 'MarkerSize', 8);
title(['Original Signal: cos(2\pi * ' num2str(f_orig) 't) - Sampled at ' num2str(fs_demo) ' Hz']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Original', 'Samples');

subplot(4,1,2);
plot(t_demo, x_alias1, 'g-'); hold on;
plot(t_sample_demo, x_sample_demo, 'ro', 'MarkerSize', 8);
title(['First Alias: cos(2\pi * ' num2str(f_alias1) 't) - Same Samples!']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Alias Signal', 'Samples');

subplot(4,1,3);
plot(t_demo, x_alias2, 'm-'); hold on;
plot(t_sample_demo, x_sample_demo, 'ro', 'MarkerSize', 8);
title(['Second Alias: cos(2\pi * ' num2str(f_alias2) 't) - Same Samples!']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Alias Signal', 'Samples');

subplot(4,1,4);
plot(t_demo, x_alias3, 'c-'); hold on;
plot(t_sample_demo, x_sample_demo, 'ro', 'MarkerSize', 8);
title(['Third Alias: cos(2\pi * ' num2str(f_alias3) 't) - Same Samples!']);
xlabel('Time (s)'); ylabel('Amplitude'); grid on;
legend('Alias Signal', 'Samples');{}