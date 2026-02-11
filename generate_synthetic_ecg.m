% create_synthetic_ecg.m
% Script to create synthetic ECG data for normal and tachycardia conditions


% Parameters
duration = 10;     % 10 seconds of data
fs = 360;          % 360 Hz sampling rate (common for ECG)

% Create normal ECG data (60-100 bpm)
fprintf('Generating normal ECG at 72 bpm...\n');
[normal_ecg, t_normal] = generate_ecg(72, duration, fs);  % 72 bpm

% Create tachycardia ECG data (>100 bpm)
fprintf('Generating tachycardia ECG at 140 bpm...\n');
[tachy_ecg, t_tachy] = generate_ecg(140, duration, fs);   % 140 bpm

% Save to MAT file
fprintf('Saving ECG data to ecg_data.mat...\n');
save('ecg_data.mat', 'normal_ecg', 't_normal', 'tachy_ecg', 't_tachy');

% Plot for comparison
figure;
subplot(2,1,1);
plot(t_normal, normal_ecg);
title('Normal ECG (72 bpm)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t_tachy, tachy_ecg);
title('Tachycardia ECG (140 bpm)');
xlabel('Time (s)'); ylabel('Amplitude');
grid on;

