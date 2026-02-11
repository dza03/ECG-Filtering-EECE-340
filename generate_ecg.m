%This code is taken from publishe MATLAB libraries and GitHub (Referenced in report)
function [ecg_signal, t] = generate_ecg(heart_rate, duration, fs)
    % GENERATE_ECG Creates synthetic ECG with specified heart rate
    % Inputs:
    %   heart_rate - Heart rate in beats per minute (bpm)
    %   duration - Duration of signal in seconds
    %   fs - Sampling frequency in Hz
    %
    % Outputs:
    %   ecg_signal - Synthetic ECG signal
    %   t - Time vector

    % Create time vector
    t = 0:1/fs:duration;
    
    % Set default parameters
    li = 30/heart_rate;  % Adjusts RR interval based on heart rate
    
    % P wave parameters
    a_pwav = 0.25;
    d_pwav = 0.09;
    t_pwav = 0.16;
    
    % Q wave parameters
    a_qwav = 0.025;
    d_qwav = 0.066;
    t_qwav = 0.166;
    
    % QRS complex parameters
    a_qrswav = 1.6;
    d_qrswav = 0.11;
    
    % S wave parameters
    a_swav = 0.25;
    d_swav = 0.066;
    t_swav = 0.09;
    
    % T wave parameters
    a_twav = 0.35;
    d_twav = 0.142;
    t_twav = 0.2;
    
    % U wave parameters
    a_uwav = 0.035;
    d_uwav = 0.0476;
    t_uwav = 0.433;
    
    % Generate individual waves
    pwav = p_wav(t, a_pwav, d_pwav, t_pwav, li);
    qwav = q_wav(t, a_qwav, d_qwav, t_qwav, li);
    qrswav = qrs_wav(t, a_qrswav, d_qrswav, li);
    swav = s_wav(t, a_swav, d_swav, t_swav, li);
    twav = t_wav(t, a_twav, d_twav, t_twav, li);
    uwav = u_wav(t, a_uwav, d_uwav, t_uwav, li);
    
    % Combine all waves to form ECG signal
    ecg_signal = pwav + qwav + qrswav + swav + twav + uwav;
end