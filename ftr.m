% Numerical Continuous-Time Fourier Transform
function [f, xf, W] = ftr(t, xt, T)
    % Inputs:
    % t - time vector
    % xt - time domain signal
    % T - period to assume
    
    % Outputs:
    % f - frequency vector
    % xf - frequency domain signal (CTFT)
    % W - frequency range
    
    dt = t(2) - t(1);           % Time step
    N = length(t);              % Number of time samples
    
    % Calculate frequency range (W) and frequency vector (f)
    df = 1/T;                   % Frequency resolution
    W = 1/dt;                   % Total frequency range (related to sampling rate)
    f = (-N/2:N/2-1) * df;      % Centered frequency vector
    
    % Compute CTFT numerically using integration
    xf = zeros(size(f));
    for i = 1:length(f)
        xf(i) = sum(xt .* exp(-1j * 2 * pi * f(i) * t)) * dt;
    end
end