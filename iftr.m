% Numerical Inverse Continuous-Time Fourier Transform
function [t, xt_rec, T] = iftr(f, xf, W)
    % Inputs:
    % f - frequency vector
    % xf - frequency domain signal
    % W - frequency range
    
    % Outputs:
    % t - reconstructed time vector
    % xt_rec - reconstructed time domain signal
    % T - period
    
    df = f(2) - f(1);           % Frequency step
    N = length(f);              % Number of frequency samples
    
    % Calculate time range and time vector
    T = 1/df;                   % Period (time range)
    dt = 1/W;                   % Time resolution
    t = (-N/2:N/2-1) * dt;      % Centered time vector
    
    % Compute inverse CTFT numerically
    xt_rec = zeros(size(t));
    for i = 1:length(t)
        xt_rec(i) = sum(xf .* exp(1j * 2 * pi * f * t(i))) * df;
    end
    
    % Take real part assuming real signal
    xt_rec = real(xt_rec);
end