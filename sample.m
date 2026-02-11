function [t_sample, x_sample] = sample(t, xt, fs)
    % Inputs:
    % t - time vector of original signal
    % xt - original continuous signal
    % fs - sampling frequency in Hz
    
    % Outputs:
    % t_sample - time vector of sampled points
    % x_sample - sampled signal values
    
    % Calculate sampling period
    Ts = 1/fs;
    
    % Find min and max time values
    t_min = min(t);
    t_max = max(t);
    
    % Create time vector for samples
    t_sample = t_min:Ts:t_max;
    
    % Initialize sampled signal
    x_sample = zeros(size(t_sample));
    
    % Interpolate the original signal at the sampling instants
    x_sample = interp1(t, xt, t_sample, 'linear');
end