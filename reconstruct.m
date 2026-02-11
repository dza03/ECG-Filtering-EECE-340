function xrcon = reconstruct(t, x_sample, t_sample, fs)
    % Inputs:
    % t - time vector for reconstruction points
    % x_sample - sampled signal values
    % t_sample - time vector of sampled points
    % fs - sampling frequency in Hz
    
    % Outputs:
    % xrcon - reconstructed signal at points specified by t
    
    % Sampling period
    Ts = 1/fs;
    
    % Initialize reconstructed signal
    xrcon = zeros(size(t));
    
    % Sinc reconstruction (Whittaker–Shannon interpolation)
    for i = 1:length(t)
        % Sum of weighted sinc functions
        for k = 1:length(t_sample)
            % Normalized sinc function
            sinc_val = sinc((t(i) - t_sample(k))/Ts);
            xrcon(i) = xrcon(i) + x_sample(k) * sinc_val;
        end
    end
end