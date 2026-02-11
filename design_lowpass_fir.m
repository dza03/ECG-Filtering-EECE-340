function h = design_lowpass_fir(cutoff_freq, fs, filter_order)
    % Inputs:
    % cutoff_freq - cutoff frequency in Hz
    % fs - sampling frequency in Hz
    % filter_order - number of taps (should be odd)
    %
    % Output:
    % h - filter impulse response (FIR)

    % Normalized cutoff (0 to 1, where 1 = Nyquist freq)
    wc = cutoff_freq / (fs/2);

    % Ensure filter_order is odd
    if mod(filter_order,2) == 0
        error('Filter order must be odd.');
    end

    % Sinc filter (ideal LPF) * window (Hamming)
    n = -(filter_order-1)/2:(filter_order-1)/2;
    h = wc * sinc(wc * n);  % ideal LPF

    % Apply Hamming window (To smooth out sinc sharp edges (Idea suggested by ChatGPT))
    h = h .* hamming(filter_order)';

    % Normalize for unity gain at DC
    h = h / sum(h);
end
