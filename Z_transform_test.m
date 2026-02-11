% Z-transform of FIR filter
syms z;
H_z = 0;
for k = 1:length(h)
    H_z = H_z + h(k) * z^(-(k-1));
end
disp('System function H(z):');
pretty(H_z)

% Stability: FIR => always stable (finite sum)
% Causality: Check if h(n) = 0 for n < 0 (which is true if you center properly)
