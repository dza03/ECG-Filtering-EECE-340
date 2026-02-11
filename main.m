addpath('functions');
% 1. Sampling & Reconstruction
run('tests/sample_test.m');
run('tests/reconstruction_test.m');

% 2. Aliasing & Noise
run('tests/aliasing.m');
run('tests/noise_robustness.m');

% 3. Fourier Analysis
run('tests/ftr_iftr_test.m');
run ('tests/connect.m');
run('tests/Z_transform_test.m');
run('tests/ffs_test.m');
run('tests/Filter_Test.m');
