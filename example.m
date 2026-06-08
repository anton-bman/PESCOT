close all
clear all
clc

%% Generating the signal

fs = 8000; %Sampling frequency
pitches = [197, 240, 272];
nActivePitches = length(pitches);
nHarmonics = randi([3 10], 1, nActivePitches);
N = 250; % number of samples
inharm = 0.01; % inharmonicity of the signal
SNR = 10; % SNR of the signal

[y, ~, inharmHarmEst] = generateAlmostHarmonic(pitches, nHarmonics, N, fs, inharm, SNR);
realPitch = inharmHarmEst(:,1);


%% PESCOT algorithm, finding the pitches



%%%%% Parameter values
epsilon = 1e-7; % entropic regularization parameter
eta = 1e0; % OT sparsity regularization parameter
zeta = 8e3; % OT regularization parameter
beta = 8e-2; % l1 regularization parameter

%%%%% Parameter values, stochastic model
eta_s = 1e-1; % OT sparsity regularization parameter
zeta_s = 1e1; % OT regularization parameter
epsilon_s = 1e-6; % entropic regularization parameter
beta_s = 1.5e-2;

%%%%% Parameter values, deterministic model
eta_1 = 1e-2; % OT sparsity regularization parameter
zeta_1 = 1e2; % OT regularization parameter
epsilon_1 = 1e-6; % entropic regularization parameter
beta_1 = 8e-2;

%%%%% Optional inputs
nPitches = nActivePitches; % Assumed number of active pitches
max_iter = 1000; % Maximum iterations, default 1000
pitchLim = [50, 500]/fs; % Minimum and maximum evaluated pitch
maxHarmonics = 10; % Assumed maximum number of harmonics
% Number of grid points for the pitch grid
nPitchGrid = 226; 
freqGrid = linspace(50,4000,nPitchGrid*maxHarmonics)/fs;
doPrint = 1; % Plotting the OT estimate of pitch distribution.
% Note: having doPrint set to integer > 0 slows everything down.



%% ICASSP 2025 estimator


%%% Running the actual algorithm
if rand >= 0.5 % :)
    % Using only the necessary inputs (others are set to default values)
    f0s = PESCOT(y, epsilon, zeta, eta, beta);
else
    % Using all optional inputs
    [f0s, alpha, pitchGrid, pitchDist, M] = PESCOT(y, epsilon, zeta, eta, beta, nPitches, ...
        max_iter, pitchLim, maxHarmonics, nPitchGrid, freqGrid, doPrint);
end


%% Stochastic model
f0s_s = PESCOT_s(y, epsilon_s, zeta_s, eta_s, beta_s, nPitches, ...
        100, pitchLim, maxHarmonics, nPitchGrid, freqGrid, doPrint);

%% Deterministic model

f0s_d = PESCOT_1(y, epsilon_1, zeta_1, eta_1, beta_1, nPitches, ...
        max_iter, pitchLim, maxHarmonics, nPitchGrid, freqGrid, doPrint);


%% Presenting the results

sort_f0s = sort(f0s(1:min(nPitches, nActivePitches)), 'descend')*fs;
sort_f0s_s = sort(f0s_s(1:min(nPitches, nActivePitches)), 'descend')*fs;
sort_f0s_d = sort(f0s_d(1:min(nPitches, nActivePitches)), 'descend')*fs;

sort_pitch = sort(realPitch, 'descend');

fprintf('\n')
fprintf(' Most probable pitch candidates \n');
fprintf('  real pitch    PESCOT-2    PESCOT-1    PESCOT-s \n');
fprintf('  ----------------------------------------------\n');
for i = 1:nActivePitches
    pitch_value = '';
    f0_value = '';
    if i <= length(sort_pitch)
        pitch_value = sort_pitch(i);
    end
    if i <= length(sort_f0s)
        f0_value = sort_f0s(i);
    end
    if i <= length(sort_f0s_d)
        f0_d_value = sort_f0s_d(i);
    end
    if i <= length(sort_f0s_s)
        f0_s_value = sort_f0s_s(i);
    end
    fprintf('  %8.3f %12.3f %11.3f  %10.3f\n', pitch_value, f0_value, f0_s_value, f0_d_value);
end
fprintf('  ----------------------------------------------\n');
fprintf('             Grid spacing = %f         \n', ((pitchLim(2)-pitchLim(1))/(nPitchGrid-1))*fs);
