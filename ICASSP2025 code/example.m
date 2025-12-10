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



%%%%% Optional inputs
nPitches = nActivePitches; % assumed number of active pitches
max_iter = 1000; % maximum iterations, default 1000
pitchLim = [50, 500]/fs; %Minimum and maximum evaluated pitch
maxHarmonics = 10; %Assumed maximum number of harmonics
% Number of grid points for pitch and frequency grids
nPitchGrid = 226; 
nFreqsGrid = maxHarmonics*nPitchGrid;
doPrint = 1; % printing the OT estimate of pitch distribution.
% Note: having doPrint set to integer > 0 slows everything down.


%%% Running the actual algorithm
if rand >= 0.5 % :)
    % Using only the necessary inputs (others are set to default values)
    f0s = PESCOT(y, epsilon, zeta, eta, beta);
else
    % Using all optional inputs
    [f0s, alpha, pitchGrid, pitchDist, M] = PESCOT(y, epsilon, zeta, eta, beta, nPitches, ...
        max_iter, pitchLim, maxHarmonics, nPitchGrid, nFreqsGrid, doPrint);
end




%% Presenting the results

sort_f0s = sort(f0s(1:min(nPitches, nActivePitches)), 'descend')*fs;
sort_pitch = sort(realPitch, 'descend');

fprintf('\n')
fprintf(' Most probable pitch candidates \n');
fprintf('  real pitch         est pitch \n');
fprintf('  ----------------------------\n');
for i = 1:nActivePitches
    pitch_value = '';
    f0_value = '';
    if i <= length(sort_pitch)
        pitch_value = sort_pitch(i);
    end
    if i <= length(sort_f0s)
        f0_value = sort_f0s(i);
    end
    fprintf('  %8.5f  %17.5f\n', pitch_value, f0_value);
end
fprintf('  ----------------------------\n');
fprintf('     Grid spacing = %f         \n', ((pitchLim(2)-pitchLim(1))/(nPitchGrid-1))*fs);
