function [y, x, inharmHarmEst, harmonics, harmonicAmp, harmonicFreq] = generateAlmostHarmonic(pitch, nHarmonics, N, fs, inharm, SNR)
% Generates a(n almost) harmonic signal.
%
% INPUT
% pitch         -  vector consisting of the fundamental frequencies
% nHarmonics    -  vector consisting of number of harmonics for each
%                  fundamental frequency
% N             -  number of samples
% inharmParam   -  inharmonicity of each frequency
% SNR           -  the SNR of the signal, y
%
% OUTPUT
% y             -  the noisy signal
% x             -  the noise-less signal
% inharmHarmEst -  all frequencies in the signal if inharm = 0, otherwise
%                  the frequencies corresponding to the closest fit to the
%                  harmonic model of the signal
% harmonics     -  all harmonics of the signal
% harmonicAmp   -  all amplitudes of the corresponding true (potentially 
%                  inharmonic) frequencies
% harmonicFreq  -  all frequencies of the true signal


help
% time vector
t = 0:1/fs:1/fs*(N-1);

harmonics = zeros(length(pitch),size(t,2));

for pitchIdx = 1:size(pitch,2)
    for harmonicIdx = 1:nHarmonics(pitchIdx)

        % The amplitude for each frequency, randomized between 0.7 and 1
        harmonicAmp(pitchIdx,harmonicIdx) = (0.7+rand*0.3);

        % The frequency, generated based on the inharmonicity parameter
        harmonicFreq(pitchIdx,harmonicIdx) = pitch(pitchIdx)*harmonicIdx + unifrnd(-pitch(pitchIdx)*harmonicIdx*inharm,pitch(pitchIdx)*harmonicIdx*inharm);

        % Random phase
        harmonicPhase(pitchIdx,harmonicIdx) = 2j*unifrnd(0,pi);

        % The final harmonic signal
        harmonics(pitchIdx,:) = harmonics(pitchIdx,:) + harmonicAmp(pitchIdx,harmonicIdx)*exp(2j*pi*(harmonicFreq(pitchIdx,harmonicIdx)).*t + harmonicPhase(pitchIdx,harmonicIdx));


    end
end
% Adding all harmonic signals together, generating a multi-pitch signal
x = sum(harmonics,1).';
harmonics = harmonics.';

% Estimating the fundamental frequency for inharmonic/harmonic case
for pitchIdx2 = 1:size(pitch,2)
    normRealFreq = harmonicFreq(pitchIdx2,:).'/fs*2*pi;
    realAmp = harmonicAmp(pitchIdx2,:).';
    realPhase = harmonicPhase(pitchIdx2,:).'/(1j);
    sample_t = t.'*fs;

    sigma2_tilde = 10^(-SNR/10)*sum(harmonicAmp(pitchIdx2,:).^2);

    if abs(inharm) ~= 0
        harmonicVec = [1:nHarmonics(pitchIdx2)].';
        normSearchInterval = [pitch(pitchIdx2)-pitch(pitchIdx2)*inharm*3, pitch(pitchIdx2)+pitch(pitchIdx2)*inharm*3]/fs*2*pi;
        nGridSearch = 500;
        [omega0,r_pseudo,phi_pseudo] = pseudotrue_omega0(normRealFreq,realAmp,realPhase,sample_t,harmonicVec,normSearchInterval,nGridSearch);
        mu = exp(2i*pi*t.'*(omega0*fs/2/pi*harmonicVec)')*(r_pseudo.*exp(1i*phi_pseudo));
        inharmHarmEst(pitchIdx2,1:length(harmonicVec)) = omega0*fs/2/pi*harmonicVec;
    else
        harmonicVec = zeros(1, max(nHarmonics));
        harmonicVec(1:nHarmonics(pitchIdx2)) = [1:nHarmonics(pitchIdx2)].';
        inharmHarmEst(pitchIdx2,1:length(harmonicVec)) = normRealFreq(1)*fs/2/pi*harmonicVec;
    end

end


% Power of the signal
power = sum(harmonicAmp.^2,'all');
% Power of the noise
p_noise = power*10^(-SNR/10);
% Signal with noise
y = x + normrnd(0,sqrt(p_noise/2), size(x)) + 1j*normrnd(0,sqrt(p_noise/2), size(x));

end