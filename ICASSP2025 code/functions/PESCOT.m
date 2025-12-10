function [f0s, alpha, pitchGrid, pitchDist, M] = PESCOT(y, epsilon, zeta, eta, beta, varargin)
%%%%
%   PESCOT - Multi-Pitch Estimation with Sparse Clustering through Optimal
%   Transport
%
%   [f0s, alpha, M] = PESCOT(y, epsilon, zeta, eta, beta, nPitches,...
%   max_iter, pitchLim, maxHarm, nPitchGrid, nFreqsGrid, doPrint)
%
%   Algorithm solving the problem min_alpha 1/N||y-A*alpha||_2^2 +
%   beta||alpha||_1 + zeta S(alpha), where S(alpha) is a modified optimal
%   transport problem.
%
%   INPUTS: 
%   y           multipitch input signal (or single-pitch)
%   epsilon     entropic regularization parameter, epsilon > 0 
%   zeta        sparsity regularization parameter, zeta > 0 
%   eta         OT problem regularization parameter, eta > 0 
%   beta        l1 norm regularization parameter, beta > 0
%
%   OPTIONAL INPUTS:
%
%   nPitches    known or assumed number of active pitches in the signal
%   max_iter    maximum number of iterations (default: 1e3) 
%   pitchLim    minimum and maximum evaluated normalized pitch
%               (default: [50, 500]/8000)
%   maxHarm     highest assumed harmonic order, H (default 10)
%   nPitchGrid  number of grid points for the pitch grid, G (default 226)
%   nFreqsGrid  number of grid points for the frequency grid, F (default
%               2260)
%   doPrint     integer value - if doPrint > 0, a figure(doPrint) will 
%               show the updates of the pitch and frequency distributions
%               (default 0)
%   
%
%   OUTPUTS: 
%   f0s         Vector consisting of the (normalized) fundamental
%               frequencies, in order of pitch with most assigned power.
%               Only returns the pitches with most power that collectively
%               sum up to at least 95% of the total power in the pitch
%               distribution, if nPitches is not specified. If there are
%               fewer active pitches than nPitches, the additional pitches
%               will be identified as 0.
%   alpha       complex valued amplitudes that solve the problem 
%   pitchGrid   pitch grid 
%   pitchDist   pitch distribution 
%   M           transport plan from the last iteration of the proximal
%               operator
%
%   EXAMPLE 1: 
%   f0s = PESCOT(y, epsilon, zeta, eta, beta)
%
%   EXAMPLE 2:
%   [f0s, alpha, pitchGrid, pitchDist, M]] = PESCOT(y, epsilon, zeta, ...
%   eta, beta, nPitches, max_iter, pitchLim, maxHarmonics, nPitchGrid, ...
%   nFreqsGrid, doPrint)
%
%
% Reference: "Robust Multi-Pitch Estimation via Optimal Transport
% Clustering", submitted to ICASSP 2025.
%
% Implemented by: Anton Björkman
% Date: October 11, 2024

if isrow(y)
    y = y.';
end

stdy = std(y);
y = y/stdy;

% default values
max_iter = 1000;
pitchLim = [50, 500]/8000;
nPitchGrid = 226;
maxHarm = 10;
nFreqsGrid = nPitchGrid*10;
doPrint = 0;
nPitches = [];

if nargin >= 6 && ~isempty(varargin{1})
    nPitches = varargin{1};
end
if nargin >= 7 && ~isempty(varargin{2})
    max_iter = varargin{2};
end
if nargin >= 8 && ~isempty(varargin{3})
    pitchLim = varargin{3};
end
if nargin >= 9 && ~isempty(varargin{4})
    maxHarm = varargin{4};
end
if nargin >= 10 && ~isempty(varargin{5})
    nPitchGrid = varargin{5};
    if maxHarm ~= Inf
        nFreqsGrid = nPitchGrid*maxHarm;
    else
        nFreqsGrid = nPitchGrid*10;
    end
end
if nargin >= 11 && ~isempty(varargin{6})
    nFreqsGrid = varargin{6};
end
if nargin >= 12 && ~isempty(varargin{7})
    doPrint = varargin{7};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that can be changed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter_prox_operator = 1000;
nesterov = 1; % "nesterov acceleration", speeds up convergence
tol = 1e-3; % convergence tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dictionary and cost matrix generation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pitchGrid = linspace(pitchLim(1),pitchLim(2),nPitchGrid);
if pitchLim(2)*maxHarm > 1/2
    freqGrid = linspace(pitchLim(1),1/2,nFreqsGrid);
else
    freqGrid = linspace(pitchLim(1),pitchLim(2)*maxHarm,nFreqsGrid);
end
t = 0:1:1*(length(y)-1);
A = exp(2j*pi*t(:)*freqGrid);
C = generateCostMatrix(freqGrid, pitchGrid, maxHarm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if doPrint > 0
    figure(doPrint)
end

N = length(y);
L = (norm(A)^2)/N;

gamma = 1/L;
alpha = A'*y/N;
alpha_old = alpha;

dual0 = zeros(length(alpha),1);
dual1 = zeros(size(C));

tk = 1;
for iter = 1:max_iter
    % "improved Nestrov acceleration"
    if nesterov == 1
        tk_new = (1+sqrt(1+4*tk^2))/2;
        v = alpha + (tk-1)/tk_new*(alpha-alpha_old);
        tk = tk_new;
    else
        v = alpha;
    end

    % Gradient of LS term
    gradient = -A'*(y-A*v)/N;
    [alpha_new,dual0,dual1] = proximal_operator_square_l1(C,v-gamma*gradient,epsilon, zeta, eta, gamma,max_iter_prox_operator,dual0,dual1, beta);


    if sum(isnan(alpha_new)) > 0
        alpha = alpha_old;
        fprintf("Numerical instability \n")
        return
    end



    if doPrint
        logM = -C/epsilon + dual1./(zeta*gamma*epsilon) + dual0/(zeta*gamma*epsilon)*ones(1,size(dual1,2));
        M = exp(logM);
        fprintf("Iteration %.d \n ",iter)
        subplot(1,2,1)
        stem(pitchGrid*8000, sum(squeeze(M),1))
        % ylim([0 max(0.05/stdy^2/3, max(sum(squeeze(M),1)))])
        xlabel('Pitch (Hz)')

        subplot(1,2,2)
        stem(freqGrid*8000, sum(squeeze(M),2))
        xlabel('Frequency (Hz)')

        drawnow;
    end

    if iter > 1
        if norm(alpha_new-alpha)/norm(alpha_new) < tol
            break;
        end
    end

    alpha_old = alpha;
    alpha = alpha_new;


end

logM = -C/epsilon + dual1./(zeta*gamma*epsilon) + dual0/(zeta*gamma*epsilon)*ones(1,size(dual1,2));
M = exp(logM);

M = M*stdy.^2;
alpha = alpha*stdy;

pitchDist = sum(squeeze(M),1);
[pitchamp_sorted, idx] = sort(pitchDist, 'descend');
pitchGrid_sorted = pitchGrid(idx);


if isempty(nPitches)
    cumuSum = cumsum(pitchamp_sorted);
    totalSum = sum(pitchDist);
    thresh = 0.95 * totalSum;
    idx95 = find(cumuSum >= thresh, 1);
    pitchGrid_tot = pitchGrid_sorted(1:idx95);
else
    pitchGrid_tot = pitchGrid_sorted(1:nPitches);
    pitchGrid_tot = pitchGrid_tot.*(~pitchamp_sorted(1:nPitches) == 0);
    % pitchGrid_tot = pitchGrid_tot(pitchamp_sorted(1:nPitches) > 0);
end
f0s = (pitchGrid_tot).';
end