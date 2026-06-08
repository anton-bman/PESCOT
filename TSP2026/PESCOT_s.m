function [f0s, x_k, pitchGrid, pitchDist, M] = PESCOT_s(y, epsilon, zeta, eta, beta, varargin)
%%%%
%   PESCOT - Multi-Pitch Estimation with Sparse Clustering through Optimal
%   Transport
%
%   [f0s, alpha, M] = PESCOT_s(y, epsilon, zeta, eta, beta, nPitches,...
%   max_iter, pitchLim, maxHarm, nPitchGrid, nFreqsGrid, doPrint)
%
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
%   freqGrid    frequency grid (default uniform grid from 50/8000 to 1/2, 
%               2260 grid points)
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
% Reference: "Inverse harmonic clustering for multi-pitch estimation:
% an optimal transport approach", TSP 2026.
%
% Implemented by: Anton Björkman
% Date: June 8, 2026

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
    % if maxHarm ~= Inf
    %     nFreqsGrid = nPitchGrid*maxHarm;
    % else
    %     nFreqsGrid = nPitchGrid*10;
    % end
end

if nargin >= 11 && ~isempty(varargin{6})
    freqGrid = varargin{6};
    % nFreqsGrid = varargin{6};
else
    nFreqsGrid = nPitchGrid*10;
    if pitchLim(2)*maxHarm > 1/2
        freqGrid = linspace(pitchLim(1),1/2,nFreqsGrid);
    else
        freqGrid = linspace(pitchLim(1),pitchLim(2)*maxHarm,nFreqsGrid);
    end
end
if nargin >= 12 && ~isempty(varargin{7})
    doPrint = varargin{7};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters that can be changed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter_prox_operator = 1000;
nesterov = 1; % "nesterov acceleration", speeds up convergence
tol = 4e-3; % convergence tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dictionary and cost matrix generation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pitchGrid = linspace(pitchLim(1),pitchLim(2),nPitchGrid);

Ns = length(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% new stuff %%%%%%%%%%%%%%%%%%%
nDelays = round(Ns*2/3);
delays = 1:nDelays;
R_y = xcorr(y, 'unbiased');
r2_all = (R_y(Ns+1:Ns+nDelays));
rfull = [real(r2_all); imag(r2_all)];
A_full = exp(2j*pi*delays(:)*freqGrid);
Afull = [real(A_full); imag(A_full)];

[U,S,V] = svd(Afull,'econ');
rank_A = rank(Afull);
U = U(:, 1:rank_A);
S = S(1:rank_A, 1:rank_A);
V = V(:, 1:rank_A);

r = U'*rfull/std(rfull);
A = S*V';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = generateCostMatrix(freqGrid, pitchGrid, maxHarm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if doPrint > 0
    figure(doPrint)
end

N = length(r);
L = (norm(A)^2)/N;

gamma = 1/(0.2*L);
x_k = abs((A'/N)*r);
theta = 1;
stepsize_gamma = 2;
z_k = x_k;
z_knew = x_k;


dual0 = zeros(length(x_k),1);
dual1 = zeros(size(C));

tk = 1;
for iter = 1:max_iter
    % "improved Nestrov acceleration"
    if nesterov == 1
        yk = (1-theta)*x_k + theta*z_knew;
        gradient = -A'*(r-A*yk)/N + beta;
        gamma = theta^(1-stepsize_gamma)/L;
    else
        yk = z_k;
        gradient = -A'*(r-A*yk)/N + beta;
        gamma = 1;
    end
    
    [z_knew,dual0,dual1] = prox_op_spectrum(C,z_k,gradient, epsilon, zeta, eta, gamma,max_iter_prox_operator,dual0,dual1);

    if nesterov == 1
        x_k = (1-theta)*x_k + theta*z_knew;
        theta = find_theta_next(theta, stepsize_gamma);
    else
        x_k = z_knew;
    end

    rel_error(iter) = norm(z_knew-z_k)/norm(z_knew);

    z_k = z_knew;








    if doPrint > 0
        logM = -C/epsilon + dual1./(zeta*gamma*epsilon) + dual0/(zeta*gamma*epsilon)*ones(1,size(dual1,2));
        M = exp(logM);
        fprintf("Iteration %.d \n ",iter)
        subplot(1,2,1)
        stem(pitchGrid,sum(squeeze(M),1))
        xlabel('Pitch (cycles/sample)')

        subplot(1,2,2)
        stem(freqGrid,sum(squeeze(M),2))
        xlabel('Frequency (cycles/sample)')
        drawnow;
    end

    if iter > 1
        if rel_error(iter) < tol || isnan(rel_error(iter))
            break;
        end
    end



    if sum(isnan(x_k)) > 0
        fprintf("Numerical instability \n")
        return
    end


end

logM = -C/epsilon + dual1./(zeta*gamma*epsilon) + dual0/(zeta*gamma*epsilon)*ones(1,size(dual1,2));
M = exp(logM);

M = M*stdy.^2;
x_k = x_k*stdy;

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
end
f0s = (pitchGrid_tot).';
end


function theta_k_plus_1 = find_theta_next(theta_k, gamma)
    % inequality function
    ineq_func = @(theta_next) (1-theta_next)./(theta_next.^gamma)-1./(theta_k^gamma);
    
    % find theta_k_plus_1 using linspace
    theta_candidates = linspace(1e-3, 1-1e-3, 10000);  % Avoid boundaries (0 and 1)
    values = ineq_func(theta_candidates);
    
    % find where the sign changes
    idx = find(values(1:end-1).*values(2:end)<0,1);  % Look for sign change
    
    if isempty(idx)
        error('No valid theta found in (0, 1)');
    end

    % interpolation for a better estimate
    theta1 = theta_candidates(idx);
    theta2 = theta_candidates(idx + 1);
    val1 = values(idx);
    val2 = values(idx + 1);
    
    % estimate root using interpolation
    theta_k_plus_1 = theta1 - val1 * (theta2 - theta1) / (val2 - val1);
end