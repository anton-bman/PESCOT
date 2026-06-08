function C = generateCostMatrix(freqGrid, pitchGrid, varargin)
%%%%
%   generateCostMatrix: generates a cost matrix based on frequencies,
%   pitch candidates and maximum number of harmonics per pitch
%
%   cost = generateCostMatrix(freqs, pitch_cand, maxHarmonics)
%
%
%   INPUTS: 
%   freqGrid    frequencies considered 
%   pitchGrid   pitches considered
%
%   OPTIONAL INPUTS:
%   maxHarm     maximum harmonic order (default: 10)
%   
%
%   OUTPUTS: 
%   C           Cost matrix 

freqGrid = freqGrid(:);
pitchGrid = pitchGrid(:).';

maxHarm = 10;
if nargin >= 3 && ~isempty(varargin{1})
    maxHarm = varargin{1};
end

    n = round(freqGrid ./ pitchGrid);
    n = max(n, 1);
    n = min(n, maxHarm);
    
    C = (freqGrid ./ pitchGrid - n).^2;
end