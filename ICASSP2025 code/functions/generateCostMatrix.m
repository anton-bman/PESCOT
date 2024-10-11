function C = generateCostMatrix(freqGrid, pitchGrid, varargin)


nFreqsGrid = length(freqGrid);
nPitchGrid = length(pitchGrid);

doInf = 0;
if nargin >= 3 && ~isempty(varargin{1})
    maxHarmonics = varargin{1};
    if maxHarmonics == Inf
        doInf = 1;
    end
end


if doInf == 0
    C = zeros(nFreqsGrid, nPitchGrid);
    for l = 1:nPitchGrid
        C(:,l) = min((freqGrid(:)./pitchGrid(l) - (1:maxHarmonics)).^2,[],2);
    end
else
    C = zeros(nFreqsGrid, nPitchGrid);
    for l = 1:nPitchGrid
        temp1 = zeros(nFreqsGrid,1);
        temp2 = zeros(nFreqsGrid,1);
        normFreqs = freqGrid(:) ./ pitchGrid(l);
        octaveNbr = round(normFreqs);
        octaveNbrGz = (octaveNbr == 0);

        if sum(octaveNbrGz,'all') > 0
            temp1(octaveNbrGz) = abs(normFreqs(octaveNbrGz) - 1).^2;
            temp1(~octaveNbrGz) = (mod(normFreqs(~octaveNbrGz), round(normFreqs(~octaveNbrGz)))).^2;

            temp2(octaveNbrGz) = ((octaveNbr(octaveNbrGz) - mod(normFreqs(octaveNbrGz), octaveNbr(octaveNbrGz)))).^2*inf;
            temp2(~octaveNbrGz) = ((octaveNbr(~octaveNbrGz) - mod(normFreqs(~octaveNbrGz), octaveNbr(~octaveNbrGz)))).^2;
        else
            temp1 = (mod(normFreqs, octaveNbr)).^2;
            temp2 = ((octaveNbr - mod(normFreqs, octaveNbr))).^2;
        end
        %find min distance to nearest octave nbr
        minDist = min(temp1, temp2);
        C(:, l) = minDist;
    end

end
end