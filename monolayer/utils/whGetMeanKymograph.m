% Input: cell array of kymographs
% Output: average kymograph
function [meanKymograph] = whGetMeanKymograph(kymographs)
N = length(kymographs);
assert(N > 0);

sumKymograph = kymographs{1};

for i = 2 : N
    sumKymograph = sumKymograph + kymographs{i};
end

meanKymograph = sumKymograph ./ N;
end