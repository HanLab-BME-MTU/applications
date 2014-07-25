function sigLevel = corrSigLevel(c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
sigLevel=1.96./sqrt(sum(~isnan(c)));

end

