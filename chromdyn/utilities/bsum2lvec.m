function logicalVector=bsum2lvec(binsum)
%BSUM2LVEC returns the binary representation of a number in a vector 
%
% SYNOPSIS  binvector=bsum2lvec(binsum)
%
% INPUT binsum : scalar, sum of powers of 2
%
% OUTPUT binvector : binary representation of a number in a vector, e.g. [1 0 0 1]' from an input of 9
% 
% DEPENDS ON --
%   
% c: 02/03	Jonas

%test input
if any(size(binsum)>1)
    error('input for bsum2lvec needs to be scalar')
end

%convert to binary string
b=dec2bin(binsum);

%convert to binary vector
logicalVector=str2num(b');
logicalVector=logicalVector(end:-1:1);