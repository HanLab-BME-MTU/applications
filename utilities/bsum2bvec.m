function binvector=bsum2bvec(binsum)
%BSUM2BVEC returns the vector containing the addends of a sum of powers of two, i.e. 4->[4], 6->[2 4]', 7->[1 2 4]'
% 
%
% SYNOPSIS  binvector=bsum2bvec(binsum)
%
% INPUT binsum : scalar, sum of powers of 2
%
% OUTPUT binvector : vector containing the addents of the input, e.g. [1 2 4]' from an input of 7
% 
% DEPENDS ON --
%
% CALLED BY spotID
%   
% c: 26/09/02	Jonas

%test input
if any(size(binsum)>1)
    error('input for bsum2bvec needs to be scalar')
end

%convert to binary string
b=dec2bin(binsum);

%write binvector in an elegant, yet nontrivial way
binvector=nonzeros((b([end:-1:1])=='1').*2.^[0:length(b)-1]);