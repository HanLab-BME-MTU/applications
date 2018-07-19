function binsum=lvec2bsum(logicalVector)
%BSUM2LVEC returns the binary representation of a number in a vector 
%
% SYNOPSIS  binsum=lvec2bsum(logicalVector)
%
% INPUT  vector of 1 and 0
%
% OUTPUT binsum: number represented by the binary vector
% 
% DEPENDS ON --
%   
% c: 02/03	Jonas

%test input
if any(logicalVector)>1|any(logicalVector<0)
    error('input has to be 0 or 1')
end

%convert to binary string
b=num2str(logicalVector)';
b=b(end:-1:1);

binsum=bin2dec(b);
