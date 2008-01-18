function [M2]=swapMaskValues(M1,oldValues,newValues)
% SWAPMASKVALUES swaps or replaces values of a mask/image
%
% DESCRIPTION: In the case of a matrix with only two values, the function
%              simply swaps the values if oldValues and newValues are not
%              given.  If they are given, the function replaces
%              oldValues(i) with newValues(i).
%
%              e.g.) say you have a binary mask and want to replace the 0's
%              with 2's and the 1's by NaN's.  Then you can run
%              swapMaskValues(M1,[0 1],[2 NaN]).  This function
%              conveniently deals with multiple replacements.
%
% SYNOPSIS: [M2]=swapMaskValues(M1,oldValues,newValues)
%
% INPUT: M1             : starting mask (or any array)
%        oldValues (opt): n-vector containing values to change
%        newValues (opt): n-vector containing corresponding new values
%
% OUTPUT: M2: matrix updated with new values
%
% MATLAB VERSION (originally written on): 7.2.0.232 (R2006a) Windows_NT
%
% USERNAME: kathomps
% DATE: 13-Jun-2007
%
%

if nargin==2 || nargin>3 % need 1 or 3 arguments
    error('Wrong number of input arguments');
end
if nargin==1 % assume user wants a simple swap of 2 values
    U=unique(M1); % need to ensure mask only has 2 values to swap
    UNaNs=sum(isnan(U)); % NaN's get counted separately with unique.m
    if length(U)==2 || length(U)-UNaNs==1 % either 2 numbers, or 1 number and a nan
        oldValues=[U(1) U(2)];
        newValues=[U(2) U(1)];
    else
        error('Not simple swap: M1 contains one or more than two values')
    end
end
if nargin==3
    if ~isvector(oldValues) || ~isvector(newValues)
        error('oldValues and newValues must both be in vector format and be nonempty');
    end
    if length(oldValues)~=length(newValues)
        error('oldValues and newValues must have same length')
    end
end

M2=double(M1); % initialize output matrix

nSwaps=length(oldValues);

for i=1:nSwaps
    if ~isnan(oldValues(i))
        indexList=find(M1==oldValues(i));
    else
        indexList=find(isnan(M1));
    end
    M2(indexList)=newValues(i);
end