function out = returnRightVector(vector, nFix, rowOrCol)
%RETURNRIGHTVECTOR checks if n-by-x has the correct size and transposes it if necessary
% 
% For many functions, it is necessary that inputs are of the form n-by-x or
% x-by-n, where x is a fixed number of columns/rows. RETURNRIGHTVECTOR is a
% utility function for checking such input and transposing it if necessary.
% It will return an error:bad input size if the input can not be corrected
% and warning:square input matrix, if n=x.
%
% SYNOPSIS out = returnRightVector(vector, nFix, rowOrCol)
%
% INPUT     vector: input to check
%           nFix  : (opt, {1}) number of fixed rows/cols in the output (x)
%           rowOrCol : (opt, ['r'/{'c'}]) whether output should be in rows
%                      (x-by-n) or columns (n-by-x).
%
% OUTPUT    out   : output vector
%
%c: 9/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=====================
% CHECK INPUT
%=====================

% defaults
defNFix = 1;
defRowOrCol = 'c';

% check input
if nargin == 0 || isempty(vector)
    error('specify at least input argument for RETURNRIGHTVECTOR');
end

if nargin > 1 && ~isempty(nFix)
    % take nFix
else
    % use default
    nFix = defNFix;
end

if nargin > 2 && ~isempty(rowOrCol) && ischar(rowOrCol)
    % take input
else
    % use default
    rowOrCol = defRowOrCol;
end
% turn rowOrCol into index
switch rowOrCol
    case 'r'
        rocIndex = 1;
    case 'c'
        rocIndex = 2;
    otherwise
        error('rowOrCol has to be either ''r'' or ''c''!')
end

%========================



%========================
% FORMAT VECTOR
%========================

% get vector size
vSize = size(vector);
% test vector size
if length(vSize) > 2
    error('RETURNRIGHTVECTOR can only handle 2-D arrays')
end
if ~any(vSize == nFix)
    error('at least one dimension has to have the fixed lenght')
end
if all(vSize == nFix)
    warning('RETURNRIGHTVECTOR:SQUAREINPUT',...
        'Ambiguous input: The number of rows matches the number of columns. \n The input might have to be transposed')
end

% transpose if necessary
if vSize(rocIndex) ~= nFix
    out = vector';
else
    out = vector;
end

