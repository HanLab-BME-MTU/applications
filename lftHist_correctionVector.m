function [corrVect] = lftHist_correctionVector(N)
% Correction vector for the lifetime histogram. Longer lifetimes are more likely
% to be cut off at the beginning or end which will result in them not being counted
%
% The weighting function is N/(N-x), where N is the movie length.
% To do: current function truncates the sequence, should be fixed.
% 
% Last modified by Francois Aguet, 02/18/2010

corrVect = N./(N-2:-1:1);