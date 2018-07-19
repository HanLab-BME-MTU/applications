function weights = lftHist_correctionVector(N)
% Weighting vector for lifetime histogram. Longer lifetimes are less likely to be 
% observed in their entirety over the movie duration. This introduces a bias in the
% lifetime histogram, which can be corrected by the weighting function N/(N-t-1),
% where N is the longest observable lifetime, i.e., the lifetime histogram extends
% from 1 to MovieLength-2 frames.
% Single frame events are thus observable in all N frames, two-frame events in N-1/N frames etc.
% 
% Last modified by Francois Aguet, 09/30/2010
weights = N./(N:-1:1);