function weights = lftHist_correctionVector(N)
% Weighting vector for lifetime histogram. Longer lifetimes are less likely to be 
% observed in their entirety over the movie duration. This introduces a bias in the
% lifetime histogram, which can be corrected by a weighting function, given by
% N/(N-t-1), where N is the movie length.
% Single frame events are observable in N-2/N frames, two-frame events in N-3/N frames etc.
% Events of duration N-1 and N cannot be observed, and are set to zero.
% 
% Last modified by Francois Aguet, 09/30/2010
weights = [N./(N-2:-1:1) 0 0];