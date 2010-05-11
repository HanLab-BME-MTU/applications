function weights = lftHist_correctionVector(N)
% Weighting vector for lifetime histogram. Longer lifetimes are less likely to be 
% observed in their entirety over the movie duration. This introduces a bias in the
% lifetime histogram, which can be corrected by a weighting function, given by
% N/(N-t-1), where N is the movie length.
% The time vector is defined as t in [-1,0,...,N-2]. The first frame yields no
% useable events, time t=0 corresponds to the second frame.
% 
% Last modified by Francois Aguet, 05/11/2010

weights = N./(N:-1:1);