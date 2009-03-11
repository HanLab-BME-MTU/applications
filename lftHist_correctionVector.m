function[bincorvec]=lftHist_correctionVector(length)
% calculates the correction vector for the lifetime histogram, based on the
% fact that longer lifetimes are more likely to be cut off at the beginning
% or end which will result in them not being counted

% vector of possible lifetimes
redlen = length-2;
taus = [1:redlen];

% detection probability of a given lifetime
dp_tau = (redlen-taus+1) / length;

% bincorvec
bincorvec = 1./dp_tau;

end % of function