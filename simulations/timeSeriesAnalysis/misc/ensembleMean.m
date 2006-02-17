function [eMean,errFlag] = ensembleMean(traj)
%ENSEMBLEMEAN calculates the ensemble average of a series at each observation point
%
%SYNOPSIS [eMean,errFlag] = ensembleMean(traj)
%
%INPUT  traj   : Observations of time series to be averaged. An array of 
%                structures traj(1:nTraj).observations, where "observations"
%                has 2 colmuns: value + std. Missing points should be
%                indicated with NaN.
%
%OUTPUT eMean  : Ensemble average of trajectories.
%       errFlag: 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errFlag = 0;
eMean = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments were used when function was called
if nargin ~= 1
    disp('--ensembleMean: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
numTraj = length(traj);
if numTraj <= 1
    disp('--ensembleMean: traj should have more than one entry!');
    errFlag = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mean calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get maximum trajectory length
maxLength = 0;
for i=1:numTraj
    maxLength = max([maxLength size(traj(i).observations,1)]);
end

%construct matrix of trajectories
observationsMat = NaN*ones(maxLength,numTraj);
for i=1:numTraj
    observationsMat(1:length(traj(i).observations(:,1)),i) = traj(i).observations(:,1);
end

%calculate ensemble average at each observation point
eMean = nanmean(observationsMat')';


%%%%% ~~ the end ~~ %%%%%
