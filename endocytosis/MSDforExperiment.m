function [experiment] = MSDforExperiment(experiment)

% MSDforExperiment adds the mean squared dsiplacement of all tracks in 
% movie to experiment structure
%
% INPUT:   experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%    
% OUTPUT
%           clusterResults.msd =  msd from first frame to each consecutive
%           frame for each track in movie
%
% Uses:
%       MSDfromMPM_general
%
% Daniel Nunez, updated March 27, 2009

od = pwd;

for iexp = 1:length(experiment)

    waitHandle = waitbar(iexp/length(experiment),['running movie ' num2str(iexp) ' out of ' num2str(length(experiment))]);

    %LOAD TRACKINFO
    cd([experiment(iexp).source filesep 'TrackInfoMatrices'])
    load('trackInfo.mat')

    %CREATE MPM FROM TRACKINFO
    indX = 1:8:size(trackInfo,2);
    indY = 2:8:size(trackInfo,2);
    ind = sort([indX indY]);
    mpm = trackInfo(:,ind);

    %CALCULATE MSD FROM MPM
    [msd] = MSDfromMPM_general(mpm);
    experiment(iexp).msd = msd;

    close(waitHandle)

end %for each experiment

cd(od)

end %of function