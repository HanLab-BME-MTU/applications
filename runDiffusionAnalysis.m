function [] = runDiffusionAnalysis;

% runDiffusion analyses diffusion of all tracks for each chosen movie
%
% INPUT
%                    
%
% OUTPUT
%
% REMARKS   The function does not overwrite previous results; it uses
%           secureSave. It makes a new folder called diffusionAnalysis
%           under movie folder. It uses the function
%           trackDiffusionAnalysis1 to do the diffusion analysis.
%
% Uses: trackDiffusionAnalysis1
%       loadIndividualMovies
%       plotPopulationProbabilityDistribution
%
% Daniel Nunez, July 9, 2008
%   updated September 1, 2008 to no longer restrict tracks analyzed

%%
%LOAD MOVIES
[experiment] = loadIndividualMovies();

%FOR EACH MOVIE
for iexp = 1:length(experiment)

    %print movie number
    display(['running movie number ' num2str(iexp) ' out of ' num2str(length(experiment))]);

    %LOAD TRACKS
    cd(experiment(iexp).source);
    cd('TrackInfoMatrices');
    load trackInfo.mat;
    if exist('trackInfo')
        trackInfoMat = trackInfo;
    end

    %PREPARE TRACKS FOR ANALYSIS
    %zeros for position are interpreted not as missing tracks, but as point
    %on the vbertex; these are therefore made NaNs. The other non-position
    %columns left as zeros. This assumes that zeros for positions are
    %missing tracks.
    tracks = full(trackInfoMat);
    tracks(tracks == 0) = nan;
    tracks(:,3:8:end-5) = 0;
    tracks(:,5:8:end-3) = 0;
    tracks(:,6:8:end-2) = 0;
    tracks(:,7:8:end-1) = 0;
    tracks(:,8:8:end) = 0;

    %RUN ANALYSIS
    alphaValues = [0.05 0.1];
    checkAsym = 1;
    [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,1,2,checkAsym,alphaValues,1);

    %SAVE RESULTS
    mkdir([experiment(iexp).source filesep 'diffusionAnalysis']);
    secureSave([experiment(iexp).source filesep 'diffusionAnalysis' filesep 'diffusionAnalysisResultsForAllTracks'],'diffAnalysisRes','errFlag','alphaValues','checkAsym');
    
    %clear trackInfo so that it doesn't interfere with the trackInfo of the
    %next movie (because there are two names for this variable depending on
    %the time the movie was tracked, it is possiblt for the analysis to
    %continue with the trackInfo from a previous movie..giving repeated
    %results and skiping the movie)
    clear trackInfo;
    clear trackInfoMat;
end %of for each movie

end %of function