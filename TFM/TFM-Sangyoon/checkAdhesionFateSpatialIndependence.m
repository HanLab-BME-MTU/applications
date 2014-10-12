function [popStat,failNeiStat,matureNeiStat] = checkAdhesionFateSpatialIndependence(tracksNA, r, outputPath)
% [popStat,failNeiStat,matureNeiStat] = checkAdhesionFateSpatialIndependence(tracksNA, r, outputPath)
% checks how indenpendent the neighboring adhesions' fates are against
% overall ratio of all adhesion tracks between maturing vs. failing
% adhesions. separateMatureAdhesionTracks should be run beforehand.
% input
%           tracksNA:           tracks that have member of .maturing
%           (created by separateMatureAdhesionTracks after colocalizationAdhesionsWithTFM)
%           r:                       neighbor radius
%           
% Basically for neighbors of every adhesion with a certain radius, we store
% the number of maturing adhesions, the number of failing adhesions, and
% the ratio between them in a structure.
% output
%           popStat:            statistics about whole popoulation containing,
%                                   .nMature: the number of maturing adhesions
%                                   .nFail: the number of failing adhesions
%                                   .mRatio: nMature/(nMature+nFail)
%           failNeiStat:        statistics about neighbors of failing tracks containing,
%                                   .id: corresponding adhesion id in tracksNA
%                                   .neiIDs: neighboring adhesions in tracksNA (cell array).
%                                   .neiDist: neighbor distance
%                                   .nMature: the number of maturing adhesions in neighbors.
%                                   .nFail: the number of failing adhesions
%                                   .mRatio: nMature/(nMature+nFail)
%           matureNeiStat:        statistics about neighbors of maturing tracks
%           
% Sangyoon Han October 2014

%% Set up the output file path
dataPath = [outputPath filesep 'adhesionFateAnalysis'];
if ~exist(dataPath,'dir') 
    mkdir(dataPath);
end
%% Statistics about population
% get the statistics from each time point
nFrames = tracksNA(1).iFrame(end);
matureID = arrayfun(@(x) (x.maturing==1),tracksNA);
failID = arrayfun(@(x) (x.maturing==0),tracksNA);
for ii = 1:nFrames
    % maturing adhesion at ii-th point
    nMature = sum(arrayfun(@(x) (x.presence(ii)==true),tracksNA(matureID)));
    nFail = sum(arrayfun(@(x) (x.presence(ii)==true),tracksNA(failID)));
    popStat.nMature(ii) = nMature;
    popStat.nFail(ii) = nFail;
    popStat.mRatio(ii) = nMature/(nMature+nFail);
end
% popStat.nMature = sum(matureID);
% popStat.nFail = sum(failID);
% popStat.mRatio = popStat.nMature/(popStat.nMature+popStat.nFail);
%% Neighbor sample stati
% for failing adhesions
p = 0;
for ii = find(failID)'
    p = p+1;
    failNeiStat.id(p) = ii;
    % emerging frame
    emergingFrame = tracksNA(ii).emergingFrame;
    % find the adhesions that are present
    presentID = arrayfun(@(x) (x.presence(emergingFrame)==true),tracksNA);
    % find the neighbors in r
    [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(emergingFrame),tracksNA(presentID)),...
        arrayfun(@(x) x.yCoord(emergingFrame),tracksNA(presentID))],...
        [tracksNA(ii).xCoord(emergingFrame),tracksNA(ii).yCoord(emergingFrame)],r);
    if ~isempty(idx{1})
        % inspect their fates
        nNei = length(idx{1});
        presentIDind = find(presentID);
        neiMatureID = presentIDind(idx{1});
        nMature = sum(arrayfun(@(x) x.maturing==1,tracksNA(neiMatureID)));
        nFail = nNei - nMature;
        failNeiStat.neiIDs{p} = neiMatureID;
        failNeiStat.neiDist{p} = dist{1};
        failNeiStat.nMature(p) = nMature;
        failNeiStat.nFail(p) = nFail;
        failNeiStat.mRatio(p) = nMature/nNei;
    else
        disp(['No neighbors in r=' num2str(r) ' pixel for an adhesion at x=' num2str(tracksNA(ii).xCoord(emergingFrame)) ', y=' tracksNA(ii).yCoord(emergingFrame) ' of id (' num2str(ii) ')'])
        failNeiStat.neiIDs{p} = [];
        failNeiStat.neiDist{p} = [];
        failNeiStat.nMature(p) = 0;
        failNeiStat.nFail(p) = 0;
        failNeiStat.mRatio(p) = NaN;
    end
end
% for matruing adhesions
p = 0;
for ii = find(matureID)'
    p = p+1;
    matureNeiStat.id(p) = ii;
    % emerging frame
    emergingFrame = tracksNA(ii).emergingFrame;
    % find the adhesions that are present
    presentID = arrayfun(@(x) (x.presence(emergingFrame)==true),tracksNA);
    % find the neighbors in r
    [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(emergingFrame),tracksNA(presentID)),...
        arrayfun(@(x) x.yCoord(emergingFrame),tracksNA(presentID))],...
        [tracksNA(ii).xCoord(emergingFrame),tracksNA(ii).yCoord(emergingFrame)],r);
    if ~isempty(idx{1})
        % inspect their fates
        nNei = length(idx{1});
        presentIDind = find(presentID);
        neiMatureID = presentIDind(idx{1});
        nMature = sum(arrayfun(@(x) x.maturing==1,tracksNA(neiMatureID)));
        nFail = nNei - nMature;
        matureNeiStat.neiIDs{p} = neiMatureID;
        matureNeiStat.neiDist{p} = dist{1};
        matureNeiStat.nMature(p) = nMature;
        matureNeiStat.nFail(p) = nFail;
        matureNeiStat.mRatio(p) = nMature/nNei;
    else
        disp(['No neighbors in r=' num2str(r) ' pixel for an adhesion at x=' num2str(tracksNA(ii).xCoord(emergingFrame)) ', y=' tracksNA(ii).yCoord(emergingFrame) ' of id (' num2str(ii) ')'])
        matureNeiStat.neiIDs{p} = [];
        matureNeiStat.neiDist{p} = [];
        matureNeiStat.nMature(p) = 0;
        matureNeiStat.nFail(p) = 0;
        matureNeiStat.mRatio(p) = NaN;
    end
end

save([dataPath filesep 'stats.mat'], 'tracksNA', 'matureNeiStat','failNeiStat','popStat')

end