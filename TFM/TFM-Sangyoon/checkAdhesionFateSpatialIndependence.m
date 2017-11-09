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
% nFrames = tracksNA(1).iFrame(end);
matureID = arrayfun(@(x) (x.maturing==1),tracksNA);% | x.maturing==3),tracksNA);
failID = arrayfun(@(x) (x.maturing==0),tracksNA);% | x.maturing==2),tracksNA);
% % collect from only emerging frames
% % equal sampling between from matureID and failID
% nMatureTotal = sum(matureID);
% nFailTotal = sum(failID);
% [equalSampleN,matureOrfail] = min([nMatureTotal,nFailTotal]);
% if matureOrfail==1 %nMatureTotal<nFailTotal
%     sampleMatureID = matureID;
%     sampleFailID = false(size(failID));
%     failIDidx = find(failID);
%     q=0;
%     unequal = true;
%     while unequal
%         idx = round(nFailTotal*rand());
%         if idx>0 && ~sampleFailID(failIDidx(idx))
%             sampleFailID(failIDidx(idx))=true;
%             q=q+1;
%             if q>equalSampleN
%                 unequal = false;
%             end
%         end
%     end
% else %nMatureTotal>nFailTotal
%     sampleMatureID = false(size(matureID));
%     sampleFailID = failID;
%     matureIDidx = find(matureID);
%     q=0;
%     unequal = true;
%     while unequal
%         idx = round(nMatureTotal*rand());
%         if idx>0 && ~sampleMatureID(matureIDidx(idx))
%             sampleMatureID(matureIDidx(idx))=true;
%             q=q+1;
%             if q>equalSampleN
%                 unequal = false;
%             end
%         end
%     end
% end
% 
% nEmergingFrame = arrayfun(@(x) x.emergingFrame,tracksNA(sampleMatureID | sampleFailID));
% nEmergingFrame = unique(nEmergingFrame);

% there was a temporal variation. now we have to collect the data from
% frames overlapping existance of both maturing and failing adhesions
% These will be used for neighbor sampling too.
% get the frames associated with matureID
nMatureFrame = [];
matureIDidx = find(matureID);
for ii=matureIDidx'
    NAFrames = find(cellfun(@(x) strcmp(x,'NA'),tracksNA(ii).state));%,'UniformOutput', false);
    nMatureFrame = [nMatureFrame, NAFrames];
end
nMatureFrame = unique(nMatureFrame);
% get the frames associated with failID
nFailFrame = [];
failIDidx = find(failID);
for ii=failIDidx'
    NAFrames = find(cellfun(@(x) strcmp(x,'NA'),tracksNA(ii).state));%,'UniformOutput', false);
    nFailFrame = [nFailFrame, NAFrames];
end
nFailFrame = unique(nFailFrame);
commonFrames = intersect(nMatureFrame,nFailFrame);

p=0;
% for ii = nEmergingFrame'%1:nFrames
for ii = commonFrames%1:nFrames
    p=p+1;
    % maturing adhesion at ii-th point
    nMature = sum(arrayfun(@(x) (x.presence(ii)),tracksNA(matureID)));
    nFail = sum(arrayfun(@(x) (x.presence(ii)),tracksNA(failID)));
    popStat.nMature(p) = nMature;
    popStat.nFail(p) = nFail;
    popStat.mRatio(p) = nMature/(nMature+nFail);
end
% popStat.nMature = sum(matureID);
% popStat.nFail = sum(failID);
% popStat.mRatio = popStat.nMature/(popStat.nMature+popStat.nFail);
%% Neighbor sample statistics - failing adhesions
% for failing adhesions
p = 0;
usedFrames=[];
for ii = find(failID)'
    p = p+1;
    failNeiStat.id(p) = ii;
    % emerging frame
%     emergingFrame = tracksNA(ii).emergingFrame;
    % all frames with NA states
    NAFrames = find(cellfun(@(x) strcmp(x,'NA'),tracksNA(ii).state));%,'UniformOutput', false);
%     if isempty(emergingFrame) % maturing==2 (long NA)
%         emergingFrame = find(tracksNA(ii).presence,1);
%     end
    % find the adhesions that are present
    nNei = 0;
    failNeiStat.neiIDs{p} = [];
    failNeiStat.neiDist{p} = [];
    failNeiStat.nMature(p) = 0;
    failNeiStat.nFail(p) = 0;
    failNeiStat.mRatio(p) = NaN; 
    NAFrames = setdiff(intersect(commonFrames,NAFrames),usedFrames);
    for jj=NAFrames
%         presentID = arrayfun(@(x) (x.presence(NAFrames)),tracksNA);
        presentID = arrayfun(@(x) (logical(x.presence(jj))),tracksNA);
    %     presentID(ii) = false; % exclude the self
        % find the neighbors in r
        [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(jj),tracksNA(presentID)),...
            arrayfun(@(x) x.yCoord(jj),tracksNA(presentID))],...
            [tracksNA(ii).xCoord(jj),tracksNA(ii).yCoord(jj)],r);
%         [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(NAFrames),tracksNA(presentID)),...
%             arrayfun(@(x) x.yCoord(NAFrames),tracksNA(presentID))],...
%             [tracksNA(ii).xCoord(NAFrames),tracksNA(ii).yCoord(NAFrames)],r);
        if ~isempty(idx{1})
            % inspect their fates
    %         nNei = length(idx{1});
            presentIDind = find(presentID);
            neiMatureID = presentIDind(idx{1});
            nMature = sum(arrayfun(@(x) x.maturing==1,tracksNA(neiMatureID)));
            nFail = sum(arrayfun(@(x) (x.maturing==0),tracksNA(neiMatureID)));% | x.maturing==2),tracksNA(neiMatureID)));
            nNei = nNei + nMature + nFail;
            failNeiStat.neiIDs{p} = [failNeiStat.neiIDs{p}; neiMatureID];
            failNeiStat.neiDist{p} = [failNeiStat.neiDist{p}; dist{1}];
            failNeiStat.nMature(p) = failNeiStat.nMature(p) + nMature;
            failNeiStat.nFail(p) = failNeiStat.nFail(p) + nFail;
            failNeiStat.mRatio(p) = failNeiStat.nMature(p)/nNei;
        else
            display(['No neighbors in r= ' num2str(round(r)) ' pixel for an adhesion at x= ' num2str(round(tracksNA(ii).xCoord(jj))) ', y= ' num2str(round(tracksNA(ii).yCoord(jj))) ' of id (' num2str(ii) ')'])
        end
    end
    usedFrames = [usedFrames, NAFrames];
end
%% for matruing adhesions
p = 0;
usedFrames=[];
for ii = find(matureID)'
    p = p+1;
    matureNeiStat.id(p) = ii;
    % emerging frame
%     NAFrames = tracksNA(ii).emergingFrame;
%     if isempty(NAFrames) % maturing==3 (already FC or FA)
%         NAFrames = find(tracksNA(ii).presence,1);
%     end
    NAFrames = find(cellfun(@(x) strcmp(x,'NA'),tracksNA(ii).state));%,'UniformOutput', false);
%     if isempty(emergingFrame) % maturing==2 (long NA)
%         emergingFrame = find(tracksNA(ii).presence,1);
%     end
    % find the adhesions that are present
    nNei = 0;
    matureNeiStat.neiIDs{p} = [];
    matureNeiStat.neiDist{p} = [];
    matureNeiStat.nMature(p) = 0;
    matureNeiStat.nFail(p) = 0;
    matureNeiStat.mRatio(p) = NaN; 
    NAFrames = setdiff(intersect(commonFrames,NAFrames),usedFrames);
    for jj=NAFrames
%         presentID = arrayfun(@(x) (x.presence(NAFrames)),tracksNA);
        presentID = arrayfun(@(x) logical(x.presence(jj)),tracksNA);
    %     presentID(ii) = false; % exclude the self
        % find the neighbors in r
        [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(jj),tracksNA(presentID)),...
            arrayfun(@(x) x.yCoord(jj),tracksNA(presentID))],...
            [tracksNA(ii).xCoord(jj),tracksNA(ii).yCoord(jj)],r);
%         [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(NAFrames),tracksNA(presentID)),...
%             arrayfun(@(x) x.yCoord(NAFrames),tracksNA(presentID))],...
%             [tracksNA(ii).xCoord(NAFrames),tracksNA(ii).yCoord(NAFrames)],r);
        if ~isempty(idx{1})
            % inspect their fates
    %         nNei = length(idx{1});
            presentIDind = find(presentID);
            neiMatureID = presentIDind(idx{1});
            nMature = sum(arrayfun(@(x) x.maturing==1,tracksNA(neiMatureID)));
            nFail = sum(arrayfun(@(x) (x.maturing==0),tracksNA(neiMatureID)));% | x.maturing==2),tracksNA(neiMatureID)));
            nNei = nNei + nMature + nFail;
            matureNeiStat.neiIDs{p} = [matureNeiStat.neiIDs{p}; neiMatureID];
            matureNeiStat.neiDist{p} = [matureNeiStat.neiDist{p}; dist{1}];
            matureNeiStat.nMature(p) = matureNeiStat.nMature(p) + nMature;
            matureNeiStat.nFail(p) = matureNeiStat.nFail(p) + nFail;
            matureNeiStat.mRatio(p) = matureNeiStat.nMature(p)/nNei;
        else
            display(['No neighbors in r= ' num2str(round(r)) ' pixel for an adhesion at x= ' num2str(round(tracksNA(ii).xCoord(jj))) ', y= ' num2str(round(tracksNA(ii).yCoord(jj))) ' of id (' num2str(ii) ')'])
        end
    end
%     % find the adhesions that are present
%     presentID = arrayfun(@(x) (x.presence(NAFrames)),tracksNA);% & x.maturing<2),tracksNA);
% %     presentID(ii) = false; % exclude the self
%     % find the neighbors in r
%     [idx, dist] = KDTreeBallQuery([arrayfun(@(x) x.xCoord(NAFrames),tracksNA(presentID)),...
%         arrayfun(@(x) x.yCoord(NAFrames),tracksNA(presentID))],...
%         [tracksNA(ii).xCoord(NAFrames),tracksNA(ii).yCoord(NAFrames)],r);
%     if ~isempty(idx{1})
%         % inspect their fates
% %         nNei = length(idx{1});
%         presentIDind = find(presentID);
%         neiMatureID = presentIDind(idx{1});
%         nMature = sum(arrayfun(@(x) x.maturing==1,tracksNA(neiMatureID)));
%         nFail = sum(arrayfun(@(x) (x.maturing==0),tracksNA(neiMatureID)));% | x.maturing==2),tracksNA(neiMatureID)));
%         nNei = nMature + nFail;
% %         nFail = nNei - nMature;
%         matureNeiStat.neiIDs{p} = neiMatureID;
%         matureNeiStat.neiDist{p} = dist{1};
%         matureNeiStat.nMature(p) = nMature;
%         matureNeiStat.nFail(p) = nFail;
%         matureNeiStat.mRatio(p) = nMature/nNei;
%     else
%         display(['No neighbors in r= ' num2str(round(r)) ' pixel for an adhesion at x= ' num2str(round(tracksNA(ii).xCoord(NAFrames))) ', y= ' num2str(round(tracksNA(ii).yCoord(NAFrames))) ' of id (' num2str(ii) ')'])
%         matureNeiStat.neiIDs{p} = [];
%         matureNeiStat.neiDist{p} = [];
%         matureNeiStat.nMature(p) = 0;
%         matureNeiStat.nFail(p) = 0;
%         matureNeiStat.mRatio(p) = NaN;
%     end
    usedFrames = [usedFrames, NAFrames];
end

save([dataPath filesep 'stats.mat'], 'tracksNA', 'matureNeiStat','failNeiStat','popStat')

end