function [idlist] = linkerFuseTags(idlist,dataProperties, maxTagIdx, constants, verbose, intAx)
%LINKERFUSETAGS fuses tags that are closer than a given distance
%
% SYNOPSIS: [idlist] = linkerFuseTags(idlist,dataProperties, constants, verbose, intAx)
%
% INPUT idlist: idlist
%		dataProperties: dataProperties
%       maxTagIdx: number of tags
%		constants: constants-structure
%		verbose: 1 if plotting enabled
%		intAx axes to the intensity figure
%
% OUTPUT idlist: idlist with fused tags
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 09-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find frames with estimated tags. For every estimated tag: find vector to
% closest other position, store it along with tag to fuse to, and actual
% distance.
% For every timepoint: check that the fused-to tag is unique. If not, only
% take closest estimate for fusion (change in the future??)
% For all tags: calculate Rayleigh-limit. If a tag is closer than 2xRL,
% fuse it if also the tag intensity is higher than the estimated intensity.
% Mark primary fusions as "3", secondary fusions as "4"
% Fusions are considered real positions. Therefore, they get the "correct"
% Q, a spotNumber etc.

% defaults
fuseRatio = constants.fuseRatio;



%% FIND FRAMES WITH ESTIMATED TAGS

% linklists = cat(3,idlist.linklist);
% goodTimes = squeeze(linklists(1,1,:));
% nGoodTimes = length(goodTimes);
%
% % estimatedTags is a length(goodTimes)-by-nGoodTags array
% goodTags = linklists(:,5,1) < 2;
% estimatedTags = squeeze(linklists(goodTags,3,:) == 1)';
%
% % loop through every timepoint

% collect all distances
[distance, distanceUnitVectors, dummy, dummy, dummy, idxLists, intensity] = ...
    idlist2distMat(idlist, dataProperties, 0, 1);
nGoodTags = length(idxLists.goodTagIdx);

% estimatedTagList: 1-t, 2-iTag (among good tags), 3-jTag (among good tags),
% 4-distance, 5-7-distanceVector
estimatedTagList = zeros(sum(idxLists.estimatedTag(:)),7);
eIdx = 1; % index into estimatedTagList

% get distances to estimated tags
for t=idxLists.goodTimes'

    % find all estimated tags for this frame
    estimatedIdx = find(idxLists.estimatedTag(t,:));
    nEstimates = length(estimatedIdx);

    % read list of possible targets
    targetIdx = missingIndices(estimatedIdx,nGoodTags);

    % are there targets?
    if ~isempty(targetIdx)

        % for all targets: calculate actual and projected intensity, find good
        % targets
        expValue = exp(idlist(1).stats.intFit.xFit(end) * t);
        targetIntensities = expValue * ...
            idlist(1).stats.intFit.tagFactor(idxLists.goodTagIdx(targetIdx));
        % good target: intensity larger than expected intensity
        goodTargets = intensity(t==idxLists.goodTimes,targetIdx) > targetIntensities;

        % are there any good targets?
        if any(goodTargets)

            % loop through estimated tags for distances
            for i=1:nEstimates
                % distance is a nGoodTags-by-nGoodTags array with distances below
                % the diagonal and the corresponding variances above the diagonal
                allDist1 = distance(estimatedIdx(i),:,t);
                allDist1 = allDist1(:);
                allDist2 = distance(:,estimatedIdx(i),t);
                allDist = [allDist1(1:estimatedIdx(i)-1);NaN;...
                    allDist2(estimatedIdx(i)+1:end)];

                % find the minimum distance. Mask estimated poitions, first
                allDist(estimatedIdx) = NaN;
                [minDistance, minTagIdx] = nanmin(allDist);

                % check whether minTagIdx points to a good target (intensity-wise)
                % careful: goodTargets is relative to targetIdx
                if goodTargets(minTagIdx==targetIdx)

                    % store stuff
                    rowIdx = max(estimatedIdx(i),minTagIdx);
                    colIdx = min(estimatedIdx(i),minTagIdx);
                    estimatedTagList(eIdx,:) = ...
                        [t, estimatedIdx(i), minTagIdx, minDistance,...
                        squeeze(distanceUnitVectors(rowIdx, colIdx, t, :))'];
                    eIdx = eIdx + 1;


                    % check "recursively" for same fusion target, and remove the one
                    % with larger distance
                    sameTarget = find(estimatedTagList(:,3) == ...
                        estimatedIdx(i) & estimatedTagList(:,1) == t);
                    if length(sameTarget) > 1
                        % remove target with larger distance
                        [dummy, maxIdx] = max(estimatedTagList(sameTarget,4));
                        estimatedTagList(sameTarget(maxIdx),:) = [];
                        eIdx = eIdx - 1;
                    end
                end % check whether we would fuse to a good target
            end % loop through estimated tags
        end % any good targets?
    end % any targets?
end % loop time

% remove empty entries in estimatedTagList
estimatedTagList(eIdx:end,:) = [];


%=================================
%% FUSE PROXIMAL TAGS
%=================================

% calculate Rayleigh-ratio
rayleighRatio = estimatedTagList(:,4)./...
    rayleighFromOri(estimatedTagList(:,5:7),dataProperties.WVL,...
    dataProperties.NA);

fuseTags = find(rayleighRatio < fuseRatio);

% loop through fuseTags and fuse. Fusions are considered real positions.
% Since we need to update Q from the new spotNumbers, revert Q here
idlist = linkerRevertQmatrices(idlist, idxLists.goodTimes);

for i = fuseTags'
    t = estimatedTagList(i,1);
    % primary tag: Tag to which another is being fused
    primaryTag = estimatedTagList(i,3);
    % secondary tag: Tag which is fused to another
    secondaryTag = estimatedTagList(i,2);
    
    % set fusion flag
    if idlist(t).linklist(primaryTag,3) == 2
        idlist(t).linklist(primaryTag,3) = 5;
    else
        idlist(t).linklist(primaryTag,3) = 3;
    end
    idlist(t).linklist(secondaryTag,3) = 4;
    
    % adjust position, chi2, of secondary tag
    idlist(t).linklist(secondaryTag,9:12) = ...
        idlist(t).linklist(primaryTag,9:12);
    
    % set spotNumber
    idlist(t).linklist(secondaryTag,2) = ...
        idlist(t).linklist(primaryTag,2);
    
    % remove trackInit
    trackInitIdx = idlist(t).trackInit(:,1) == secondaryTag;
    if any(trackInitIdx)
        idlist(t).trackInit(trackInitIdx,:) = [];
    end
        
    % if verbose, plot 'F' over the fused intensity
    if verbose
        axes(intAx)
        text(t,intensity(t,secondaryTag),'F');
    end
end

% rebuild Q-matrices
idlist = linkerWriteQmatrices(idlist,idxLists.goodTimes);

%=====================


%========================
%% DON'T RE-ESTIMATE POSITIONS
%========================

% the user should be able to see the difference between "too far away"-tags
% and close tags. 
% If we re-estimate now, we could fuse some more, etc.

% % prepare idlist for recalc first
% idlist(1).stats.recalc{1} = 2;
% [idlist] = linkerReadIdlist(idlist, constants);
% [idlist] = linkerEstimatePositions(idlist, maxTagIdx, ...
%     idxLists.goodTimes, constants, verbose, intAx);
% % write Q again, because readIdlist reverts, and finishIdlist only cares
% % about frames with fewer tags!!
% idlist = linkerWriteQmatrices(idlist,idxLists.goodTimes);








 

