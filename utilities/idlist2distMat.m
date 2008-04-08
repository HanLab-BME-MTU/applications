function [distance, distanceUnitVectors, displacement, displacementUnitVectors, idlist, idxLists, intensity] = idlist2distMat(idlist, dataProperties, correctDisplacement, allowEstimatedTags, tagOrder)
%IDLIST2DISTMAT calculates distance matrices for idlists.
%
% SYNOPSIS: [distance, sigmaDistance, unitVectors,...
%               displacement, displacementSigma, displacementUnitVectors...
%               idlist, idxLists, intensity] ...
%               = idlist2distMat(idlist, dataProperties, ...
%                   correctDisplacement, allowEstimatedTags)
%
% INPUT idlist: idlist2
%       dataProperties: dataProperties-structure (see defaultDataProperties
%                       for information)
%       correctDisplacement (opt): correct the displacement
%                       {0} use the setting in dataProperties.linker_useCOM
%                           If there is no setting, use COM
%                        1  use COM
%                       -1  don't correct
%       allowEstimatedTags (opt): whether or not to include estimated tags
%                       in the calculations. [{0}/1]
%       tagOrder (opt): Because the output matrices contain both values
%                       and variances, reordering will be difficult. Use
%                       tagOrder to indicate how you want the tags to be
%                       ordered in the output. Default: (1:nTags)
%
% OUTPUT distance: nTags x nTags x nTimepoints distance array (estimated
%               spots will not be counted). Below the diagonal is the
%               distance, above the diagonal the uncertainty of the
%               distance
%		 distanceUnitVectors: nTags x nTags x nTimepoints x 3 array of unit
%               vectors for every distance measurement
%        displacement : (nTimepoints-1) x nTags x 2 array of displacements
%               (:,:,1), and uncertainties (:,:,2)
%        displacementUnitVectors : nTimepoints x nTags x 3 unit vectors of
%               displacements
%        idlist: idlist without "bad" tags (= no single occurences)
%        idxLists: structure of lists with fields:
%                   - estimatedTag: nTimepoints x nTags array with 1
%                      wherever a tag has been estimated 
%                   - goodTimes: list of timePoints of non-deleted frames
%                   - isGoodTime: nTimepoints x 1 array with ones wherever
%                      there are goodTimes (goodTimes = find(isGoodTime))
%                   - isTracked: nTimepoints x nTags array with 1 for
%                       tracked tags
%                   - isFusion : nTimepoints x nTags array with tag numbers
%                      indicating to which other tag the current tag is
%                      fused. !!!! Currently, this doesn't support multiple
%                      fusions, and it doesn't distinguish between primary
%                      and secondary fusions.
%                   - goodTagIdx : index into the tags that have been used
%                      to calculate distances
%        intensity: list of tag intensities
%
% REMARKS Estimated spots, deleted frames are all converted into NaNs!
%           (I know that displacements could be packed in the diagonal of
%           the distance matrix. As of now, it's too much of a pain to read
%           the diagonal of a n-d array.
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 27-Apr-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==================================
%% TEST INPUT AND READ IDLIST DATA
%==================================

% check input
if nargin < 1 || isempty(idlist) || ~checkIdlist(idlist,1)
    error('missing or old idlist')
end
if nargin < 2 || isempty(dataProperties)
    error('missing dataProperties')
end
if nargin < 3 || isempty(correctDisplacement)
    if isfield(dataProperties,'linker_useCOM')
        correctDisplacement = dataProperties.linker_useCOM;
    else
        correctDisplacement = 0;
    end
elseif correctDisplacement == -1
    correctDisplacement = 0;
end
if nargin < 4 || isempty(allowEstimatedTags)
    allowEstimatedTags = 0;
end

% read some idlist data
nTimepoints = length(idlist);
linklists = cat(3,idlist.linklist);
goodTimes = squeeze(linklists(1,1,:));
centroids = cat(1,idlist(goodTimes).centroid);

% remove bad tags - also from the Q-matrices. Remember goodTagIdx
badTagIdx = find(ismember(linklists(:,5,1),[2,3]));
goodTagIdx = missingIndices(badTagIdx,size(linklists,1));
linklists(badTagIdx,:,:) = [];

idlist = LG_deleteTag(idlist, badTagIdx, goodTimes);
% count tags only after deleting tags!!
nTags = size(idlist(goodTimes(1)).linklist,1);

% check tagOrder
if nargin < 5 || isempty(tagOrder)
    tagOrder = (1:nTags)';
end

%% read Q-matrix diagonals

% find which Q to use. Have primary and secondary qName, in case there is a
% frame that has not been successfully tracked from any source
if isfield(idlist(goodTimes(1)).info,'totalQ_Pix') &&...
        ~isempty(idlist(goodTimes(1)).info.totalQ_Pix)
    qName{1} = 'totalQ_Pix';
else
    qName{1} = 'detectQ_Pix';
end
qName{2} = 'detectQ_Pix'; % there's no alternative for detectQ

% get pix2mu2 - squared conversion factor from pixels to microns. Repmat so
% that we can directly multiply all Q-diag-entries at once
pix2mu2 = repmat([dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z].^2,nTags,1);

% preassing qMatrixDiags so that we only have to loop through goodTimes
qMatrixDiags = repmat(NaN,[nTimepoints,nTags,3]);
for t = goodTimes'
    if isfield(idlist(t).info,qName{1}) && ~isempty(idlist(t).info.(qName{1}))
    qMatrixDiags(t,:,:) = reshape(reshape(...
        full(diag(idlist(t).info.(qName{1}))) .*...
        repeatEntries(idlist(t).linklist(:,12),3), [3,nTags])'.*...
        pix2mu2,[1,nTags,3]);
    else
        % if no qName1, use secondary Q
        qMatrixDiags(t,:,:) = reshape(reshape(...
        full(diag(idlist(t).info.(qName{2}))) .*...
        repeatEntries(idlist(t).linklist(:,12),3), [3,nTags])'.*...
        pix2mu2,[1,nTags,3]);
    end
end



% order tags
linklists = linklists(tagOrder,:,:);
qMatrixDiags = qMatrixDiags(:,tagOrder,:);



% read intensity
intensity = squeeze(linklists(:,8,:))';



%===============================

%===============================
%% CALCULATE DISTANCE MATRIX
%===============================

% loop through all the pairs of tags. For each pair, read coordinates and
% covariances, and get the data from deltaCoordinates

% preassign stuff. Displacement has one fewer entry!
distance = repmat(NaN,[nTags, nTags, nTimepoints]);
distanceUnitVectors = repmat(NaN,[nTags, nTags, nTimepoints, 3]);
displacement = deal(repmat(NaN,[nTimepoints-1, nTags, 2]));
displacementUnitVectors = repmat(NaN,[nTimepoints-1, nTags, 3]);

coco = repmat(NaN, [nTimepoints, 3]);
points([1,2]) = struct('coordinates',coco,'covariances',coco);

% preassign idxLists
idxLists = struct('estimatedTag',false(nTimepoints, nTags),...
    'isTracked',false(nTimepoints,nTags), 'goodTimes',goodTimes,...
    'fusedTag',zeros(nTimepoints,nTags), ...
    'isGoodTime',false(nTimepoints,1),...
    'goodTagIdx',goodTagIdx);


% loop. First: calculate displacement. Then assign other tag and calculate
% distance
for tag1 = 1:nTags

    % clear point 1
    [points(1).coordinates, points(1).covariances] = deal(coco);

    % read coordinates and covariances of first tag
    points(1).coordinates(goodTimes,:) = ...
        permute(linklists(tag1,9:11,:),[3,2,1]);
    points(1).covariances = ...
        squeeze(qMatrixDiags(:,tag1,:));

    % remove estimated entries if requested, and store data in idxLists
    estimatedIdx = goodTimes(linklists(tag1,3,:) == 1);
    idxLists.estimatedTag(estimatedIdx,tag1) = true;
    if ~allowEstimatedTags
        points(1).coordinates(estimatedIdx,:) = NaN;
        points(1).covariances(estimatedIdx,:) = NaN;
    end
    
    % check for tracked tags
    isTrackedIdx = goodTimes(ismember(linklists(tag1,3,:),[2,5]));
    idxLists.isTracked(isTrackedIdx,tag1) = true;

    % don't calculate displacement quite yet - we might want to correct
    % the positions!

    % calculate distance
    for tag2 = 1:tag1-1

        % clear point 2
        [points(2).coordinates, points(2).covariances] = deal(coco);

        % read coordinates and covariances of second tag
        points(2).coordinates(goodTimes,:) = ...
            permute(linklists(tag2,9:11,:),[3,2,1]);
        points(2).covariances = ...
            squeeze(qMatrixDiags(:,tag2,:));

        % remove estimated entries if necessary. We don't update idxList
        % here, because in the first loop, we go through all tags already
        if ~allowEstimatedTags
            estimatedIdx = goodTimes(linklists(tag2,3,:) == 1);
            points(2).coordinates(estimatedIdx,:) = NaN;
            points(2).covariances(estimatedIdx,:) = NaN;
        end

        % calculate distance
        [d, sigmaD, unitVector] =  deltaCoordinates(points);

        % fill in variables
        distance(tag1, tag2, :) = permute(d,[2,3,1]);
        distance(tag2, tag1, :) = permute(sigmaD,[2,3,1]);
        distanceUnitVectors(tag1, tag2,:,:) = permute(unitVector,[3,4,1,2]);
        distanceUnitVectors(tag2, tag1,:,:) = - distanceUnitVectors(tag1,tag2,:,:);
    end

    % correct Pos1
    switch correctDisplacement
        case 1 % COM
            points(1).coordinates(goodTimes,:) = points(1).coordinates(goodTimes,:) - centroids;
        case 0
            % don't correct
        otherwise
            error('sorry, correction %i has not been implemented yet',...
                correctDisplacement);
    end


    % calculate displacement
    [d, sigmaD, unitVector] =  deltaCoordinates(points(1));

    % fill in variables
    displacement(:,tag1,1) = d;
    displacement(:,tag1,2) = sigmaD;
    displacementUnitVectors(:,tag1,:) = permute(unitVector,[1,3,2]);

end

%==============================
%% COLLECT INDICES
%==============================

% estimatedTag, goodTimes has already been collected. Also isTracked

% for fusedTag: Loop through distanceMatrices, and find zero distance. I
% don't quite know what the uncertainty will be. Therefore, I check for
% every whether row or col has a zero.
for i=1:nTags
    [tagNum,t] = find(squeeze(distance(i,:,:) == 0) | squeeze(distance(:,i,:) == 0));
    idxLists.fusedTag(t,i) = tagNum;
end

% isGoodTime is a logical representation of goodTimes
idxLists.isGoodTime(goodTimes) = true;