function [idlist,positions,sigmaPositions, nsp] = slist2idlist(slist, projectNumber, dataProperties)
%SLIST2IDLIST generates an idlist from a slist with known positions

if isscalar(projectNumber)
truePositions = loadPositions(projectNumber);
else
    truePositions = projectNumber;
end

% read constants etc from input
sp = size(truePositions);
nTimepoints = sp(1);
if size(sp)<3
    nTags = 1;
else
    nTags = sp(3);
end

% augment slist with field COM
meanPos = mean(truePositions(:,1:3,:),3);
goodTimes = zeros(nTimepoints,1);
for t = 1:nTimepoints
    if ~isempty(slist(t).sp)
        goodTimes(t) = 1;
        slist(t).COM = meanPos(t,:);
    end
end
goodTimes = find(goodTimes);

% run linker

% set expected tags
dataProperties.MAXSPOTS = 2;
% set linker properties
dataProperties.linkerProperties.maxDistance = 50;
idlist = linker(slist, dataProperties,0);
% add randomState
idlist(1).info.randomState = slist(1).randomState;
% add labelcolor
labelcolor = [repmat('tag_',[idlist(1).stats.maxColor,1]),num2str([1:idlist(1).stats.maxColor]')]; %string tag_  1 - tag_999
labelcolor = cellstr(labelcolor); % conversion here, because matlab7 regexprep does not work on string arrays
labelcolor = regexprep(labelcolor,' ',''); %labelcolor = tag_1 - tag_999
idlist(1).stats.labelcolor = labelcolor;

if 0
    % show linking
    inputDataProperties.movieSize = [35,35,24,225];
                movieNumber = [num2str(projectNumber), '_A1'];
    [rawMovie] = ...
        generateProject(movieNumber, truePositions,inputDataProperties);
    LG_loadAllFromOutside(rawMovie,'',[],dataProperties,idlist,'idlist');
end

% make sure that the spot-assignment is correct. Get positions of first
% time and compare with positions from truePositions
% pixel 2 micron transformation
pixel2microns = [dataProperties.PIXELSIZE_XY * ones(1,2),...
    dataProperties.PIXELSIZE_Z];
ll1 = idlist(goodTimes(1)).linklist;
firstPos = ll1(:,[10,9,11]);
goodIdx = find(ll1(:,5)<2);
firstPos = firstPos(goodIdx,:);
firstPosTrue = squeeze(truePositions(goodTimes(1),1:3,:))';
firstPosTrue = firstPosTrue .* repmat(pixel2microns,nTags,1);

distance = distMat2(firstPos,firstPosTrue);

% LAP
nTags = length(goodIdx);
[pos2true, true2pos] = lap(distance, -1, 0, 1);

% reAssign. Create pdValue-matrix
pdValues = [goodIdx,ll1(goodIdx(pos2true(1:nTags)),4),zeros(size(goodIdx))];
[idlist, success, dataProperties] = ...
    LG_reAssignIdlist(idlist,goodTimes(1),...
    goodTimes, pdValues, 1, 0,dataProperties);
%

% debug = 0;
% 
% % read constants etc from input
% sp = size(truePositions);
% nTimepoints = sp(1);
% if size(sp)<3
%     nTags = 1;
% else
% nTags = sp(3);
% end
% 
% % pixel 2 micron transformation
% pixel2microns = [dataProperties.PIXELSIZE_XY * ones(1,2),...
%     dataProperties.PIXELSIZE_Z];
% 
% 
% idlist(1:nTimepoints) = struct('linklist',zeros(nTags,12),...
%     'centroid',zeros(1,3),'stats',[],...
%     'info',struct('detectQ_Pix',[],'trackQ_Pix',[]));
% 
% labelcolor = [repmat('tag_',[nTags,1]),num2str([1:nTags]')]; %string tag_  1 - tag_999
% labelcolor = cellstr(labelcolor); % conversion here, because matlab7 regexprep does not work on string arrays
% labelcolor = regexprep(labelcolor,' ',''); %labelcolor = tag_1 - tag_999
% idlist(1).stats.labelcolor = labelcolor;
% idlist(1).stats.labellist = labelcolor;
% idlist(1).weight = repmat(0.5,nTimepoints,1);
% idlist(1).stats.status = {};
% idlist(1).stats.maxColor = 2^(nTags+1);
% idlist(1).stats.created = date;
% idlist(1).info.randomState = slist(1).randomState;
% 
% [positions,sigmaPositions] = deal(repmat(NaN,[nTimepoints,4,nTags]));
% nsp = zeros(nTimepoints,1);
% goodTimeCt = 0;
% 
% 
% for t = 1:nTimepoints
%     
%     %==========
%     % DEBUG
%     %==========
%     if debug
%         if t==179
%             1;
%         end
%     end
%     % check that there are spots at all
%     if ~isempty(slist(t).sp)
% 
%         % Find nuimber of spots. If there are more spots than tags, we only
%         % consider the best nTag spots for linking.
%         nsp(t) = length(slist(t).sp);
% 
%         % initialize temp data. Already keep space for sigmaAmp.
%         % tmpPos, tmpSigmaPos, tmpChi2 can be longer than nTags
%         [tmpPos, tmpQ] = deal(zeros(nsp(t),4));
%         tmpChi2 = zeros(nsp(t),1);
% 
%         % set time
%         idlist(t).linklist(:,1) = t;
%         % tag color
%         idlist(t).linklist(:,4) = 2.^[0:nsp(t)-1]';
%         
% 
% 
%         % read the positions; then assign them to true
%         % positions with LAP. If less spots than true number of
%         % spots, we will fuse towards the "closest" spot (the
%         % one with the lowest cost)
%         tmpTruePos = permute(truePositions(t,:,:),[3,2,1]);
% 
%         tmpChi2 = slist(t).statistics.chi;
%         for i=1:nsp(t)
%             tmpPos(i,[2,1,3,4]) = [slist(t).sp(i).cord,...
%                 slist(t).sp(i).amp];
%             tmpQ(i,[2,1,3]) = ...
%                 diag(slist(t).statistics.Q(...
%                 (i-1)*3+1:i*3,(i-1)*3+1:i*3));
%             % add sigmaAmp
%         end
% 
%         % LAP. Cost = x^2+y^2+z^2+c*a^2
%         % constant is bigger the bigger the ratio between the
%         % true amplitudes. Use for now std(amp)/0.8
% 
%         idxR = permute(reshape(repmat([1:(nsp(t)*4)],nTags,1),...
%             nTags,nsp(t),4),[2,1,3]);
%         idxT = reshape(repmat([1:(nTags*4)],nsp(t),1),...
%             nsp(t),nTags,4);
%         weightMat = repmat(reshape([1,1,1,std(tmpTruePos(:,4))],...
%             1,1,4), nsp(t), nTags);
% 
%         costBlock = sum((...
%             tmpTruePos(idxT) - tmpPos(idxR)).^2 ...
%             .* weightMat, 3);
%         if nsp(t) == nTags
%             cost = costBlock;
%         else
%         cost = [costBlock, 1000 * eye(nsp(t)) - ones(nsp(t));...
%             1000 * eye(nTags) - ones(nsp(t)), zeros(nTags,nsp(t))];
%         end
% 
%         [r2t,t2r] = lap(cost);
% 
%         for tagNum=1:nTags
%             % time and color has been set already
% 
%             if t2r(tagNum) > nsp(t)
%                 % fusion tag. Fuse to closest spot by
%                 % searching the costBlock
%                 [dummy,spotNum] = min(costBlock(:,tagNum));
% 
%                 % amplitude is theoretical amp for now
%                 idlist(t).linklist(tagNum,8) = tmpTruePos(tagNum,4);
% 
%                 % flag
%                 idlist(t).linklist(tagNum,5) = 1;
% 
% 
%             else
% 
%                 % spot number is t2r(i), while tagNum is the tagNumber
%                 spotNum = t2r(tagNum);
% 
%                 % amplitude is the spot amplitude
%                 idlist(t).linklist(tagNum,8) = ...
%                     tmpPos(spotNum,4);
% 
%                 % write already spotColor (might be changed below)
%                 idlist(t).linklist(tagNum,3) = ...
%                     idlist(t).linklist(tagNum,4);
% 
%                 % only write sigmaPosition, position if the tag was
%                 % actually found
%                 % write distances
%                 positions(t,:,tagNum) = tmpPos(spotNum,:);
%                 sigmaPositions(t,:,tagNum) = sqrt(...
%                     tmpQ(spotNum, :) * tmpChi2(spotNum));
% 
% 
%             end
% 
%             % position. Rearrange.
%             idlist(t).linklist(tagNum,[10,9,11]) = ...
%                 tmpPos(spotNum,1:3) .* pixel2microns;
% 
%             % spot number
%             idlist(t).linklist(tagNum,2) = spotNum;
% 
%             % Q-matrix. Augment
%             idlist(t).info.detectQ_Pix = blkdiag(...
%                 idlist(t).info.detectQ_Pix,...
%                 diag(tmpQ(spotNum, [2,1,3])));
%             idlist(t).info.trackQ_Pix = [];
%             
%             % chi2
%             idlist(t).linklist(tagNum,12) = tmpChi2(spotNum);
% 
%         end
% 
%         % take care of fusions
%         if nsp(t) < nTags
%             % loop through spots, assign the right colors
%             for i = 1:max(idlist(t).linklist(:,2))
%                 if idlist(t).linklist(i,3) == 0
%                     spotIdx = find (idlist(t).linklist(i,2)==...
%                         idlist(t).linklist(:,2));
%                     idlist(t).linklist(spotIdx,3) = ...
%                         sum(idlist(t).linklist(spotIdx,4));
%                 end
%             end
%         end
% 
% 
% 
%         % idlist: linkup/linkdown
%         if goodTimeCt > 0
%             idlist(t).linklist(:,6) = idlist(goodTimeCt).linklist(:,3);
%             idlist(goodTimeCt).linklist(:,7) = idlist(t).linklist(:,3);
%         end
%         % centroid
%         idlist(t).centroid = mean(idlist(t).linklist(:,9:11),1);
%         goodTimeCt = t;
% 
%     else % if no spot found
%         idlist(t).linklist = [];
%     end % if ~isempty
% 
% end
% 
% 
