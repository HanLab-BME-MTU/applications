function [data,orientation,positions,sigmaZero,dataProperties,snrMax,isTracked,rayleighLimit,fusedTag] = calculateTrajectoryFromIdlist(idlist,dataProperties,tag1,tag2,opt)
%CALCULATETRAJECTORYFORMIDLIST calculates distance trajectories from an idlist
%
%SYNOPSIS [data, orientation, positions, sigmaZero,dataProperties,...
%           snrMax,isTracked,rayleighLimit] = ...
%           calculateTrajectoryFromIdlist(idlist,dataProperties,tag1,tag2,opt)
%
%INPUT    idlist         : any type of idlist
%         dataProperties : the corresponding data properties
%         tag1,2         : number (1:nTag) of first and second tag, OR string from labelcolor
%   if the first four arguments are omitted, the user is asked for input
%         opt            : optional options structure
%                           .info :   struct that should be returned in the
%                                     output field info.
%                           .calc2d : [{0}/1/2] depending on whether the
%                                     normal 3D-data or the maxProjection
%                                     or the in-focus-slice data should be
%                                     used.
%                           .nanList: [{0}/1] whether to give data in the
%                                     form of nanList (converted with convertTrajectoryData)
%                           .oldIdlist: if specified, the program
%                                       calculates snrMax
%                           .realTime : [0/{1}] whether to use real time or
%                                       rounded time
%
%OUTPUT   data           : structure containing trajectory data
%                           .time       = [time in sec, sigmaTime] -
%                           sigmaTime is 1/4 the time between the first and
%                           the last slice in a frame
%                           .timePoints = corresponding timepoint in raw data movie
%                           .distance   = [distance in um (tag1-tag2), distanceSigma]
%                           .info: substructure with fields helping to
%                                  describe the data
%                               .tags   cell containing the two tags used
%                                       calculate the trajectory {'tag1';'tag2'}
%
%         orientation    : 1-by-n matrix of angles [rad] from the xy-plane
%
%         positions      : 1-by-2 struct with fields
%           .coordinates : ntp-by-3 array of coordinates. Deleted frames
%                           are replaced by NaNs
%           .covariances : 3-by-3-by-ntp array of covariances
%
%         sigmaZero      : average chi2 (noise) of the two tags combined
%
%         dataProperties : dataProperties
%
%         snrMax         : SNR of the higher intensity tag (not supported
%                           for idlist2 yet)
%
%         isTracked      : whether the frame has been tracked or not
%
%         rayleighLimit  : Rayleigh limit calculated from the orientation
%
%         fusedTag       : Whether any of the two tags is fused to another
%                          tag. 
%
%
%REMARKS  fusions will not be considered as valid timepoints
%         Warning: the length of fusedTag etc will be equal to the length
%         of timepoints 
%
%c: 11/03 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================
%% ------TEST INPUT--------
%=========================

%=========================
%check nargin (we do not really check for correct idlist, dataProperties)
if nargin < 2 || isempty(idlist) || isempty(dataProperties)
    [idlist,dataProperties] = loadProjectData;
end

if nargin < 4 || isempty(tag1) || isempty(tag2)

    % get tags
    tagList = idlist(1).stats.labelcolor;

    h = helpdlg('Please choose first tag','');
    uiwait(h);
    tagNum1 = chooseFileGUI(tagList);
    tag1 = tagList{tagNum1};
    tagList(tagNum1)=[];

    h = helpdlg('Please choose second tag','');
    uiwait(h);
    tagNum2 = chooseFileGUI(tagList);
    tag2 = tagList{tagNum2};

    clear tagNum1 tagNum2 tagList
end
%=========================

%check options structure
%set defaults
info = [];
calc2d = 0;
nanList = 0;
oldIdlist = [];
realTime = 1; % if zero, time is rounded

%go through fields of opt, change defaults with input values if they exist
if nargin < 5 || isempty(opt)
    % keep defaults
else
    if isfield(opt,'info')
        info = opt.info;
    end
    if isfield(opt,'calc2d')
        calc2d = opt.calc2d;
    end
    if isfield(opt,'nanList')
        nanList = opt.nanList;
    end
    if isfield(opt,'oldIdlist') && nargout > 5
        oldIdlist = opt.oldIdlist;
    end
    if isfield(opt,'realTime')
        realTime = opt.realTime;
    end
end

%check tags and convert to numbers
try
    labelcolor = idlist(1).stats.labelcolor;
    if ischar(tag1)
        tag1 = find(strcmp(tag1,labelcolor));
    end
    if ischar(tag2)
        tag2 = find(strcmp(tag2,labelcolor));
    end
    if isempty(tag1)||isempty(tag2)||any([tag1;tag2]>length(labelcolor))
        error('no tags found')
    else
        % read 'tags'
        tags = labelcolor([tag1,tag2]);
    end
    %define indexList for indexing Q-matrix
    tagIdxList1 = [(tag1-1)*3+1:tag1*3];
    tagIdxList2 = [(tag2-1)*3+1:tag2*3];
catch
    error(['bad idlist or pointer to nonexistent tag: ',lasterr])
end



%calc pix2mu-matrix
pix2muMat = diag([dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_Z]);

%=========================
%---END TEST INPUT--------
%=========================

%=========================
%% CALC 2D-CASE IF SELECTED
%=========================
switch calc2d
    case 1 %maximum projection
        [outputData3S,outputDataMP] = threeVsTwoD(idlist, dataProperties);
        idlist = outputDataMP.idlist;
    case 2 %in-focus slices
        [outputData3S,outputDataMP] = threeVsTwoD(idlist, dataProperties);
        idlist = outputData3S.idlist;
    otherwise
        %keep normal idlist
end

if isempty(idlist)
    error('can not calculate trajectory: idlist is empty')
end

%=========================
%---END CALC 2D-CASE IF SELECTED
%=========================


%===================================
%% CALCULATE PARAMETERS FROM IDLIST
%===================================

if checkIdlist(idlist,1)
    % calculate data from new idlist
    [distance, distanceVectors, dummy, dummy, dummy, idxLists] = ...
        idlist2distMat(idlist,dataProperties);
    
    % remove any fusions to each other, and estimated tags from the
    % list
    
    % isGoodTime has ones wherever there is a good time
    isGoodTime = idxLists.isGoodTime;
    
    % estimatedIdx has ones wherever there is an estimated tag
    estimatedIdx = ...
        idxLists.estimatedTag(:, tag1) | idxLists.estimatedTag(:, tag2);
    
    % ownFusion has ones wherever the tags are fused to each other
    ownFusion = idxLists.fusedTag(:,tag1) == tag2 & ...
        idxLists.fusedTag(:,tag2) == tag1;
    fusedTag1 = idxLists.fusedTag(:,tag1) ~= 0;
    fusedTag2 = idxLists.fusedTag(:,tag2) ~= 0;
    
    % create list of good timepoints by removing timepoints with deleted
    % frames, with one or both tags estimated, or both tags fused one
    % another
    timePoints = find(isGoodTime & ~estimatedIdx & ~ownFusion);
    
    % correct fusedTag
    fusedTag = fusedTag1(timePoints) + 2*fusedTag2(timePoints);
    
    
    % now read distance into old variables. Careful with the order of tags
    tagRow = max(tag1,tag2);
    tagCol = min(tag1,tag2);
    distanceN = squeeze(distance(tagRow, tagCol, timePoints));
    distanceSigma = squeeze(distance(tagCol, tagRow, timePoints));
    normedVectors = squeeze(distanceVectors(tagRow, tagCol, timePoints, :));
    
    % snrMax is not returned by idlist2distMat yet. I'll add that should I
    % ever need it.
    snrMax = [];
    positions = []; % you can get all the relevant info from idlist2distMat
    sigmaZero = [];
    
    % isTracked will be 1 if any of the two tags has been found by the
    % tracker
    isTracked = idxLists.isTracked(:,tag1) | idxLists.isTracked(:,tag2);
   
else
    % calculate from old idlist
    [distanceN, distanceSigma, normedVectors, timePoints, snrMax, isTracked, fusedTag] = ...
        calculateFromOldIdlist;
end
    






%=========================
%% ------READ TIME AND TIMESIGMA
%=========================

if realTime
    timeAll = dataProperties.frameTime;
    time = mean(timeAll(timePoints,:),2);
    %timeSigma is one fourth the time between (lastCol - firstCol of frameTime)
    timeSigma = (timeAll(timePoints,end)-timeAll(timePoints,1))/4;

else
    % start at t=1*rounded timeLapse
    avgTimeLapse = round(dataProperties.timeLapse);
    time = (timePoints) * avgTimeLapse;
    timeSigma = zeros(size(time));

    if abs(dataProperties.timeLapse - avgTimeLapse) > 0.11
        warning('CALCULATETRAJECTORYFROMIDLIST:wrongTime',...
            'Possible source of error: timeLapse %f used instead of %f\n in %s',avgTimeLapse,dataProperties.timeLapse,dataProperties.name);
    end
end


%=========================
%---END READ TIME AND TIMESIGMA
%=========================


%=========================
%% ORIENTATION & RL
%=========================
if nargout > 1
    [rayleighLimit, orientation] =...
        rayleighFromOri(normedVectors, ...
        dataProperties.WVL, dataProperties.NA, []);
else
    orientation = [];
end

%=========================
%% -------ASSIGN OUTPUT
%=========================
data.distance   = [distanceN,distanceSigma];
data.time       = [time,timeSigma];
data.timePoints = timePoints;
data.info       = info;
data.info.tags  = tags;

if nanList && length(timePoints)>0
    data = convertTrajectoryData(data);
end


%---END ASSIGN OUTPUT


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin nested function
function [distanceN, distanceSigma, normedVectors, timePoints, snrMax, isTracked, fusedTag] = calculateFromOldIdlist


%=========================
%% READ COORDS AND SIGMAZERO AND Q
%=========================

%cat linklist and extract coords, chi2, timePoints
linkLists = cat(3,idlist.linklist);


% coords
coords = linkLists([tag1,tag2],9:11,:);



% check for fused tags. Distinguish which tag is a fusion tag, and remember
% whether they are fused to each other
nTags = size(linkLists,1);
for t = size(coords,3):-1:1
    removeFusedTag(t,1) = isequal(coords(1,:,t),coords(2,:,t));
    
    fusedTag1(t,1) = any(any(repmat(coords(1,:,t),nTags,1) == linkLists(:,9:11,t),2));
    fusedTag2(t,1) = any(any(repmat(coords(2,:,t),nTags,1) == linkLists(:,9:11,t),2));
    
end

% remove entries where tags are fused to each other from coords and
% linkLists. FusedTags is 0,1,2,3 depending on whether there is no fusion,
% fusion of tag 1 or fusion of tag 2, or fusion of both (to different
% tags!)
linkLists(:,:,removeFusedTag) = [];
coords(:,:,removeFusedTag) = [];
fusedTag1(removeFusedTag) = [];
fusedTag2(removeFusedTag) = [];
fusedTag = fusedTag1 + 2*fusedTag2;

% timepoints
timePoints = squeeze(linkLists(1,1,:));
maxTime = max(timePoints);


% sigmaZero (be backward compatible)
if size(linkLists,2)<12
    fillSigma = 1;
    sigma0 = zeros(length(timePoints),2);
else
    fillSigma = 0;
    sigma0 = squeeze(linkLists([tag1,tag2],12,:))';
end


% read Q-Matrices, isTracked
[covariance1,covariance2] = deal(repmat(NaN,[3,3,maxTime]));
isT = repmat(NaN,maxTime,1);
for t = timePoints'
    if isfield(idlist(t).info,'trackQ_Pix') && ~isempty(idlist(t).info.trackQ_Pix) && isfield(idlist(t).info,'trackerMessage')
        % tracked frame. Add covariance of source and target
        sourceT = str2double(idlist(t).info.trackerMessage.source);
        covariance1(:,:,t) = pix2muMat*(idlist(t).info.trackQ_Pix(tagIdxList1,tagIdxList1)+...
            idlist(sourceT).info.detectQ_Pix(tagIdxList1,tagIdxList1))*pix2muMat;
        covariance2(:,:,t) = pix2muMat*(idlist(t).info.trackQ_Pix(tagIdxList2,tagIdxList2)+...
            idlist(sourceT).info.detectQ_Pix(tagIdxList2,tagIdxList2))*pix2muMat;
        isT(t) = 1;
    else
        % only detected
        covariance1(:,:,t) = pix2muMat*(idlist(t).info.detectQ_Pix(tagIdxList1,tagIdxList1))*pix2muMat;
        covariance2(:,:,t) = pix2muMat*(idlist(t).info.detectQ_Pix(tagIdxList2,tagIdxList2))*pix2muMat;
        isT(t) = 0;
    end



    if fillSigma
        sigma0(t,:) = idlist(t).info.noise([tag1,tag2]);
    end
    
    % check whether the tag is fused
    
end

% fill coords and sigmaZero with NaNs for deleted times
[pos1,pos2] = deal(repmat(NaN,maxTime,3));
sigmaZero   = pos1(:,1);

pos1(timePoints,:)    = squeeze(coords(1,:,:))';
pos2(timePoints,:)    = squeeze(coords(2,:,:))';
if fillSigma
    sigmaZero(timePoints) = mean(sigma0(timePoints,:),2);
else
    sigmaZero(timePoints) = mean(sigma0,2);
end

% take care of nans in isTracked
isTracked = isT(timePoints);

% calculate snrMax if selected
if ~isempty(oldIdlist)
    % loop to read amplitudes
    intData = repmat(NaN,maxTime,2);
    for t = timePoints' % everything that's in lastResult is also in idlist_L
        rowIdx1 = find(oldIdlist(t).linklist(tag1,2)==oldIdlist(t).linklist(:,2));
        rowIdx2 = find(oldIdlist(t).linklist(tag2,2)==oldIdlist(t).linklist(:,2));
        % intData = [totalAmp1, totalAmp2]
        intData(t,:) = ...
            [sum(oldIdlist(t).linklist(rowIdx1,8)),sum(oldIdlist(t).linklist(rowIdx2,8))];
    end
    % decide for a tag
    meanInt = nanmean(intData);
    higherInt = 1 + (meanInt(1) < meanInt(2));

    snrMax = intData(timePoints,higherInt)./sqrt(sigma0(:,higherInt));

else
    snrMax = [];
end


%=========================
%---END READ COORDS AND SIGMAZERO AND Q
%=========================



%=========================
%% CALC DISTANCE & DISTANCESIGMA
%=========================

% prepare data for distance calculation
positions(1:2) = struct('coordinates',{pos1,pos2},'covariances',{covariance1,covariance2});

[distanceN, distanceSigma, normedVectors] = deltaCoordinates(positions,sigmaZero);

% throw away data with zero distance or NaN
badIdx = find(distanceN == 0 | isnan(distanceN));

distanceN(badIdx) = [];
distanceSigma(badIdx) = [];
normedVectors(badIdx,:) = [];

% if no zero distance, all badIdx are not in timePoints
badTpIdx = (ismember(timePoints,badIdx));
timePoints(badTpIdx) = [];


%=========================
%---END CALC DISTANCE AND DISTANCESIGMA
%=========================

end % subfunction

end % main function