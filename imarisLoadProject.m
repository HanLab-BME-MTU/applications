function imarisApplication = imarisLoadProject(addOrReplace,centerOrNot)
%IMARISLOADPROJECT loads, transforms and displays chromdyn-data in imaris
%
% SYNOPSIS imarisApplication = imarisLoadProject(addOrReplace,centerOrNot)
%
% INPUT    the code asks for project data and a filtered movie
%          addOrReplace : whether to add to or replace the current surpass 
%                           scene 'add'/{'replace'}
%          center       : whether to center the data on the centroid [{0}/1]
%
% OUTPUT   imarisApplication: Handle to imaris application
%          (the code will launch imaris and show the image data there)
%
% c: 04/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============
% TEST INPUT
%=============
replace = 1;
center  = 0;
% find handle to imaris figure if exists
imFigH = findall(0,'Name','Imaris figure');
if isempty(imFigH)
    imFigHandles = [];
else
    % get handles
    imFigHandles = guidata(imFigH);
end


if nargin < 1 || isempty(addOrReplace)
    % default
else
    switch addOrReplace
        case 'add'
            if ~isempty(imFigH)
                replace = 0;
            else
                warning('No imaris figure found to add to - launching new imaris')
            end
        case 'replace'
            replace = 1;
        otherwise
            error('please specify either ''add'' or ''replace''!')
    end
end

if nargin < 2 || isempty(centerOrNot)
    % default
else
    switch centerOrNot
        case 1
            center = 1;
        case 0
            center = 0;
        otherwise
            error('Please specify 1 or 0 for the second input variable centerOrNot')
    end
end
%========================
% LOAD DATA
%========================

[idlist,dataProperties,projectProperties,dummy,filteredMovie] = loadProjectData;

projectName = dataProperties.name;

%============================
% data loaded:
%   - idlist
%   - dataProperties
%   - filteredMovie
%   - projectName

 

%===========================
% calculate imaris data
%===========================

% loop throught timepoints. Collect all spots for display, and create a
% track for every tag.

tMax = length(idlist);
oneVoxel = [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z];

spotCt = 0;
edgeCt = 0; 
connectionCt = 0;
nTags  = length(idlist(1).stats.labelcolor);


for t = 1:tMax
    if ~isempty(idlist(t).linklist)
        
        % find list of spots
        spotNumVector = idlist(t).linklist(:,2);
        [spotNum, singleIdx] = unique(spotNumVector);
        spotCtNew = spotCt + spotNum(end);
        
        
        
        % fill spotList : [tp, Xspot, Yspot, Zspot, Amp]
        % careful! we have to correct for the fact that a matlab matrix
        % starts with [1,1,1] -> the um-coordinates are shifted by 1 voxel
        for iSpot = spotNum'
            spotList(spotCt+iSpot,:) = [t,...
                    idlist(t).linklist(singleIdx(iSpot),9:11)-oneVoxel,...
                    sum(idlist(t).linklist(find(spotNumVector==iSpot),8))];
        end
        
        % center if selected: subtract the mean of the spots
        if center
            spotList(spotCt+1:spotCtNew,2:4) = spotList(spotCt+1:spotCtNew,2:4) - ...
                repmat(mean(spotList(spotCt+1:spotCtNew,2:4),1),spotNum(end),1);
        end

        % fill edgeList - [prevSpot,currentSpot], as ntpx2xntag array
        
        linkup = idlist(t).linklist(:,6);
        % we do not want to write linkup for the very first linklist - but
        % we do have to make sure that we are not victim of a sloppy update
        % of the linklist after deleting the first frame in a movie
        if all(linkup) & spotCt ~= 0 
            edgeCt = edgeCt +1;
            edgeList(edgeCt,:,:) = [reshape(linkup + spotCtOld, [1,1,nTags]),...
                    reshape(spotNumVector + spotCt, [1,1,nTags])];
        end
        
        % fill connectionList. Make connections from tag 4 to tags 1-3,
        % from tag 3 to 1-2 and from 2 to 1
        connCt = 0;
        for iTag = nTags-1:-1:1
            tagNum1 = iTag + 1;
            for tagNum2 = 1:iTag
                connectionCt = connectionCt + 1;
                connectionList(connectionCt,:,connCt + tagNum1-tagNum2) = ...
                    [spotNumVector(tagNum1) spotNumVector(tagNum2)] + spotCt;
            end
            connCt = connCt + tagNum1;
        end
        
        % count total # of spots including this frame
        spotCtOld = spotCt;
        spotCt = spotCtNew;
        
       
    end
end

if center
% create pseudo-movie if centered
minMovieSizeUM = max(spotList(:,2:4),[],1)-min(spotList(:,2:4),[],1);

% increase movie size by 10%
movieSizeUM = minMovieSizeUM * 1.1;

% calc in pixels
movieSizePIX = ceil(movieSizeUM./oneVoxel);
% one channel, ntp time points
filteredMovie = zeros([movieSizePIX,1,size(filteredMovie,5)]);

% shift zero of coordinates to center of movie
spotList(:,2:4) = spotList(:,2:4) + repmat(movieSizeUM/2, spotCt, 1);
end

%================================

%===================================
% Launch imaris and feed with data
%====================================

% get handle to imaris application. Start new if necessary
if isempty(imFigH)
    imaApplication = actxserver('Imaris.Application');
else
    imaApplication = imFigHandles.imarisData.imaApplication;
end


%--------- load movie into data 
if center
    imarisMovieSize = movieSizeUM;

else
    matlabMovieSize = size(filteredMovie);


    % convert every single frame to [0 1]
    for t = 1:matlabMovieSize(5)
        currImg = filteredMovie(:,:,:,:,t);
        currImg = currImg - min(currImg(:));
        filteredMovie(:,:,:,:,t) = currImg/max(currImg(:));
    end

    % switch x,y
    filteredMovie = permute(filteredMovie,[2,1,3,4,5]);

    % calculate extent; do not forget to subtract 1!
    imarisMovieSize = ...
        (matlabMovieSize(1:3)-1) .* oneVoxel;

end

% create dataSet
imaDataSet = imaApplication.mFactory.CreateDataSet;

% convert double to single and set
imaDataSet.SetData(single(filteredMovie));

% set extent of movie 

imaDataSet.mExtendMaxX = imarisMovieSize(1);
imaDataSet.mExtendMaxY = imarisMovieSize(2);
imaDataSet.mExtendMaxZ = imarisMovieSize(3);

% load dataSet into Imaris
imaApplication.mDataSet = imaDataSet;


%-------- create surpass scene
if replace
    imaSurpassScene = imaApplication.mFactory.CreateDataContainer;
    imaApplication.mSurpassScene = imaSurpassScene;
else
    imaSurpassScene = imaApplication.mSurpassScene;
end

% add light, frame if necessary
if replace
    imaLightSource = imaApplication.mFactory.CreateLightSource;
    imaSurpassScene.AddChild(imaLightSource);
    imaFrame = imaApplication.mFactory.CreateFrame;
    imaSurpassScene.AddChild(imaFrame);
end
    
    
%-------- load spots ---------
% create a spots component
imaSpots = imaApplication.mFactory.CreateSpots;
% name the spots after the project
imaSpots.mName = projectName;

% put them in - don't forget the -1
%([x,y,z], t, radius)
imaSpots.Set(spotList(:,2:4), spotList(:,1)-1, spotList(:,5)*100);
% add spots to the surpass scene as last child
imaSurpassScene.AddChild(imaSpots);



% sequentially load the different tags, name them, and set up the tracks
for i = 1:nTags
    
    % create track component
    imaTrack = imaApplication.mFactory.CreateTrack;
    
    % name it according to the tag
    imaTrack.mName = idlist(1).stats.labelcolor{i};
    
    % set spots possible for this track (=all of them)
    imaTrack.SetSpots(imaSpots)
    
    % set edges. Rember -1
imaTrack.SetEdges(edgeList(:,:,i)-1)

% name timepoints with intensity (good idea?)
    for nsp = 1:size(spotList,1)
        imaTrack.GetChild(nsp-1).mName = ['t ' num2str(spotList(nsp,1)) ' int ' num2str(spotList(nsp,5))];
        imaTrack.GetChild(nsp-1).mVisible = 0;
    end
    
    % here is where I would love to set the line thickness
    
    
    % add track to the surpass scene as last child
    imaSurpassScene.AddChild(imaTrack);
end
    

% now load all the connections & name them
% Make connections from tag 4 to tags 1-3,
% from tag 3 to 1-2 and from 2 to 1
connCt = 0;
for iTag = nTags-1:-1:1
    tagNum1 = iTag + 1;
    for tagNum2 = 1:iTag
        
        % create new track component
        imaTrack = imaApplication.mFactory.CreateTrack;
        
        % name it
        imaTrack.mName = ['conn. ' idlist(1).stats.labelcolor{tagNum1} idlist(1).stats.labelcolor{tagNum2}];
        
        % and load the list. Don't forget -1
        imaTrack.SetSpots(imaSpots)
        imaTrack.SetEdges(connectionList(:,:,connCt + tagNum1-tagNum2)-1);
        
        % make all the spots invisible
        for nsp = 1:size(spotList,1)
            imaTrack.GetChild(nsp-1).mVisible = 0;
        end
        
        % finally, add connection to surpass scene
        imaSurpassScene.AddChild(imaTrack);
        
    end
    connCt = connCt + tagNum1;
end

% make sure that the view fits and set it to surpass
imaApplication.mSurpassCamera.Fit;
imaApplication.mViewer = 'eViewerSurpass';

%=================================
% LAUNCH "IMARIS-FIGURE"
%=================================

% launch figure
if replace
    imFigH = findall(0,'Name','Imaris figure');
    delete(imFigH)
    imFigH = figure('Name','Imaris figure','NumberTitle','off');
end

% get handles
imFigHandles = guidata(imFigH);

% store project data
imFigHandles.projectData.idlist = idlist;
imFigHandles.projectData.dataProperties = dataProperties;
%imFigHandles.projectData.projProperties = projProperties;

% store imaris data
imFigHandles.imarisData.imaApplication = imaApplication;
imFigHandles.imarisData.spotList = spotList;
imFigHandles.imarisData.edgeList = edgeList;

guidata(imFigH,imFigHandles)

if nargout > 0
    imarisApplication = imaApplication;
end
