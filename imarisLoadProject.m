function imarisLoadProject(addOrReplace,centerOrNot)
%imarisLoadProject loads, transforms and displays chromdyn-data in imaris
%
% SYNOPSIS imarisLoadProject
%
% INPUT    the code asks for project data and a filtered movie
%          addOrReplace : whether to add to or replace the current surpass 
%                           scene 
%          center       : whether to center the data on the centroid [{0}/1]
%
% OUTPUT   the code will launch imaris and show the image data there
%
% c: 04/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=============
% TEST INPUT
%=============
replace = 1;
center  = 0;

if nargin < 1 | isempty(addOrReplace)
    % default
else
    switch addOrReplace
        case 'add'
            % check whether there is an imarisFigure open
            imFigH = findall(0,'Name','Imaris figure');
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

if nargin < 1 || isempty(centerOrNot)
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

% load imaris
vImarisApplication = actxserver('Imaris.Application');


vImarisDataSet = vImarisApplication.mDataSet;

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

% convert double to single and set
vImarisDataSet.SetData(single(filteredMovie));



% set extent of movie 

vImarisDataSet.mExtendMaxX = imarisMovieSize(1);
vImarisDataSet.mExtendMaxY = imarisMovieSize(2);
vImarisDataSet.mExtendMaxZ = imarisMovieSize(3);


%-------- create surpass scene
if replace
    vImarisSurpassScene = actxserver('Imaris.DataContainer');
    vImarisApplication.mSurpassScene = vImarisSurpassScene;
else
    imFigHandles = guidata(imFigH);
    vImarisSurpassScene = imFigHandles.vImarisSurpassScene;
end

% add light, frame if necessary
if replace
    vImarisLightSource = actxserver('Imaris.LightSource');
    vImarisSurpassScene.AddChild(vImarisLightSource);
end
    vImarisFrame = actxserver('Imaris.Frame');
    vImarisSurpassScene.AddChild(vImarisFrame);
    
% if ~replace
%     % problem: if I add a second track the first frame gets weird -
%     % therefore add a second (good) frame and delete the first one
%     nKids = vImarisSurpassScene.GetNumberOfChildren;
%     kid = nKids-2;
%     frameDeleted = 0;
%     vImarisTypeConvert = actxserver('Imaris.TypeConvert');
%     while ~frameDeleted & kid >= 0
%         
%         k = vImarisSurpassScene.GetChild(kid);
%         if vImarisTypeConvert.IsFrame(k)
%             vImarisSurpassScene.RemoveChild(k);
%             frameDeleted = 1;
%         else
%             kid = kid-1;
%         end
%     end
% end
%-------- load spots ---------
% create a spots component
vImarisSpots = actxserver('Imaris.Spots');
% name the spots after the project
vImarisSpots.mName = projectName;
% set size to 3 pixels (intensity does not work yet)
vImarisSpots.mRadius = dataProperties.PIXELSIZE_XY*3;
% put them in - don't forget the -1
vImarisSpots.SetPos(spotList(:,2:4), spotList(:,1)-1);
% add spots to the surpass scene as last child
vImarisSurpassScene.AddChild(vImarisSpots);



% sequentially load the different tags, name them, and set up the tracks
for i = 1:nTags
    
    % create track component
    vImarisTrack = actxserver('Imaris.Track');
    
    % name it according to the tag
    vImarisTrack.mName = idlist(1).stats.labelcolor{i};
    
    % put in edges. Do not forget -1
    vImarisTrack.SetTrack(vImarisSpots, edgeList(:,:,i)-1);
    
    for nsp = 1:size(spotList,1)
        vImarisTrack.GetChild(nsp-1).mName = ['t ' num2str(spotList(nsp,1)) ' int ' num2str(spotList(nsp,5))];
        vImarisTrack.GetChild(nsp-1).mVisible = 0;
    end
    
    % add track to the surpass scene as last child
    vImarisSurpassScene.AddChild(vImarisTrack);
    
end

% now load all the connections & name them
% Make connections from tag 4 to tags 1-3,
% from tag 3 to 1-2 and from 2 to 1
connCt = 0;
for iTag = nTags-1:-1:1
    tagNum1 = iTag + 1;
    for tagNum2 = 1:iTag
        
        % create new track component
        vImarisTrack = actxserver('Imaris.Track');
        
        % name it
        vImarisTrack.mName = ['conn. ' idlist(1).stats.labelcolor{tagNum1} idlist(1).stats.labelcolor{tagNum2}];
        
        % and load the list. Don't forget -1
        vImarisTrack.SetTrack(vImarisSpots, connectionList(:,:,connCt + tagNum1-tagNum2)-1);
        
        % make all the spots invisible
        for nsp = 1:size(spotList,1)
            vImarisTrack.GetChild(nsp-1).mVisible = 0;
        end
        
        % finally, add connection to surpass scene
        vImarisSurpassScene.AddChild(vImarisTrack);
        
    end
    connCt = connCt + tagNum1;
end

vImarisApplication.mViewer = 'eViewerSurpass';

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
imFigHandles.imarisData.spotList = spotList;
imFigHandles.imarisData.edgeList = edgeList;

imFigHandles.vImarisSurpassScene = vImarisSurpassScene;

guidata(imFigH,imFigHandles)
