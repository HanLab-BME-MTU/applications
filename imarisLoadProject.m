function imarisLoadProject
%imarisLoadProject loads, transforms and displays chromdyn-data in imaris
%
% SYNOPSIS imarisLoadProject
%
% INPUT    the code asks for project data and a filtered movie
%
% OUTPUT   the code will launch imaris and show the image data there
%
% c: 04/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%========================
% LOAD DATA - make seperate function, replace in label_LoadSlist
%========================

%check if default biodata-dir exists and cd if exist
mainDir = cdBiodata(2);

%get project data file
[fileName,pathName] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files'},'select project data file');

if fileName==0;
    % not defined yet
    error('no data loaded')
    return
else
    cd(pathName);
    data = load(fileName); %loads everything into the structure data
    if ~isfield(data,'dataProperties')
        error('No dataProperties in project data: corrupt data file');
    else
        dataProperties = data.dataProperties;
    end
    
    %load projectProperties
    if ~isfield(data,'projProperties')
        error('No projProperties in project data!');
        uiwait(h);
        return %end evaluation here
    else
        projProperties = data.projProperties;
    end
    
    
    
    %--------------try to load filtered movie
    %try to find filenames in the path from which projectData has been loaded
    filteredMovieName = chooseFile('filtered_movie',[],'new');
    altFilteredMovieName = chooseFile('moviedat',[],'new');
    if isempty(filteredMovieName)
        if isempty(altFilteredMovieName) %to ensure compatibility with earlier versions
            disp('no filtered movie found. load unfiltered movie instead')
            if findstr(projProperties.dataPath(end-10:end),'crop')|findstr(projProperties.dataPath(end-10:end),'corr')
                %cropped movie
                moviename = chooseFile('.r3c');
                filteredMovie  =  readmat(moviename);
            else
                %normal movie
                moviename = chooseFile('.r3d');
                filteredMovie  =  r3dread(moviename);
            end
        else
            filteredMovie = readmat(altFilteredMovieName);
        end
    else 
        filteredMovie = readmat(filteredMovieName);
    end;
    
    %test if everything correctly loaded
    if ~exist('filteredMovie','var')
        error('no movie found')
        return
    end
    %---let the user choose which idlist to load
    
    %find which idlists there are
    dataFieldNames = fieldnames(data);
    idnameListIdx = strmatch('idlist',dataFieldNames);
    idnameList = dataFieldNames(idnameListIdx);
    
    %have the user choose, if there is more than one entry left
    switch length(idnameList)
        case 0 %no idlist loaded. continue w/o loading
            
            idname = '[]';
            error('no idlist found in data')
            
        case 1 %only one idlist loaded. Continue
            
            idname = char(idnameList);
            
        otherwise %let the user choose
            idSelect = chooseFileGUI(idnameList);
            if isempty(idSelect)
                idname = '';
            else
                idname = idnameList{idSelect};
            end
    end
end

if isempty(idname)
    error('no idlist loaded')
else
    idlist = eval(['data.' idname ';']);
end

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
        
        
        
        % fill spotList : [tp, Xspot, Yspot, Zspot, Amp]
        % careful! we have to correct for the fact that a matlab matrix
        % starts with [1,1,1] -> the um-coordinates are shifted by 1 voxel
        for iSpot = spotNum'
            spotList(spotCt+iSpot,:) = [t,...
                    idlist(t).linklist(singleIdx(iSpot),9:11)-oneVoxel,...
                    sum(idlist(t).linklist(find(spotNumVector==iSpot),8))];
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
        spotCt = spotCt + spotNum(end);
        
       
    end
end


%================================

%===================================
% Launch imaris and feed with data
%====================================

% load imaris
vImarisApplication = actxserver('Imaris.Application');



% connect to Imaris - done at the beginning of the program

vImarisDataSet = vImarisApplication.mDataSet;

%--------- load movie into data 
matlabMovieSize = size(filteredMovie);


% convert every single frame to [0 1]
for t = 1:matlabMovieSize(5)
    currImg = filteredMovie(:,:,:,:,t);
    currImg = currImg - min(currImg(:));
    filteredMovie(:,:,:,:,t) = currImg/max(currImg(:));
end

% convert double to single and switch x,y
filteredMovie = permute(filteredMovie,[2,1,3,4,5]);
vImarisDataSet.SetData(single(filteredMovie));



% set extent of movie do not forget to subtract 1!
imarisMovieSize = ...
    [matlabMovieSize(1:3)-1] .* [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z];

vImarisDataSet.mExtendMaxX = imarisMovieSize(1);
vImarisDataSet.mExtendMaxY = imarisMovieSize(2);
vImarisDataSet.mExtendMaxZ = imarisMovieSize(3);


%-------- create surpass scene
vImarisSurpassScene = actxserver('Imaris.DataContainer');
vImarisApplication.mSurpassScene = vImarisSurpassScene;

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
imFigH = figure('Name','Imaris figure','NumberTitle','off');

% get handles
imFigHandles = guidata(imFigH);

% store project data
imFigHandles.projectData.idlist = idlist;
imFigHandles.projectData.dataProperties = dataProperties;
%imFigHandles.projectData.projProperties = projProperties;

% store imaris data
imFigHandles.imarisData.spotList = spotList;
imFigHandles.imarisData.edgeList = edgeList;
