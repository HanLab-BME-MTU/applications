function [inputList,imgScalingFactor,maxColor,doImage,doSpots] = make3DspotImageBuildInputList(idlist,image,pixelSize)
%creates input data for make3DspotImage3BuildInputList
%
%SYNOPSIS [inputList,imgScalingFactor,maxColor,doImage,doSpots] = make3DspotImageBuildInputList(idlist,image,pixelSize)
%
%INPUT    idlist    : (opt) idlist
%         image     : (opt) movie that hopefully corresponds to the idlist
%         pixelSize : (opt) pixelsize in micron
%
%         if any of the inputs are missing, the program will ask
%         WARNING: input checks are very rudimentary!
%
%OUTPUT   inputList        : inputList for make3DspotImage
%         imgScalingFactor : one of the fields of inputProperties for
%                            make3DspotImage
%         maxColor
%         doImage          : whether to plot images
%         doSpots          : whether to plot spots
%
%c: 01/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define defaults
ask4image = 1;
ask4idlist = 1;
ask4pixelSize = 1;

spotSize = 4; %could be adapted in the future to reflect psf or intensity

%--test input

if nargin > 0 & ~isempty(idlist)
    ask4idlist = 0;
    idlist = [];
end

if nargin > 1 & ~isempty(image)
    ask4image = 0;
    image = [];
end

if nargin > 2 & ~isempty(pixelSize)
    ask4pixelSize = 0;
    pixelSize = [];
end

%store old path    
oldPath = pwd;                
                
 %--end test input               
                
%load idlist              
if ask4idlist    
    %check if default biodata-dir exists and cd if exist
    cdBiodata(2);
    
    %1) load idlist to calculate inputList.spots
    
    %tell the user what's going on
    h = helpdlg('Please load an idlist or a trajectory. If you do not want to, press cancel','load idlist/trajectory');
    uiwait(h);
    
    %let the user choose. filterIdx will return 0 if cancel
    [fileName,pathName,filterIdx] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files';...
            'idlist*','idlists'},... %add option to load saved trajectory here
        'load idlist/trajectory');
    
    switch filterIdx
        case 0
            %no spots
            doSpots = 0;
            idlist = [];
            
        case 1
            
            cd(pathName);
            
            %dataFile. have the user choose an idlist
            dataStruct = load([pathName,fileName]);
            dataFieldNames = fieldnames(dataStruct);
            idnameListIdx = strmatch('idlist',dataFieldNames);
            idnameList = dataFieldNames(idnameListIdx);
            
            %have the user choose, if there is more than one idlist
            switch length(idnameList)
                case 0 %no idlist loaded. continue w/o loading
                    idname = '[]';
                case 1 %only one idlist loaded. Continue
                    idname = char(idnameList);
                otherwise %let the user choose
                    idSelect = chooseFileGUI(idnameList);
                    if isempty(idSelect)
                        idname = '[]';
                    else
                        idname = idnameList{idSelect};
                    end
            end
            
            if ~isempty(idname)
                %assign idlist
                eval(['idlist = dataStruct.',idname,';']);
                doSpots = 1;
            else
                %complain & no doSpots
                h = errordlg('No idlist found in project data - no trajectory will be displayed','Warning!');
                uiwait(h);
                doSpots = 0;
            end
            
            if ask4pixelSize
                %we're loading dataProperties - we can get the pixelsize
                %immediately
                try
                    pixelSize = [dataStruct.dataProperties.PIXELSIZE_XY,...
                            dataStruct.dataProperties.PIXELSIZE_XY,...
                            dataStruct.dataProperties.PIXELSIZE_Z];
                    ask4pixelSize = 0;
                catch
                    h = warndlg('there was no dataProperties in the data file!','WARNING')
                    uiwait(h);
                end
            end
            
        case 2
            
            cd(pathName);
            
            %load idlist
            dataStruct = load([pathName,fileName]);
            dataName = fieldnames(dataName);
            
            %assign idlist
            eval(['idlist = dataStruct.',char(dataName),';']);
            doSpots = 1;
            
    end %switch filterIdx
    
end %if ask4idlist

if doSpots %ask for filler for idlist
    %1) load idlist to calculate inputList.spots
    
    %tell the user what's going on
    h = helpdlg('Please load a second idlist or a trajectory to fill in deleted frames. If you do not want to, press cancel',...
        'load filler idlist/trajectory');
    uiwait(h);
    
    %let the user choose. filterIdx will return 0 if cancel
    [fileName,pathName,filterIdx] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files';...
            'idlist*','idlists'},... %add option to load saved trajectory here
        'load idlist/trajectory');
    
    switch filterIdx
        case 0
            %no filler
            fillerIdlist = [];
            
        case 1
            
            cd(pathName);
            
            %dataFile. have the user choose an idlist
            dataStruct = load([pathName,fileName]);
            dataFieldNames = fieldnames(dataStruct);
            idnameListIdx = strmatch('idlist',dataFieldNames);
            idnameList = dataFieldNames(idnameListIdx);
            
            %have the user choose, if there is more than one idlist
            switch length(idnameList)
                case 0 %no idlist loaded. continue w/o loading
                    idname = '[]';
                case 1 %only one idlist loaded. Continue
                    idname = char(idnameList);
                otherwise %let the user choose
                    idSelect = chooseFileGUI(idnameList);
                    if isempty(idSelect)
                        idname = '[]';
                    else
                        idname = idnameList{idSelect};
                    end
            end
            
            if ~isempty(idname)
                %assign idlist
                eval(['fillerIdlist = dataStruct.',idname,';']);
            else
                %complain & die
                error('No idlist found in project data - no trajectory will be displayed','Warning!');
            end
            
            
            
        case 2
            
            cd(pathName);
            
            %load idlist
            dataStruct = load([pathName,fileName]);
            dataName = fieldnames(dataName);
            
            %assign idlist
            eval(['fillerIdlist = dataStruct.',char(dataName),';']);
            
            
    end %switch filterIdx
end %ask for filler for idlist

if ask4image
    
    %2) load image to calculate inputList.image
    
    %tell the user what's going on
    h = helpdlg('Please load an image file. If you do not want to, press cancel','load idlist/trajectory');
    uiwait(h);
    
    %let the user choose. filterIdx will return 0 if cancel
    [fileName,pathName,filterIdx] = uigetfile({'*.r3c;*.fim;*.r3d','movie files'},... 
        'load image file');
    
    switch filterIdx
        case 0
            %no image
            doImage = 0;
            
        case 1
            
            cd(pathName);
            
            %load movie file
            try
                if isempty(findstr(fileName,'.r3d'))
                    [image] = readmat([pathName,fileName]); %works for *.fim and *.r3c and moviedat
                else
                    [image] = r3dread([pathName,fileName]);
                end
                doImage = 1;
            catch
                h = errordlg(['no file loaded - ',lasterr]);
                uiwait(h);
                doImage = 0
            end
            
    end %switch filterIdx
    
end %ask4image

if ask4pixelSize
    
    %have the user load a data file to load dataProperties
    
    %tell the user what's going on
    h = helpdlg('Please load dataProperties from a data file or a temporary dataProperties file. If you do not want to, press cancel','load scaling factor');
    uiwait(h);
    
    %let the user choose. filterIdx will return 0 if cancel
    [fileName,pathName,filterIdx] = uigetfile({'*-data-??-???-????-??-??-??.mat;tmpDataProperties.mat','data/dataProperties files'},...
        'load scaling factor');
    
    switch filterIdx
        case 0
            imgScalingFactor = [];
        case 1
            dpStruct = load([pathName,fileName]);
            try
                pixelSize = [dpStruct.dataProperties.PIXELSIZE_XY,...
                            dpStruct.dataProperties.PIXELSIZE_XY,...
                            dpStruct.dataProperties.PIXELSIZE_Z];
            catch
                error(['no pixelSize found in ',pathName, fileName])
            end
    end
    
end %ask4pixelSize

%check that we have at least either image or idlist
if ~(doImage | doSpots)
    error('need at least either image or idlist!')
end

%now we have all the data in place to start filling the list, which will be
%done in a loop for every timepoint.

%---prepare loop

%get number of timepoints
idlistLength = length(idlist); %returns 0 for [], size returns 1 for []
imageLength  = size(image,5); %WORKS ONLY FOR 5-D ARRAYS! (do right size after loading)
numTimePoints = max(idlistLength,imageLength);


%init vars
switch doImage + 2*doSpots
    case 1 %only image
        inputList = struct('image',[]);
    case 2 %only spots
        inputList = struct('spots',[]);
    case 3 %both
        inputList = struct('image',[],'spots',[]);
    otherwise
        %we have an errormessage above, but it can't hurt to repeat
        error('need at least either image or idlist!')
end

%start loop
for t = 1:numTimePoints
    if doSpots & t <= idlistLength
        %check whether there is data for this timePoint, then init matrix
        %and fill according to spots
        
        %get linklist - from filler if isempty
        linklist = idlist(t).linklist;
        if isempty(linklist) & ~isempty(fillerIdlist)
            linklist = fillerIdlist(t).linklist;
        end
        
        %fill in spots
        if ~isempty(linklist)
            %get # of spots
            numSpots = max(linklist(:,2));
            %init matrix
            spotMat = zeros(numSpots,5);
            %fill matrix
            for nsp = 1:numSpots
                %find where we store this spot
                spRowIdx = find(linklist(:,2)==nsp);
                %read and fill in coords
                coord = linklist(spRowIdx(1),9:11);
                spotMat(nsp,1:3) = coord./pixelSize;
                %read and fill in color
                spotMat(nsp,4) = mean(linklist(spRowIdx,4));
                %fill in size
                spotMat(nsp,5) = spotSize;
            end %or nsp = 1:numSpots
            %store spotMat
            inputList(t).spots = spotMat;
        end %~isempty(linklist)
    end %if doSpots
    
    if doImage & t <= imageLength
        %WARNING works for 5D only
        inputList(t).image = image(:,:,:,1,t);
    end
    
end %for t = 1:timePoints

%and now, create the imgScalingFactor
imgScalingFactor = pixelSize/min(pixelSize(:));

%maxColor
if doSpots
    maxColor = idlist(1).stats.maxColor;
else
    maxColor = [];
end

%switch back to old path
cd(oldPath);