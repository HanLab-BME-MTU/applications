function edit_cropMovie_PB_Callback(hObject, eventdata, handles)
%callback to allow cropping of movie in editPropertiesGUI

%-----------display figure-----------------

%init coords and filteredMovie
coords = [];
filteredMovie = [];
safetyFactor = 1.5/min(handles.sigmaCorrection); %safety to make sure that all tags are within selection
safetyFactor = max(safetyFactor,1.1);

%if no filterprm, no filtering/detecting is possible
if ~isfield(handles,'FILTERPRM')
    errordlg('Set NA first!','Insufficient Information');
    return
end

%get movie length
movieLength = str2double(get(handles.edit_movieLength_txt,'String'));
if movieLength == 1
    representativeFrames = 1;
    nRepFrames = 1;
else
    nRepFrames = floor(movieLength*.2); %take 20% representative frames
    chooseStep = max(floor(movieLength/(nRepFrames-1)),1);
    representativeFrames = zeros(1,nRepFrames);
    representativeFrames = [1:chooseStep:movieLength];
    %make sure it's nRepFrames frames (it's nrp-1 if movieLength can be divided by nrp)
    if representativeFrames(end) == 0
        representativeFrames(end) = movieLength;
    end
end

%find whether there is a filtered movie. If not, load unfiltered movie and
%filter 15 (hopefully) representative frames
if ~isempty(handles.previewData.filteredMovieName)
    %there is one. display first frame of fim
    
    %read representative frames
    fimName = handles.previewData.filteredMovieName;
    filteredMovie=readmat(fimName,representativeFrames);
    %filteredImg is the first frame
    filteredImg = filteredMovie(:,:,:,1,1);
    
else %there is no filtered image (yet)
    
    %get name of raw movie
    movieName = handles.previewData.movieName;
    fimName = movieName;
    
    filter = questdlg(['There is no filtered movie yet. Do you want to crop without knowing where the spots are',...
            ' or do you want to filter ',num2str(nRepFrames),' representative frames (takes < 1 min, very recommended)?'],...
        'Do you want to filter?','Filter','Don''t filter','Filter');
    
    switch filter
        case 'Filter'
            %filter representative frames
            
            %read movie
            switch movieName(end-3:end)
                case '.r3d' %it's a non-cropped movie
                    %read movie and do background subtraction if necessary
                    [rawMovie,movieName,handles.oldBGState] = correctBackground(movieName,handles.currentBGState,handles.oldBGState);
                    handles.previewData.movieName = movieName;
                    guidata(hObject,handles);
                    %read specific frames
                    rawMovie = rawMovie(:,:,:,:,representativeFrames);
                    
                case '.r3c' %it's an already cropped movie
                    rawMovie = readmat(movieName, representativeFrames);
            end
            
            %filter movie
            filteredMovie = filterMovie(rawMovie,[handles.FILTERPRM]);
            filteredImg = filteredMovie(:,:,:,1,1);
            
        case ''
            %user closed window -> cancel
            return
            
        otherwise % don't filter: therefore read only 1 frame (don't do any background subtraction either)
            
            %read movie
            switch movieName(end-3:end)
                
                case '.r3d' %it's a non-cropped movie
                    rawMovie = r3dread(movieName,1,1);
                    
                case '.r3c' %it's an already cropped movie
                    rawMovie = readmat(movieName, 1);
            end
            
            %filter movie
            filteredImg = rawMovie(:,:,:,1,1); 
    end
    
end

%maximum projection
maxProj=max(filteredImg(:,:,:),[],3);
maxProj=maxProj/max(maxProj(:));

%show it
cropFigH = uiViewPanel;
set(cropFigH,'Name',['Crop Figure for ',fimName])
imshow(maxProj);
hold on;

%if there is a filtered movie (now), load slist or run spotfind on
%representative frames

if handles.status > 2
    %there is a slist. Load it and extract coords
    load(handles.projData,'slist');
    tmax = length(slist);
    
    
    for t=1:tmax
        if ~isempty (slist(t).sp)
            %take coordinates; no matter which frame they belong to!
            coords = [coords;cat(1,slist(t).sp.cord)];
        end
    end
    %we don't need the z-data (no z-cropping)
    coords = coords(:,1:2);
    
elseif ~isempty(filteredMovie)  
    %there is a filtered movie, but no slist: detect spots in representative
    %frames
    
    %build dataProperties
    dataProperties.FILTERPRM = handles.FILTERPRM;
    dataProperties.CH_MAXNUMINTERV =  1000;
    CH_MAXSLOPE = str2double(get(handles.edit_maxslope_txt,'String'));
    dataProperties.CH_MAXSLOPE = CH_MAXSLOPE;
    
    %run spotfind
    slist = spotfind(filteredMovie,dataProperties);
    for t=1:nRepFrames
        if ~isempty (slist(t).sp)
            %take coordinates; no matter which frame they belong to!
            coords = [coords;cat(1,slist(t).sp.cord)];
        end
    end
    %we don't need the z-data (no z-cropping)
    coords = coords(:,1:2);
    
    %increase safetyFactor
    safetyFactor = 1.4;
else
    %user does not want to see spots: do nothing
end


%if there are representative spots, plot them on top of the image

if ~isempty(coords)
    
    %plot coords
    plot(coords(:,1),coords(:,2),'+g');
    
    %get min/max coords to generate suggested bounding rectangle
    minC = min(coords,[],1);
    maxC = max(coords,[],1);
    
    deltaC = maxC - minC;
    
    %there should be at least a space of 1 sigma around each tag position
    sigmaXY = handles.FILTERPRM(1:2)*0.61/0.21;
    
    %set new min/delta coords
    minC = minC - safetyFactor*sigmaXY;
    deltaC = deltaC + 2*safetyFactor*sigmaXY;
    
    %plot rectangle
    rectangle('Position',[minC,deltaC],'EdgeColor','b');
    
end


crop = 'Cancel (select new)';
selectedRectangleH = [];

while strcmp(crop,'Cancel (select new)')
    
    if ~isempty(selectedRectangleH)
        delete(selectedRectangleH)
    end
    
    h = helpdlg('Please select your cropping rectangle. To cancel cropping, close the cropping window.');
    uiwait(h)
    
    %this part is from help rbbox: draw rectangle
    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    selectPosMin = min(point1,point2);        % lower left corner
    selectPosMax = max(point1,point2);        % upper right corner
    
    if ~ishandle(cropFigH) %user cancelled by closing cropWindow
        delete(gca); %calling gca above created axes on the EPGUI-window
        return %end evaluation here
    end
    
    if any(selectPosMin == selectPosMax)
        %user selecte only a point/line
        h = warndlg('Please select a rectangle','Insufficient input');
        uiwait(h);
        
    else %good input
        
        %if we know some coordinates, make sure that the rectangle is large enough   
        if ~isempty(coords)
            coordsIn = coords(find( (coords(:,1) > selectPosMin(1)) & (coords(:,2) > selectPosMin(2)) &...
                (coords(:,1) < selectPosMax(1)) & (coords(:,2) < selectPosMax(2))),:);
            minCI = min(coordsIn,[],1) - safetyFactor*sigmaXY;
            maxCI = max(coordsIn,[],1) + safetyFactor*sigmaXY;
            
            %allow larger-than-necessary rectangle
            selectPosMin = min(minCI,selectPosMin);
            selectPosMax = max(maxCI,selectPosMax);
        end
        
        %check that the selection is not out of bounds
        selectPosMin = max(floor(selectPosMin),[1,1]); %minimum: 1,1
        oldXY = [handles.header.numCols,handles.header.numRows];
        selectPosMax = min(ceil(selectPosMax),oldXY); %maximum: old max
        selectDelta = abs(selectPosMax - selectPosMin)+1; %add 1 because movie goes from min:max!
        
        
        selectedRectangleH = rectangle('Position',[selectPosMin,selectDelta],'EdgeColor','r');
        
        crop = questdlg('Do you want to crop?','','Yes','No','Cancel (select new)','Cancel (select new)');
    end %test if rectangle
end %while-loop. Program loops until Yes or No are selected



switch crop
    case 'No' %cancel cropping
        close(cropFigH);
        return
    case 'Yes' %crop the movie
        %load and crop
        try
            switch movieName(end-3:end)
                
                case '.r3d' %it's a non-cropped movie
                    movie = r3dread(movieName);
                    
                case '.r3c' %it's an already cropped/corrected movie
                    movie = readmat(movieName);
            end
            
            movieLength = size(movie,5);
            stackSize = size(movie,3);
            newMovie = zeros(selectDelta(1),selectDelta(2),stackSize,1,movieLength);
            newMovie = movie(selectPosMin(2):selectPosMax(2),...
                selectPosMin(1):selectPosMax(1),:,:,:);
        catch
            %in the future: do it chunkwise if out of memory
            rethrow(lasterr);
        end
        
        
        %-----------------save movie (r3c-file) in new directory----------------
        
        %set new projName
        newProjName = [handles.previewData.movieName(1:end-4),'_crop'];
        cropNum = [];
        
        cd ..
        if exist(newProjName,'dir')
            remove = questdlg('There is already a cropped movie for this project. Should it be overwritten?',...
                'WARNING','Overwrite','Don''t overwrite','Cancel cropping','Overwrite');
            switch remove
                case 'Overwrite'
                try
                    rmdir(newProjName,'s');
                catch
                    disp('if you run Windows 98 or Millenium: remove old crop-directory manually')
                    rethrow(lasterr);
                end
                
            case 'Don''t overwrite'
                cropNum = 1;
                while exist(newProjName,'dir')
                    cropNum  = cropNum +1;
                    newProjName = [handles.previewData.movieName(1:end-4),'_crop',num2str(cropNum)];
                end
            otherwise
                return
                %end evaluation here
            end
        end
        
        %create new directory
        mkdir(newProjName);
        cd(newProjName);
        
        %set new MovieName
        newMovieName = [newProjName,'.r3c'];
        writemat(newMovieName,newMovie);
        
        %create and save new projData / dataProperties etc
        newProjDataName = [newProjName,'-data-',nowString];
        
        %create current dataProperties/projProperties etc. file
        edit_readguidata(handles);
        
        %load handles, projectProperties etc
        amgH = handles.amgHandles.AMG;
        amgHandles = guidata(amgH);
        activeJobNum = handles.activeJobNum;
        
        projProperties = amgHandles.job(activeJobNum).projProperties;
        dataProperties = amgHandles.job(activeJobNum).dataProperties;
        
        %update projProperties. If cropNum is empty, there is a resulting empty
        %string
        projProperties.status = 0;
        %make sure we don't get into troubles with corr-movies
        filesepList = findstr(projProperties.dataPath,filesep);
        projProperties.dataPath = [projProperties.dataPath(1:filesepList(end)),newProjName,num2str(cropNum)];
        projProperties.datafileName = newProjDataName;
        
        %update dataProperties
        dataProperties.crop = [selectPosMin,selectPosMax;selectDelta,oldXY];
        dataProperties.name = newProjName;
        
        %update amg-list-entry
        dirList = get(amgHandles.dirListBox,'String');
        dirList{activeJobNum} = newProjName;
        %store new dirList
        set(amgHandles.dirListBox,'String',dirList);
        
        %set createNew to 0 (we want to keep our dataFile!)
        amgHandles.job(activeJobNum).createNew = 0;
        
        %store data in amgHandles
        amgHandles.job(activeJobNum).dataProperties = dataProperties;
        amgHandles.job(activeJobNum).projProperties = projProperties;
        amgHandles.job(activeJobNum).projName = newProjName;
        amgHandles.job(activeJobNum).projData = [newProjDataName,'.mat'];
        
        %save data to disk (including header!)
        save(newProjDataName,'dataProperties');
        save(newProjDataName,'projProperties','-append');
        r3dMovieHeader = handles.header;
        r3dMovieHeader.cropInfo = dataProperties.crop;
        save('r3dMovieHeader','r3dMovieHeader');
        
        %save amgHandles
        guidata(amgH,amgHandles);
        amgHandles = guidata(amgH);
        
        %close cropFigure
        close(cropFigH);
        
        %close and relaunch EPGUI
        close(handles.EPGUI);
        analyzeMoviesGUI('amg_editProp_PB_Callback',amgHandles.amg_editProp_PB, [], amgHandles);
        
    otherwise %cancelled by closing the window
        close(cropFigH);
        return
        
end