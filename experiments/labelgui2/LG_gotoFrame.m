function LG_gotoFrame(currentFrame, currentSlice)
%LG_gotoFrame collects the necessary data for plotting and sets the correct time everywhere


% get handles
[naviHandles, movieWindowHandles] = LG_getNaviHandles;

% currentSlice is optional
if nargin < 2 || isempty(currentSlice)
    currentSlice = [];
else % if there is currentSlice, lookup frame
    if isempty(currentFrame)
        currentFrame = LG_getCurrentTime;
    end 
end

% if the movieWindow has been closed: Don't continue
if isempty(movieWindowHandles)
    return
end

% 1) set the right frame/slice everywhere
% 2) get all the necessary data
% 3) call LG_plot

% set timeSlider, timeText
set(naviHandles.LG_navi_timepointNumber_edt, ...
    'String', num2str(currentFrame));
set(naviHandles.LG_navi_timepointSlider_sli, 'Value', currentFrame);

% set slice if applicable
if ~isempty(currentSlice)
    set(naviHandles.LG_navi_zSliceNumber_edt, ...
        'String', num2str(currentSlice));
    set(naviHandles.LG_navi_zSliceSlider_sli, 'Value', currentSlice);
end

% set flag text if applicable
if ~isempty(movieWindowHandles.flagData)
    % check if the currentFrame is contained in the currentFlag-List
    flagIdx = ...
        find(currentFrame == movieWindowHandles.flagData.flaggedFrameList);
    if isempty(flagIdx)
        set(naviHandles.LG_navi_flagNumber_txt,'String','');
    else
        nFlags = length(movieWindowHandles.flagData.flaggedFrameList);
        set(naviHandles.LG_navi_flagNumber_txt,'String',...
            sprintf('%i/%i',flagIdx,nFlags));
    end
end


% get the data: Read the full frame, idlist(t)
% if there is a nonempty loadMovieStruct, load another movie-chunk if
% necessary
if ~isempty(movieWindowHandles.loadMovieStruct)
    loadMovieStruct = movieWindowHandles.loadMovieStruct;
    % check whether the currentFrame has already been loaded
    if any(currentFrame == loadMovieStruct.loadedFrames)
        % all is good
    else
        % we need to load more. Calculate frames2Load. If the number of
        % frames to load is even, we round up so that there is one less
        % frame before the currentFrame than after
        startFrame = ceil(currentFrame - (loadMovieStruct.nFrames-1)/2);
        endFrame = ceil(currentFrame + (loadMovieStruct.nFrames-1)/2);
        % correct for frame out of range
        minFrame = 1;
        maxFrame = movieWindowHandles.dataProperties.movieSize(4);
        % if endFrame = maxFrame + 1, we have to subtract 1 from startFrame.
        % There will NEVER be more nFrames than dataProperties.movieSize(4)
        if endFrame > maxFrame
            startFrame = startFrame - endFrame + maxFrame;
            endFrame = maxFrame;
        end
        if startFrame < minFrame
            endFrame = endFrame + startFrame - minFrame;
            startFrame = minFrame;
        end

        % warn the user of loading
        mHandle=myMessageBox([],...
            sprintf(...
            'Please wait \nwhile labelgui2 loads\n frames %i:%i',...
            startFrame,endFrame),'Busy');
        
        % load the necessary frames
        movieWindowHandles.movie = ...
            cdLoadMovie(loadMovieStruct.loadInfo,[],[startFrame:endFrame]);
        
        % close messagebox
        if ishandle(mHandle)
            delete(mHandle)
        end
        % and store the loadedFrames
        movieWindowHandles.loadMovieStruct.loadedFrames = ...
            [startFrame:endFrame];
    end
end
% correct currentFrame by the number of the first loaded frame before
% indexing into movie!
if ~isempty(movieWindowHandles.loadMovieStruct)
    movieFrame = movieWindowHandles.movie(:,:,:,:,currentFrame - ...
    movieWindowHandles.loadMovieStruct.loadedFrames(1)+1);
else
    % of course, if we can load the entire movie, we don't need loadStruct
    movieFrame = movieWindowHandles.movie(:,:,:,:,currentFrame);
end

% read current idlist
idlist = movieWindowHandles.idlist;
if ~isempty(idlist)
    idlist = idlist(currentFrame);
end

% get plot info: List of axes handles, pixelSizes, aspectRatios
axesH = [movieWindowHandles.xyAxesH; ...
    movieWindowHandles.yzAxesH; ...
    movieWindowHandles.xzAxesH];
imageH = [movieWindowHandles.xyImageH; ...
    movieWindowHandles.yzImageH; ...
    movieWindowHandles.xzImageH];
frameSizeMu = movieWindowHandles.frameSizeMu;
pixelSize = [movieWindowHandles.dataProperties.PIXELSIZE_XY,...
    movieWindowHandles.dataProperties.PIXELSIZE_Z];

% set plotOptions
plotOptions.currentSlice = currentSlice;
if ~isempty(idlist)
    plotOptions.maxTags = movieWindowHandles.idlistData.maxTags;
    plotOptions.labelColor = movieWindowHandles.idlistData.labelcolor;
    plotOptions.colorMap = movieWindowHandles.colorMap;
    plotOptions.tagPopupMenuH = movieWindowHandles.tagPopupMenuH;

end


% plot data. Make movieWindowHandle visible, because LG_plot is not a
% callback
set(movieWindowHandles.LG_movieWindow,'HandleVisibility','on');
LG_plot(movieFrame, idlist, axesH, imageH, frameSizeMu, pixelSize, plotOptions);


% make windows visible, restrict handle visibility
set(movieWindowHandles.LG_movieWindow,'Visible','on','HandleVisibility','Callback');

% update movieWindow
LG_showMovieData_callback(1);


% store changes to movieWindowHandles
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles);

% make navigator active
figure(naviHandles.LG_navigator);