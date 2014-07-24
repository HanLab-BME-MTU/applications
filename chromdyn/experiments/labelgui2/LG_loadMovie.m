function success = LG_loadMovie(movie, movieDir, movieInfo, dataProperties)
%LG_loadMovie is the main function to load movies in labelgui2

success = 0;

% in case the movie has been partially loaded
% make new loadMovie-structure with fields
% - nFrames: number of frames that can be loaded into memory
% - loadedFrames: list of currently loaded frames
% - loadInfo with {movieName, movieType}
if isstruct(movieInfo) && ~isempty(movieInfo.frames2load)
    loadMovieStruct.loadedFrames = movieInfo.loadedFrames;
    loadMovieStruct.nFrames = length(loadMovieStruct.loadedFrames);
    loadMovieStruct.loadInfo = ...
        {fullfile(movieInfo.moviePath, movieInfo.movieName), ...
        movieInfo.movieType};
else
    loadMovieStruct = [];
end

% launch new movie window
success = LG_launchMovieWindow(...
    movie, movieDir, dataProperties, loadMovieStruct);

if success == 0
    h = warndlg('No movie loaded','Warning');
    uiwait(h)
    return
end

[naviHandles,movieWindowHandles] = LG_getNaviHandles;

% adjust sliders on navigator. UserData of the timepointSlider contains the
% previous timepoint. Initialize here to 0 because we don't want current
% and previous timepoint to match!
set(naviHandles.LG_navi_timepointSlider_sli,'Max',dataProperties.movieSize(4)+eps,...
    'Value',1,'Min',1,'UserData',0,'SliderStep',...
    [1/dataProperties.movieSize(4),10/dataProperties.movieSize(4)]);
set(naviHandles.LG_navi_zSliceSlider_sli,'Max',dataProperties.movieSize(3),...
    'Value',1,'Min',1,'SliderStep',...
    [1/dataProperties.movieSize(3),2/dataProperties.movieSize(3)]);
% adjust other stuff on navigator
set(naviHandles.LG_navi_movieName_txt,'String',dataProperties.name);

% return success