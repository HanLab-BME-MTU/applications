function navigatorHandle = LG_loadAllFromOutside(movie,movieDir,loadMovieStruct,dataProperties,idlist,idname)
%LG_loadAllFromOutside makes all the necessary calls to load data into labelgui2
% movie: full movie or part
% movieDir : directory of moive
% loadMovieStruct : empty if entire movie fits in memory. Otherwise,
% structure with fields
%   .nFrames (number of frames already loaded)
%   .loadedFrames (list of loaded frames)
%   .loadInfo (cell with {[moviePath, movieName], movieType}
% dataProperties : dataProperties
% idlist : idlist
% idname : name of idlist (idlist, idlist_L, idlisttrack, idlisttrack_L)
%
% LG_loadAllFromOutside returns the handles to labelgui for uiwait


% first, we have to check that a labelgui2 is, in fact, open
[naviHandles] = LG_getNaviHandles(1);

% if no labelgui2: load
if isempty(naviHandles)
    labelgui2;
    [naviHandles] = LG_getNaviHandles;
end

% load movie
success = LG_loadMovie(...
    movie, movieDir, loadMovieStruct, dataProperties);

if ~success
    error('error loading movie into labelgui2')
end

% load idlist
success = LG_loadIdlist(idlist, 0, idname, movieDir, []);

if ~success
    error('error loading idlist into labelgui2')
end

% set launchedFromOutside = 1, return handle if asked for
% careful! load naviHandles again, because otherwise, we save the wrong
% handles!
% also set dataHasChanged to 1 (different from original),
% so that we can always exit by saving.
[naviHandles, movieWindowHandles] = LG_getNaviHandles;
naviHandles.launchedFromOutside = 1;
movieWindowHandles.dataHasChanged = 1;
guidata(naviHandles.LG_navigator,naviHandles);
guidata(movieWindowHandles.LG_movieWindow,movieWindowHandles)

if nargout > 0
    navigatorHandle = naviHandles.LG_navigator;
end

    