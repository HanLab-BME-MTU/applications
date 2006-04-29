function success = LG_navi_menuLoadMovie(hObject,eventdata,naviHandles)
%LG_navi_menuLoadMovie is the loader function for movie in labelgui2

% init success
success = 0;

% get naviHandles
naviHandles = LG_getNaviHandles;

% find biodata
biodataDir = cdBiodata(4);

% ask for movie. If cancelled, return
loadStruct.maxSize = 'check';
[movie, movieHeader, loadMovieStruct] = cdLoadMovie('ask',biodataDir,loadStruct);
if isscalar(movie) && movie == 0
    h = warndlg('No movie loaded','Warning');
    uiwait(h)
    return
end

%
% get movie directory
movieDir = loadMovieStruct.moviePath;

% find newest dataProperties
dataPropertiesList = searchFiles('dataproperties','',movieDir);

if isempty(dataPropertiesList)
    % create dataProperties from default
    dataProperties = defaultDataProperties(movieHeader);
    dataProperties.name = loadMovieStruct.movieName;
else
    % check for multiple dataProperties
    if size(dataPropertiesList,1) > 1
        % check whether there is tmpDataProperties!
        tmpIdx = strmatch('tmp',dataPropertiesList(:,1));
        if ~isempty(tmpIdx)
            dataPropertiesList(tmpIdx,:) = [];
        end
    end
    % check again for multiple dataProperties
    if size(dataPropertiesList,1) > 1
        % try to get movieName from loadMovieStruct. If it's a .fim, then
        % we remove the 'filtered_' and the extension. Otherwise, we have
        % the user decide
        if strcmp(loadMovieStruct.movieType,'filtered')
            dataName = loadMovieStruct.movieName;
            dataName = dataName(10:end-4);
            dataPropertiesIdx = strmatch(...
                sprintf('dataProperties_%s',dataName),dataPropertiesList(:,1));
        else
            % choose file
            dataPropertiesIdx = chooseFileGUI(dataPropertiesList(:,1));
        end

        % die if there's a problem
        if isempty(dataPropertiesIdx)
            h = errordlg('no valid dataProperties')
            uiwait(h)
            return
        end

        % load dataProperties
        load(fullfile(movieDir,dataPropertiesList{dataPropertiesIdx(1)}));

    else
        % load from file
        load(fullfile(dataPropertiesList{2},dataPropertiesList{1}))
    end
end


% load movie
success = LG_loadMovie(movie, movieDir, loadMovieStruct, dataProperties);

% if called via
if hObject == 0
    % return success - we've been called via loadAll
else
    % plot frame
    LG_gotoFrame(1);
end