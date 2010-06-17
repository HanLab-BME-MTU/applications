function movieData = getMovieDensity(movieData, batchMode)

% Remove old data process
movieData.density = [];

% Indicate that labeling was started
movieData.density.status = 0;

% Verify that the distance transform has been computed
if checkMovieBWDist(movieData)
    error('Must compute distance transform first!');
end

movieData.density.directory = fullfile(movieData.analysisDirectory, 'density');

if ~exist(movieData.density.directory, 'dir')
    mkdir(movieData.density.directory);
end

% Get the name of the 2 subfolers
nChannels = numel(movieData.fsmDirectory);

movieData.density.channelDirectory = ...
    cellfun(@(x) x(max(regexp(x,filesep))+1:end), movieData.fsmDirectory, ...
    'UniformOutput', false);

for iChannel=1:nChannels
    if ~exist(fullfile(movieData.density.directory,...
            movieData.density.channelDirectory{iChannel}), 'dir')
        mkdir(movieData.density.directory,...
            movieData.density.channelDirectory{iChannel});
    end
end

% Determine number of frames
nFrames = movieData.labels.nFrames;

% Read the list of label files
bwDistPath = movieData.bwdist.directory;
bwDistFiles = dir([bwDistPath filesep '*.tif']);

% Read the list of speckle files
specklePaths = cellfun(@(x) fullfile(x, 'tack', 'locMax'), movieData.fsmDirectory, 'UniformOutput' ,false);
speckleFiles = cellfun(@(x) dir([x, filesep, '*.mat']), specklePaths, 'UniformOutput' ,false);

%Make the string for formatting
fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

if ~batchMode
    h = waitbar(0,'Please wait, labeling window frames...');
end
 
for iFrame = 1:nFrames
    % Load bwDist
    load([bwDistPath filesep bwDistFiles(iFrame).name]);
    
    % Add virtual points on cell edges
    outline = contourc(double(distToEdge), [0, 0]);
    vX = [];
    vY = [];
    iChunk = 1;
    while iChunk < size(outline,2)
        n = outline(2,iChunk);
        vX = [vX, outline(1,iChunk+1:10:iChunk+n)];
        vY = [vY, outline(2,iChunk+1:10:iChunk+n)];
        iChunk = iChunk + n + 1;
    end    

    for iChannel = 1:nChannels
    
        % Load speckles
        load(fullfile(specklePaths{iChannel}, speckleFiles{iChannel}(iFrame).name));
        ind = find(locMax ~= 0);
    
        nPoints = numel(ind);
        [Y X] = ind2sub(size(locMax), ind);
    
        % Compute delaunay triangulation
        dt = DelaunayTri([X; vX'],[Y; vY']);
    
        % Find every index of triangle surrounding every point
        triangles = vertexAttachments(dt);
    
        % Restrict triangle list attached to real points
        triangles = triangles(1:nPoints);
    
        % Find every vertex index composing neighbors triangles (it includes itself)
        vertices = cellfun(@(iTriangles) unique(dt.Triangulation(iTriangles,:)), ...
            triangles, 'UniformOutput', false);
    
        % Remove point which index exceeds nPoints
        vIndices = cellfun(@(iVertices) iVertices(iVertices <= nPoints), vertices, 'UniformOutput', false);
        
        % Compute pairwise distances
        dist = arrayfun(@(i) sqrt((X(vIndices{i}) - repmat(X(i), numel(vIndices{i}), 1)).^2 + ...
            (Y(vIndices{i}) - repmat(Y(i), numel(vIndices{i}), 1)).^2), (1:numel(vIndices))', ...
            'UniformOutput', false);
    
        % remove 0 distance
        dist = cellfun(@nonzeros, dist, 'UniformOutput', false);
    
        % Remove lonely speckles which are linked only with virtual points.
        validDist = cellfun(@(x) ~isempty(x), dist);
        
        % create average density score
        averageDensityScore = cellfun(@mean, dist(validDist));
   
        averageDensityMap = zeros(size(L));
        averageDensityMap(ind(validDist)) = averageDensityScore; %#ok<NASGU>
    
        save(fullfile(movieData.density.directory, ...
            movieData.density.channelDirectory{iChannel}, ...
            ['averageDensityMap_', num2str(iFrame,fString), '.mat']), ...
            'averageDensityMap');
    
        % create min/max density score
        minMaxDensityScore = cellfun(@(x) (max(x) - min(x)) / ...
            (max(x) + min(x)), dist(validDist));
    
        minMaxDensityMap = zeros(size(L));
        minMaxDensityMap(ind(validDist)) = minMaxDensityScore; %#ok<NASGU>
    
        save(fullfile(movieData.density.directory, ...
            movieData.density.channelDirectory{iChannel}, ...
            ['minMaxDensityMap_', num2str(iFrame,fString), '.mat']), ...
            'minMaxDensityMap');
    end
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame/nFrames,h)
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

movieData.density.dateTime = datestr(now);
movieData.density.status = 1;

updateMovieData(movieData);

end