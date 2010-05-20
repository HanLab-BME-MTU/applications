function movieData = getMovieDensity(movieData, batchMode)

%Indicate that labeling was started
movieData.density.status = 0;

%Verify that the labels has been performed
if checkMovieLabels(movieData) && ~strcmp(movieData.labels.method,'window')
    error('Must create per-window label before create density maps.');
end

movieData.density.directory = fullfile(movieData.analysisDirectory, 'density');

if ~exist(movieData.density.directory, 'dir')
    mkdir(movieData.density.directory);
end

%Determine number of windows/bands
nBands = movieData.labels.nBands;
nSectors = movieData.labels.nSectors;
nFrames = movieData.labels.nFrames;

movieData.density.distanceValues = ...
    movieData.contours.parameters.distanceValues(1:nBands+1);

movieData.density.channels(1:2) = struct('averageMap', [], 'minMaxMap',  []);

% Read the list of label files
labelPath = movieData.labels.directory;
labelFiles = dir([labelPath filesep '*.tif']);

% Read the list of TMs speckles (channel 1)
s1Path = [movieData.fsmDirectory{1} filesep 'tack' filesep 'locMax'];
s1Files = dir([s1Path filesep '*.mat']);

% Read the list of Actin speckles (channel 2)
%s2Path = [movieData.fsmDirectory{2} filesep 'tack' filesep 'locMax'];
%s2Files = dir([s2Path filesep '*.mat']);

if ~batchMode
    h = waitbar(0,'Please wait, labeling window frames...');
end

averageMap1 = nan(nSectors,nFrames-1,nBands);
minMaxMap1 = nan(nSectors,nFrames-1,nBands);
    
for iFrame = 1:nFrames-1
    % Load label
    L = imread([labelPath filesep labelFiles(iFrame).name]);
    
    % Load channel 1 speckles
    load(fullfile(s1Path, s1Files(iFrame).name));
    ind = find(locMax);
    [I J] = ind2sub(size(locMax), ind);
    
    % Compute delaunay triangulation
    dt = DelaunayTri(I, J);
    % Find every index of triangle surrounding every point
    tIndices = vertexAttachments(dt);
    % Find every vertex index composing neighbors triangles (it includes itself)
    vIndices = cellfun(@(iTriangles) unique(dt.Triangulation(iTriangles,:)), ...
        tIndices, 'UniformOutput', false);
    % Compute distances
    dist = arrayfun(@(i) sqrt((I(vIndices{i}) - repmat(I(i), numel(vIndices{i}), 1)).^2 + ...
        (J(vIndices{i}) - repmat(J(i), numel(vIndices{i}), 1)).^2), (1:numel(vIndices))', ...
        'UniformOutput', false);
    % remove 0 distance
    dist = cellfun(@nonzeros, dist, 'UniformOutput', false);

    for iBand = 1:nBands
        % Discard first and last sector
        for iSector = 2:nSectors-1
            iWin = sub2ind([nBands, nSectors], iBand, iSector);
            
            isSpeckleInside = L(ind) == iWin;
            
            if any(isSpeckleInside)
                distInside = dist(isSpeckleInside);
                
                averageMap1(iSector,iFrame,iBand) = mean(cellfun(@mean, distInside));
                minMaxMap1(iSector,iFrame,iBand) = mean(cellfun(@(x) (max(x) - min(x)) / (max(x) + min(x)), distInside));
            end
        end
    end
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame/nFrames,h)
    end
end

movieData.density.channels(1).averageMap = averageMap1;
movieData.density.channels(1).minMaxMap = minMaxMap1;

if ~batchMode && ishandle(h)
    close(h);
end

if ~checkMovieProtrusionSamples(movieData)
    error('Must sample protrusion before labeling windows.');
end

movieData.density.dateTime = datestr(now);
movieData.density.status = 1;

updateMovieData(movieData);

end