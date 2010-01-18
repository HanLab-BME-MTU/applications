function movieData = computeFigure1(movieData, batchMode)

%Indicate that Figure 1 computing was started
movieData.output.fig1.status = 0;

%Verify that the labeling has been performed
if ~checkMovieLabels(movieData)
    error('Must label movie before computing figure 1.');
end

movieData.output.directory = [movieData.analysisDirectory filesep 'output'];

if ~exist(movieData.output.directory, 'dir')
    mkdir(movieData.output.directory);
end

movieData.output.fig1.filename = 'fig1.mat';

names = cellfun(@fliplr, strtok(cellfun(@fliplr,movieData.fsmDirectory, ...
    'UniformOutput', false), filesep), 'UniformOutput', false);

nFrames = movieData.labels.nFrames;
pixelSize = movieData.pixelSize_nm;

labelPath = movieData.labels.directory;
labelFiles = dir([labelPath filesep '*.tif']);

s1Path = [movieData.fsmDirectory{1} filesep 'tack' filesep 'cands'];
s1Files = dir([s1Path filesep '*.mat']);

s2Path = [movieData.fsmDirectory{2} filesep 'tack' filesep 'cands'];
s2Files = dir([s2Path filesep '*.mat']);

maskPath = movieData.masks.directory;
maskFiles = dir([maskPath filesep '*.tif']);

% Define the number of closest speckles away from the cell edge.
nSpeckles = 20;

D1 = cell(1, nFrames);
D2 = cell(1, nFrames);

%Go through each frame and save the windows to a file
if ~batchMode
    h = waitbar(0,'Please wait, compute figure 1 data...');
end

for iFrame = 1:nFrames
    % Load labels
    L = imread([labelPath filesep labelFiles(iFrame).name]);

    % Load speckles channel 1
    load([s1Path filesep s1Files(iFrame).name]);
    status = vertcat(cands(:).status); %#ok<NODEF>
    S1 = vertcat(cands(status == 1).Lmax);
    clear cands;
    
    if isempty(S1)
        fprintf('%s channel doesn''t contain any speckle in frame %d/%d !!!\n',...
            names{1}, iFrame, nFrames);
        continue;
    end
    
    % Load speckles channel 2
    load([s2Path filesep s2Files(iFrame).name]);
    status = vertcat(cands(:).status);
    S2 = vertcat(cands(status == 1).Lmax);
    clear cands;
    
    if isempty(S2)
        fprintf('%s channel doesn''t contain any speckle in frame %d/%n !!!\n',...
            names{2}, iFrame, nFrames);
        continue;
    end
    
    % Compute distance to the edge
    BW = imread([maskPath filesep maskFiles(iFrame).name]);
    distToEdge = bwdist(max(BW(:)) - BW) * pixelSize;
        
    % Compute linear indices of speckles
    idxS1 = sub2ind(size(distToEdge), S1(:, 1), S1(:, 2));
    idxS2 = sub2ind(size(distToEdge), S2(:, 1), S2(:, 2));
    
    maxLabel = max(L(:));
    
    D1{iFrame} = zeros(maxLabel, nSpeckles);
    D2{iFrame} = zeros(maxLabel, nSpeckles);
    
    for l = 1:maxLabel
        minD1 = sort(distToEdge(idxS1(L(idxS1) == l)));
        minD2 = sort(distToEdge(idxS2(L(idxS2) == l)));

        n1 = min(nSpeckles, numel(minD1));
        n2 = min(nSpeckles, numel(minD2));
        
        D1{iFrame}(l, 1:n1) = minD1(1:n1);
        D2{iFrame}(l, 1:n2) = minD2(1:n2);
    end
    
    if ~batchMode && ishandle(h)
        waitbar(iFrame / nFrames, h);
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

% In case some frames are empty.
D1 = D1(cellfun(@(x) ~isempty(x), D1)); %#ok<NASGU>
D2 = D2(cellfun(@(x) ~isempty(x), D2)); %#ok<NASGU>

save([movieData.output.directory filesep movieData.output.fig1.filename],...
    'names', 'D1', 'D2');

movieData.output.fig1.dateTime = datestr(now);
movieData.output.fig1.status = 1;

updateMovieData(movieData);

end

