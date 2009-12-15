function movieData = computeFigure1(movieData, batchMode)

%Indicate that Figure 1 computing was started
movieData.output.fig1.status = 0;

%Verify that the labeling has been performed
if ~checkMovieLabels(movieData)
    error('Must label movie before computing figure 1.');
end

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
    S1 = vertcat(cands(:).Lmax); %#ok<NODEF>
    clear cands;
    
    % Load speckles channel 2
    load([s2Path filesep s2Files(iFrame).name]);
    S2 = vertcat(cands(:).Lmax);
    clear cands;
    
    % Compute distance to the edge
    BW = imread([maskPath filesep maskFiles(iFrame).name]);
    distToEdge = bwdist(max(BW(:)) - BW) * pixelSize;
        
    % Compute linear indices of speckles
    idxS1 = sub2ind(size(distToEdge), S1(:, 1), S1(:, 2));
    idxS2 = sub2ind(size(distToEdge), S2(:, 1), S2(:, 2));
    
    maxLabel = max(L(:));
    
    D1{iFrame} = zeros(maxLabel, 50);
    D2{iFrame} = zeros(maxLabel, 50);
    
    for l = 1:maxLabel
        minD1 = sort(distToEdge(idxS1(L(idxS1) == l)));
        minD2 = sort(distToEdge(idxS2(L(idxS2) == l)));
        
        D1{iFrame}(l, :) = minD1(1:min(50, numel(minD1)));
        D2{iFrame}(l, :) = minD2(1:min(50, numel(minD2)));
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

movieData.output.fig1.dateTime = datestr(now);
movieData.output.fig1.status = 1;

updateMovieData(movieData);

end

