function computeActinGFPEdgesDiff(dataDirectory, analysisDirectory)

if nargin < 1 || isempty(dataDirectory)
    dataDirectory = uigetdir('', 'Select a data directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

if nargin < 2 || isempty(analysisDirectory)
    analysisDirectory = uigetdir('', 'Select an analysis directory:');

    if ~ischar(dataDirectory)
        return;
    end
end

subFolders = {'cell2', 'cell4', 'cell6', 'cell7'};
dataPaths = cellfun(@(x) [dataDirectory filesep x], subFolders,...
    'UniformOutput', false);

disp('List of directories:');

for iMovie = 1:numel(dataPaths)
    disp([num2str(iMovie) ': ' dataPaths{iMovie}]);
end

analysisPaths = cellfun(@(x) [analysisDirectory filesep x], ...
    subFolders, 'UniformOutput', false);

for iMovie = 1:numel(analysisPaths)
    if ~exist(analysisPaths{iMovie}, 'dir');
        mkdir(analysisPaths{iMovie});
    end
end

disp('Process all directories...');

nMovies = numel(dataPaths);

pixelSize = 67/1000;

for iMovie = 1:nMovies
    movieName = ['Movie ' num2str(iMovie) '/' num2str(nMovies)];
    
    dataPath = dataPaths{iMovie};
    analysisPath = analysisPaths{iMovie};
    
    % STEP 1: Extend the union of GFP and Actin masks
    pathGFP = [dataPath filesep 'GFP' filesep 'edge' filesep 'cell_mask'];
    fileGFP = dir([pathGFP filesep '*.tif']);
    fileGFP = {fileGFP(:).name};
    
    pathActin = [dataPath filesep 'actin' filesep 'edge' filesep 'cell_mask'];
    fileActin = dir([pathActin filesep '*.tif']);
    fileActin = {fileActin(:).name};
    
    if numel(fileGFP) ~= numel(fileActin)
        error([movieName ': GFP and Actin channels don''t have the same number of masks.']);
    end
    
    nFrames = numel(fileGFP);
    
    % Load GFP masks
    BW = arrayfun(@(iFrame) imread([pathGFP filesep fileGFP{iFrame}]), ...
        1:nFrames,'UniformOutput', false);
    
    % Load Actin masks and compute the union with GFP mask
    BW = arrayfun(@(iFrame) BW{iFrame} | imread([pathActin filesep ...
        fileActin{iFrame}]), 1:nFrames,'UniformOutput', false);
    
    % Dilate by 10 pixels the masks
    r = 10;
    se = strel('disk',r);
    BW = arrayfun(@(iFrame) imdilate(BW{iFrame}, se), 1:nFrames,...
        'UniformOutput', false);
    
    % STEP 2: Compute the distance transform
    D = arrayfun(@(iFrame) bwdist(1 - BW{iFrame}), 1:nFrames,...
        'UniformOutput', false);
    clear BW;
    
    % STEP 3: Load image channels
    pathGFP = [dataPath filesep 'GFP' filesep 'crop'];
    fileGFP = dir([pathGFP filesep '*.tif']);
    fileGFP = {fileGFP(:).name};
    
    pathActin = [dataPath filesep 'actin' filesep 'crop'];
    fileActin = dir([pathActin filesep '*.tif']);
    fileActin = {fileActin(:).name};
    
    if numel(fileGFP) ~= numel(fileActin)
        error([movieName ': GFP and Actin channels don''t have the same number images.']);
    end
    
    if numel(fileGFP) ~= nFrames
        error([movieName ': number of images differs from number of masks.']);
    end
        
    % Load GFP images
    imageGFP = arrayfun(@(iFrame) double(imread([pathGFP filesep fileGFP{iFrame}])), ...
        1:nFrames,'UniformOutput', false);
    
    % Load Actin images
    imageActin = arrayfun(@(iFrame) double(imread([pathActin filesep fileActin{iFrame}])), ...
        1:nFrames,'UniformOutput', false);

    % STEP 4: Sample image into bands X times    
    dContours = 0:40;
    
    sampleGFP = zeros(numel(dContours)-1, nFrames);
    sampleActin = zeros(numel(dContours)-1, nFrames);
    
    for iFrame = 1:nFrames
        % Define the bands;
        
        bands = arrayfun(@(i) D{iFrame} > dContours(i), ...
            1:numel(dContours), 'UniformOutput', false);
        
        % Get the 2*r-wide band and normalize pixels in [0-1]
        ind = find(bands{1} - bands{end});
        
        minI = min(imageGFP{iFrame}(ind));
        maxI = max(imageGFP{iFrame}(ind));
        imageGFP{iFrame}(ind) = (imageGFP{iFrame}(ind) - minI) / (maxI - minI);
        
        minI = min(imageActin{iFrame}(ind));
        maxI = max(imageActin{iFrame}(ind));
        imageActin{iFrame}(ind) = (imageActin{iFrame}(ind) - minI) / (maxI - minI);

        % Get the 1-wide bands
        ind = arrayfun(@(i) find(bands{i} - bands{i+1}), ...
            1:numel(dContours)-1, 'UniformOutput', false);
        
        % Get the average normalized intensity over each band
        sampleGFP(:,iFrame) = cell2mat(cellfun(@(i) mean(imageGFP{iFrame}(i)),...
            ind, 'UniformOutput',false));
        sampleActin(:,iFrame) = cell2mat(cellfun(@(i) mean(imageActin{iFrame}(i)),...
            ind, 'UniformOutput',false));
    end
    
    hFig = figure('Visible', 'off');

    set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');    

    x = dContours(1:end-1)*pixelSize;
    y1 = [mean(sampleGFP,2) mean(sampleActin,2)];
    y2 = [gradient(mean(sampleGFP,2)) gradient(mean(sampleActin,2))];

    [AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
        
    set(AX(1), 'FontName', 'Helvetica', 'FontSize', 14);
    set(AX(2), 'FontName', 'Helvetica', 'FontSize', 14);
    
    set(get(AX(1),'Ylabel'),'String','Normalized Intensity');
    set(get(AX(1),'Ylabel'),'String','Normalized Intensity');
    set(H1,'LineStyle','-');
    set(H1(1), 'Color', [.4 .8 .2]);
    set(H1(2), 'Color', [.6 .2 .2]);
    set(H2,'LineStyle','--');
    set(H2(1), 'Color', [.4 .8 .2]);
    set(H2(2), 'Color', [.6 .2 .2]);
    title(['Intensity Profile - ' movieName]);
    legend({'GFP', 'XRhod-Actin'});
    xlabel([char(181) 'm']);
    fileName = [analysisPath filesep 'SupplFig1_' num2str(iMovie) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
end