function makeFigure4(paths, outputDirectory)

for iTM = 1:3
    % Load Movie Data
    fileName = [paths{iTM} filesep 'windowAnalysis' filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    %Verify that the labeling has been performed
    if ~checkMovieLabels(movieData)
        error('Must label movie before computing figure 3.');
    end

    % Get the names of the 2 FSM subfolders
    names = cellfun(@fliplr, strtok(cellfun(@fliplr,movieData.fsmDirectory, ...
        'UniformOutput', false), filesep), 'UniformOutput', false);

    % Make sure TMs is the first directory in the list
    iCh = [1, 2];
    if ~strcmpi(names{2}, 'actin')
        iCh = [2, 1];
        names = fliplr(names);
    end
    % Force channel 2's name to be 'Actin'
    names{2} = 'Actin';
    
    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
    % Read the list of label files
    labelPath = movieData.labels.directory;
    labelFiles = dir([labelPath filesep '*.tif']);

    % Read the list of TMs speckles (channel 1)
    s1Path = [movieData.fsmDirectory{iCh(1)} filesep 'tack' filesep 'locMax'];
    s1Files = dir([s1Path filesep '*.mat']);

    % Read the list of Actin speckles (channel 2)
    s2Path = [movieData.fsmDirectory{iCh(2)} filesep 'tack' filesep 'locMax'];
    s2Files = dir([s2Path filesep '*.mat']);
    
    % Read the list of Actin masks
    maskPath = movieData.masks.directory;
    maskFiles = dir([maskPath filesep '*.tif']);

    % Load activity map
    fileName = [movieData.protrusion.directory filesep ...
        movieData.protrusion.samples.fileName];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);
    v = sort(protrusionSamples.averageMagnitude(:));
    val = v(ceil(.01 * numel(v)));
    activityMap = protrusionSamples.averageNormalComponent(...
        protrusionSamples.averageNormalComponent > val & ...
        protrusionSamples.averageNormalComponent < -val);    
    activityMap = activityMap * pixelSize / timeInterval;
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 4 PANEL A                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    timeScale = 0:timeInterval:(nFrames-2)*timeInterval;

    hFig = figure('Visible', 'off');    
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
    imagesc(activityMap);
    set(gca,'XTick', timeScale);
    xlabel('Time (s)');
    if iTM == 1
        title('Edge Velocity (nm.s^{-1})');
        ylabel('Window no.');
    elseif iTM == 3
        colorbar;
    end
    fileName = [outputDirectory filesep 'fig4_A' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);    
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 4 PANEL B                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    distanceMap = cell(1, length(timeScale));
    
    for iFrame = 1:nFrames-1
        % Load label
        L = imread([labelPath filesep labelFiles(iFrame).name]);
        labels = 1:max(L(:));
        
        % Load TM speckles (channel 1)
        load([s1Path filesep s1Files(iFrame).name]);
        locMax1 = locMax;
        idxS1 = arrayfun(@(l) (locMax1 .* (L == l)) ~= 0, labels, ...
            'UniformOutput', false);
        
        % Load Actin speckles (channel 2)
        load([s2Path filesep s2Files(iFrame).name]);
        locMax2 = locMax;
        idxS2 = arrayfun(@(l) (locMax2 .* (L == l)) ~= 0, labels, ...
            'UniformOutput', false);
        
        % Compute distance to the edge
        BW = imread([maskPath filesep maskFiles(iFrame).name]);        
        distToEdge = double(bwdist(1 - BW)) * (pixelSize / 1000); % in microns

        % Compute distanceMap
        distanceMap{iFrame} = cellfun(@(l) mean(distToEdge(idxS1{l})) - ...
            mean(distToEdge(idxS2{l})), labels, 'UniformOutput', false);
    end
    
    nWindows = cellfun(@length, distanceMap, 'UniformOutput', false);
    nWindows = cat(1, nWindows{:});
    maxNWindows = max(nWindows);
    distanceMap = arrayfun(@(iFrame) padarray(distanceMap{iFrame}, ...
        [0 maxNWindows - nWindows(iFrame)], 0, 'post'), 1:nFrames-1, ...
        'UniformOutput', false);
    distanceMap = cat(2, distanceMap{:});
    
    hFig = figure('Visible', 'off');    
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
    imagesc(distanceMap);
    set(gca,'XTick', timeScale);
    xlabel('Time (s)');
    title([names(iCh(1)) ' Distance to Actin Front (' char(181) 'm)']);
    if iTM == 1
        ylabel('Window no.');
    end
    fileName = [outputDirectory filesep 'fig4_B' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);
end
