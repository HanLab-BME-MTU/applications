function makeTropoFigure5(paths, outputDirectory)

for iTM = 1:3
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    %Verify that the labeling has been performed
    if ~checkMovieLabels(movieData)
        error('Must label movie before computing figure 3.');
    end
    
    % Get the names of the 2 FSM subfolders
    names = cellfun(@(x) x(max(regexp(x,filesep))+1:end),...
        movieData.fsmDirectory, 'UniformOutput', false);
    
    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
    % Read the list of label files
    labelPath = movieData.labels.directory;
    labelFiles = dir([labelPath filesep '*.tif']);

    % Read the list of TMs speckles (channel 1)
    s1Path = [movieData.fsmDirectory{1} filesep 'tack' filesep 'locMax'];
    s1Files = dir([s1Path filesep '*.mat']);

    % Read the list of Actin speckles (channel 2)
    s2Path = [movieData.fsmDirectory{2} filesep 'tack' filesep 'locMax'];
    s2Files = dir([s2Path filesep '*.mat']);

    % Read the list of distance transforms
    bwdistPath = movieData.bwdist.directory;
    bwdistFiles = dir([bwdistPath filesep '*.mat']);

    % Load activity map
    fileName = [movieData.protrusion.directory filesep ...
        movieData.protrusion.samples.fileName];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);
    % Convert activity maps to nm/s
    protrusionSamples.averageMagnitude = ...
        protrusionSamples.averageMagnitude * pixelSize / timeInterval;
    protrusionSamples.averageNormalComponent = ...
        protrusionSamples.averageNormalComponent * pixelSize / timeInterval;
    
    % Check whether the speed classification has been performed
    if ~isfield(protrusionSamples,'stateNames') || ...
            ~isfield(protrusionSamples,'states') || ...
            ~isfield(protrusionSamples,'statePersistence')
        error('Classification of edge velocity has not been performed.');
    end
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 5 PANEL A                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    states = protrusionSamples.states;
    statePersistence = protrusionSamples.statePersistence;
    
    prFrames = cell(1, max(statePersistence(states == 2)));
    reFrames = cell(1, max(statePersistence(states == 3)));
    
    for iFrame = 2:nFrames-2
        % Load label
        L = imread([labelPath filesep labelFiles(iFrame).name]);
        
        % Load TM speckles (channel 1)
        load([s1Path filesep s1Files(iFrame).name]);
        locMax1 = locMax;
        
        % Read the distance transform
        fileName = [bwdistPath filesep bwdistFiles(iFrame).name];
        load(fileName);
        distToEdge = distToEdge * (pixelSize / 1000); % in microns        

        % iterate over windows, discarding the 2 first and 2 last.
        for iWin = 3:max(L(:))-2
            ind = (locMax1 .* (L == iWin)) ~= 0;
            
            d = distToEdge(ind);
            
            p = statePersistence(iWin,iFrame);
            
            switch states(iWin,iFrame)
                case 2
                    prFrames{p} = vertcat(prFrames{p}, d);
                case 3
                    reFrames{p} = vertcat(reFrames{p}, d);
            end
        end
    end
    
    hFig = figure('Visible', 'off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');

    X = 1:numel(prFrames);
    Y = cellfun(@mean, prFrames);
    E = cellfun(@std, prFrames);
    
    errorbar(X(:),Y(:),E(:),'k');

    set(gca, 'XTick', 1:5:nFrames);
    
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x, '%.f'), (0:5:nFrames-1)*timeInterval, ...
        'UniformOutput', false));
    
    xlabel('Time (s)');
    ylabel([char(181) 'm']);

    fileName = [outputDirectory filesep 'Fig5_A' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 5 PANEL B                       %
    %                                                                 %
    %-----------------------------------------------------------------%

    hFig = figure('Visible', 'off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');

    X = 1:numel(reFrames);
    Y = cellfun(@mean, reFrames);
    E = cellfun(@std, reFrames);
    
    errorbar(X(:),Y(:),E(:),'k');
    
    set(gca, 'XTick', 1:5:nFrames);
    
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x, '%.f'), (0:5:nFrames-1)*timeInterval, ...
        'UniformOutput', false));
    
    xlabel('Time (s)');
    ylabel([char(181) 'm']);

    fileName = [outputDirectory filesep 'Fig5_B' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);    
end
