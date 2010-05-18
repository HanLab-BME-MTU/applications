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
    names = cellfun(@fliplr, strtok(cellfun(@fliplr,movieData.fsmDirectory, ...
        'UniformOutput', false), filesep), 'UniformOutput', false);
    % Force channel 2's name to be 'Actin'
    names{2} = 'Actin';
    
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
    if ~isfield(protrusionSamples,'classNames') || ~isfield(protrusionSamples,'classes')
        error('Classification of edge velocity has not been performed.');
    end
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 5 PANEL A                       %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    % Build activtyStep
    classes = protrusionSamples.classes;
    activityStep = zeros(size(classes));
    
    % Forward sweep: initialize 1st frame so that protrution / retraction
    % events which start at first frame are discarded.
    activityStep(:,1) = 0;
    
    for iFrame = 2:nFrames-1
        % Protrusion at frames iFrame and iFrame-1
        ind = classes(:,iFrame) == 1 & ...
            classes(:,iFrame-1) == 1 & ...
            activityStep(:,iFrame-1);        
        activityStep(ind,iFrame) = activityStep(ind,iFrame-1) + 1;
        
        % Retration at frames iFrame and iFrame-1
        ind = classes(:,iFrame) == 2 & ...
            classes(:,iFrame-1) == 2 & ...
            activityStep(:,iFrame-1);
        activityStep(ind,iFrame) = activityStep(ind,iFrame-1) - 1;
        
        % Protrusion at frame iFrame but not at iFrame-1
        ind = classes(:,iFrame) == 1 & classes(:,iFrame-1) ~= 1;
        activityStep(ind,iFrame) = 1;
        
        % Retraction at frame iFrame but not at iFrame-1
        ind = classes(:,iFrame) == 2 & classes(:,iFrame-1) ~= 2;
        activityStep(ind,iFrame) = -1;
    end
    
    % Backward sweep: discard protrusion / retraction events which continue
    % at last frame.
    activityStep(:,end) = 0;
    
    for iFrame = nFrames-2:-1:1
        % Protrusion or Retration at frame iFrame and iFrame+1
        ind = classes(:,iFrame) == classes(:,iFrame+1) & ...
            classes(:,iFrame) & ~activityStep(:,iFrame+1);
        activityStep(ind,iFrame) = 0;
    end
    
    prFrames = cell(1, max(activityStep(:)));
    reFrames = cell(1, -min(activityStep(:)));
    
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
            
            step = activityStep(iWin,iFrame);
            
            if step > 0
                prFrames{step} = vertcat(prFrames{step}, d);
            elseif step < 0
                reFrames{-step} = vertcat(reFrames{-step}, d);
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
