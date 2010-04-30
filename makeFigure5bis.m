function makeFigure5bis(paths, outputDirectory)

for iTM = 1:numel(paths)
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    %Verify that the labeling has been performed
    if ~checkMovieLabels(movieData)
        error('Must label movie before computing figure 5bis.');
    end
    
    % Get the names of the 2 FSM subfolders
    names = cellfun(@fliplr, strtok(cellfun(@fliplr,movieData.fsmDirectory, ...
        'UniformOutput', false), filesep), 'UniformOutput', false);
    
    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
    % Read the list of label files
    labelPath = movieData.labels.directory;
    labelFiles = dir([labelPath filesep '*.tif']);

    % Read the list of TMx speckles (channel 1)
    s1Path = [movieData.fsmDirectory{1} filesep 'tack' filesep 'locMax'];
    s1Files = dir([s1Path filesep '*.mat']);

    % Read the list of TMy speckles (channel 2)
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
    
    % Load activity class
    fileName = [movieData.labels.directory filesep 'labelClasses.mat'];
    if ~exist(fileName, 'file')
        error(['Unale to locate ' fileName]);
    end
    load(fileName);
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 5bis PANEL A                    %
    %                                                                 %
    %-----------------------------------------------------------------%
    
    % Build activtyStep
    activityStep = zeros(size(labelClasses));
    
    % Forward sweep: initialize 1st frame so that protrution / retraction
    % events which start at first frame are discarded.
    activityStep(:,1) = 0;
    
    for iFrame = 2:nFrames-1
        % Protrusion at frames iFrame and iFrame-1
        ind = labelClasses(:,iFrame) == 1 & ...
            labelClasses(:,iFrame-1) == 1 & ...
            activityStep(:,iFrame-1);        
        activityStep(ind,iFrame) = activityStep(ind,iFrame-1) + 1;
        
        % Retration at frames iFrame and iFrame-1
        ind = labelClasses(:,iFrame) == 2 & ...
            labelClasses(:,iFrame-1) == 2 & ...
            activityStep(:,iFrame-1);
        activityStep(ind,iFrame) = activityStep(ind,iFrame-1) - 1;
        
        % Protrusion at frame iFrame but not at iFrame-1
        ind = labelClasses(:,iFrame) == 1 & labelClasses(:,iFrame-1) ~= 1;
        activityStep(ind,iFrame) = 1;
        
        % Retraction at frame iFrame but not at iFrame-1
        ind = labelClasses(:,iFrame) == 2 & labelClasses(:,iFrame-1) ~= 2;
        activityStep(ind,iFrame) = -1;
    end
    
    % Backward sweep: discard protrusion / retraction events which continue
    % at last frame.
    activityStep(:,end) = 0;
    
    for iFrame = nFrames-2:-1:1
        % Protrusion or Retration at frame iFrame and iFrame+1
        ind = labelClasses(:,iFrame) == labelClasses(:,iFrame+1) & ...
            labelClasses(:,iFrame) & ~activityStep(:,iFrame+1);
        activityStep(ind,iFrame) = 0;
    end
    
    prFrames1 = cell(1, max(activityStep(:)));
    reFrames1 = cell(1, -min(activityStep(:)));

    prFrames2 = cell(1, max(activityStep(:)));
    reFrames2 = cell(1, -min(activityStep(:)));

    for iFrame = 2:nFrames-2
        % Load label
        L = imread([labelPath filesep labelFiles(iFrame).name]);
        
        % Load TMx speckles (channel 1)
        load([s1Path filesep s1Files(iFrame).name]);
        locMax1 = locMax;

        % Load TMy speckles (channel 2)
        load([s2Path filesep s2Files(iFrame).name]);
        locMax2 = locMax;

        % Read the distance transform
        fileName = [bwdistPath filesep bwdistFiles(iFrame).name];
        load(fileName);
        distToEdge = distToEdge * (pixelSize / 1000); % in microns        

        % iterate over windows, discarding the 2 first and 2 last.
        for iWin = 3:max(L(:))-2
            ind1 = (locMax1 .* (L == iWin)) ~= 0;
            ind2 = (locMax2 .* (L == iWin)) ~= 0;
            
            d1 = distToEdge(ind1);
            d2 = distToEdge(ind2);
            
            step = activityStep(iWin,iFrame);
            
            if step > 0
                prFrames1{step} = vertcat(prFrames1{step}, d1);
                prFrames2{step} = vertcat(prFrames2{step}, d2);
            elseif step < 0
                reFrames1{-step} = vertcat(reFrames1{-step}, d1);
                reFrames2{-step} = vertcat(reFrames2{-step}, d2);
            end
        end
    end
    
    hFig = figure('Visible', 'off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');

    X1 = 1:numel(prFrames1);
    Y1 = cellfun(@mean, prFrames1);
    E1 = cellfun(@std, prFrames1);

    X2 = 1:numel(prFrames2);
    Y2 = cellfun(@mean, prFrames2);
    E2 = cellfun(@std, prFrames2);
    
    errorbar(X1(:),Y1(:),E1(:),'r'); hold on;
    errorbar(X2(:),Y2(:),E2(:),'b'); hold off;
    
    legend(names);
    
    xlabel('# of frame');
    ylabel([char(181) 'm']);

    fileName = [outputDirectory filesep 'Fig5bis_A' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 5bis PANEL B                    %
    %                                                                 %
    %-----------------------------------------------------------------%

    hFig = figure('Visible', 'off');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
    set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');

    X1 = 1:numel(reFrames1);
    Y1 = cellfun(@mean, reFrames1);
    E1 = cellfun(@std, reFrames1);

    X2 = 1:numel(reFrames2);
    Y2 = cellfun(@mean, reFrames2);
    E2 = cellfun(@std, reFrames2);
    
    errorbar(X1(:),Y1(:),E1(:),'r'); hold on;
    errorbar(X2(:),Y2(:),E2(:),'b'); hold off;
    
    legend(names);
    
    xlabel('# of frame');
    ylabel([char(181) 'm']);

    fileName = [outputDirectory filesep 'Fig5bis_B' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);    
end
