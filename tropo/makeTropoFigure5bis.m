function makeTropoFigure5bis(paths, outputDirectory)

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
    
    % Check whether the speed classification has been performed
    if ~isfield(protrusionSamples,'stateNames') || ~isfield(protrusionSamples,'states')
        error('Classification of edge velocity has not been performed.');
    end
    
    %-----------------------------------------------------------------%
    %                                                                 %
    %                          FIGURE 5bis PANEL A                    %
    %                                                                 %
    %-----------------------------------------------------------------%
        
    states = protrusionSamples.states;
    statePersistence = protrusionSamples.statePersistence;

    prFrames1 = cell(1, max(statePersistence(states == 2)));
    reFrames1 = cell(1, max(statePersistence(states == 3)));

    prFrames2 = cell(1, max(statePersistence(states == 2)));
    reFrames2 = cell(1, max(statePersistence(states == 3)));

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
            
            if any(ind1(:))
                d1 = sort(distToEdge(ind1));
                d1 = d1(1);
            else
                d1 = [];
            end
            
            if any(ind2(:))
                d2 = sort(distToEdge(ind2));
                d2 = d2(1);
            else
                d2 = [];
            end
            
            p = statePersistence(iWin,iFrame);
            
            switch states(iWin,iFrame)
                case 2
                    prFrames1{p} = vertcat(prFrames1{p}, d1);
                    prFrames2{p} = vertcat(prFrames2{p}, d2);
                case 3
                    reFrames1{p} = vertcat(reFrames1{p}, d1);
                    reFrames2{p} = vertcat(reFrames2{p}, d2);
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
    
    set(gca, 'XTick', 1:5:nFrames);
    
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x, '%.f'), (0:5:nFrames-1)*timeInterval, ...
        'UniformOutput', false));
    
    xlabel('Time (s)');
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
    
    set(gca, 'XTick', 1:5:nFrames);
    
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x, '%.f'), (0:5:nFrames-1)*timeInterval, ...
        'UniformOutput', false));
    
    xlabel('Time (s)');
    ylabel([char(181) 'm']);

    fileName = [outputDirectory filesep 'Fig5bis_B' num2str(iTM) '.eps'];
    print(hFig, '-depsc', fileName);
    fixEpsFile(fileName);
    close(hFig);    
end
