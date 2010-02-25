function makeFigure1(paths, batchMode)

if batchMode
    hFig = figure('Visible', 'off');
else
    hFig = figure('Visible', 'on');
end

for iCol = 1:3
    fileName = [paths{iCol} filesep 'windowAnalysis' filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    %Verify that the labeling has been performed
    if ~checkMovieLabels(movieData)
        error('Must label movie before computing figure 1.');
    end

    names = cellfun(@fliplr, strtok(cellfun(@fliplr,movieData.fsmDirectory, ...
        'UniformOutput', false), filesep), 'UniformOutput', false);

    % Make sure Actin is the first directory
    iCh = [1, 2];
    if ~strcmpi(names{1}, 'actin')
        iCh = [2, 1];
        names = fliplr(names);
    end
    
    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    
    if ~isfinite(movieData.labels.nBandsLimit)
        error('movieData.labels.nBandsLimit must be finite.');
    end
    
    sectorLength = movieData.labels.nBandsLimit * ...
        movieData.contours.parameters.distanceValues(2) * pixelSize;
    
    labelPath = movieData.labels.directory;
    labelFiles = dir([labelPath filesep '*.tif']);

    s1Path = [movieData.fsmDirectory{iCh(1)} filesep 'tack' filesep 'locMax'];
    s1Files = dir([s1Path filesep '*.mat']);

    s2Path = [movieData.fsmDirectory{iCh(2)} filesep 'tack' filesep 'locMax'];
    s2Files = dir([s2Path filesep '*.mat']);

    maskPath = movieData.masks.directory;
    maskFiles = dir([maskPath filesep '*.tif']);

    %nPanelA = 300;
    dataPanelA = zeros(2, nFrames-1);
    dataPanelB = cell(1, nFrames-1);
    dataPanelC = cell(1, nFrames-1);
    %dataPanelC_p = cell(1, nFrames-1);
    %dataPanelC_r = cell(1, nFrames-1);
    
    if ~batchMode
        h = waitbar(0, ['Making column ' num2str(iCol) ' in figure 1...']);
    end
    
    % load protrusion sample
    fileName = [movieData.protrusion.directory filesep ...
        movieData.protrusion.samples.fileName];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);
    % We do not consider the first and last 1% of the protrusion values
    % (too small/too large to be significant)
    %protValues = protrusionSamples.averageNormalComponent;
    %nValues = numel(protValues);
    %protValuesS = sort(protValues(:));
    %protMinCutOff = protValuesS(ceil(.01 * nValues));
    %protMaxCutOff = protValuesS(ceil((1 - .01) * nValues));
    %protValues(protValues > protMaxCutOff) = NaN;
    %protValues(protValues < protMinCutOff) = NaN;
    
    for iFrame = 1:nFrames-1
        % Load speckles channel 1
        load([s1Path filesep s1Files(iFrame).name]);
        idxS1 = find(locMax); 
        clear locMax;
        
        if isempty(idxS1)
            fprintf('%s channel doesn''t contain any speckle in frame %d/%d !!!\n',...
                names{1}, iFrame, nFrames);
            continue;
        end
        
        % Load speckles channel 2
        load([s2Path filesep s2Files(iFrame).name]);
        idxS2 = find(locMax); 
        clear locMax;
        
        if isempty(idxS2)
            fprintf('%s channel doesn''t contain any speckle in frame %d/%n !!!\n',...
                names{2}, iFrame, nFrames);
            continue;
        end
        
        % Compute distance to the edge
        BW = imread([maskPath filesep maskFiles(iFrame).name]);
        distToEdge = bwdist(max(BW(:)) - BW) * pixelSize;

        % Data for panel A       
        minD1 = sort(distToEdge(idxS1));
        minD2 = sort(distToEdge(idxS2));

        n1 = find(minD1 > 5000, 1);
        n2 = find(minD2 > 5000, 1);
        
        dataPanelA(1,iFrame) = mean(minD1(1:n1));
        dataPanelA(2,iFrame) = mean(minD2(1:n2));
        
        % Data for panel B        
        L = imread([labelPath filesep labelFiles(iFrame).name]);
        labels = (1:max(L(:)))';
        
        w = arrayfun(@(l) mean(distToEdge(L == l)), labels);
        w = sectorLength / (2 * w);
        
        idxS2 = arrayfun(@(l) idxS2(L(idxS2) == l), ...
            1:max(L(:)), 'UniformOutput', false);
        dataPanelB{iFrame} = arrayfun(@(l) w(l) * mean(distToEdge(idxS2{l})),...
            (1:max(L(:)))');
        
        % We put here the Actin distance as Panel C
        idxS1 = arrayfun(@(l) idxS1(L(idxS1) == l), ...
            1:max(L(:)), 'UniformOutput', false);
        
        dataPanelC{iFrame} = arrayfun(@(l) w(l) * mean(distToEdge(idxS1{l})),...
            (1:max(L(:)))');
        
        % Data for panel C        
%         idxL_p = find(protValues(:, iFrame) > 0);
%         idxL_r = find(protValues(:, iFrame) < 0);
%                
%         dataPanelC_p{iFrame} = cell2mat(arrayfun(@(l) distToEdge(idxS2{l}) - ...
%             mean(distToEdge(idxS1{l})), idxL_p, 'UniformOutput', false));
%         
%         dataPanelC_r{iFrame} = cell2mat(arrayfun(@(l) distToEdge(idxS2{l}) - ...
%             mean(distToEdge(idxS1{l})), idxL_r, 'UniformOutput', false));
%         
         if ~batchMode && ishandle(h)
             waitbar(iFrame / (nFrames-1), h);
         end
    end

    if ~batchMode && ishandle(h)
        close(h);
    end

    set(0, 'CurrentFigure', hFig);
    
    %
    % Panel A
    %
    
    subplot(3, 3, iCol);
    plot(gca, 1:nFrames-1, dataPanelA(1,:), 'r-', 'LineWidth', 2); hold on;
    plot(gca, 1:nFrames-1, dataPanelA(2,:), 'g-', 'LineWidth', 2); hold off;
    xlabel('Frame');
    h = kstest2(dataPanelA(1,:), dataPanelA(2,:));
    if h
        ylabel('Distance to edge (nm) (*)');
    else
        ylabel('Distance to edge (nm)');
    end
    legend(names);
    
    %
    % Panel B
    %

    nSectors = cellfun(@numel, dataPanelB);
    dataPanelB = cellfun(@(x) padarray(x, [max(nSectors) - numel(x), 0], 'post'), ...
        dataPanelB, 'UniformOutput', false);
    dataPanelB = horzcat(dataPanelB{:});
    subplot(3, 3, 3 + iCol);
    imagesc(dataPanelB, 'Parent', gca); hold on;
    plot(gca, 1:numel(nSectors), nSectors, '-w', 'LineWidth', 3);
    colormap('jet');
    colorbar;
    xlabel('Frame');
    ylabel('Sector #');
  
    %
    % Panel C
    %
  
    nSectors = cellfun(@numel, dataPanelC);
    dataPanelC = cellfun(@(x) padarray(x, [max(nSectors) - numel(x), 0], 'post'), ...
        dataPanelC, 'UniformOutput', false);
    dataPanelC = horzcat(dataPanelC{:});
    subplot(3, 3, 6 + iCol);
    imagesc(dataPanelC, 'Parent', gca); hold on;
    plot(gca, 1:numel(nSectors), nSectors, '-w', 'LineWidth', 3);
    colormap('jet');
    colorbar;
    xlabel('Frame');
    ylabel('Sector #');    
    
%     dataPanelC_p = vertcat(dataPanelC_p{:});
%     dataPanelC_r = vertcat(dataPanelC_r{:});
%     subplot(3, 3, 7:9);
%     [n1, xout1] = hist(dataPanelC_p, 50);
%     [n2, xout2] = hist(dataPanelC_r, 50);
%     c = rand(3, 1);
%     bar(xout1, n1, 'FaceColor', c); hold on;
%     bar(xout2, -n2, 'FaceColor', c * .5); hold off;
%     xlabel('nm');
end

