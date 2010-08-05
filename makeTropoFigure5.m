function makeTropoFigure5(paths, outputDirectory)

nBands = 1;
dLims = [0 2 10] * 1000;
minLifetime = 3;
maxDist = 2000;

tmNames = {'TM2', 'TM4', 'TM5NM1'};

for iTM = 1:numel(paths)
    % Load Movie Data
    fileName = [paths{iTM} filesep 'movieData.mat'];
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    nFrames = movieData.labels.nFrames;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    
    % Read the list of distance transforms
    bwdistPath = movieData.bwdist.directory;
    bwdistFiles = dir([bwdistPath filesep '*.mat']);

    % Read distance transforms
    distToEdge = zeros(movieData.imSize(1),movieData.imSize(2),nFrames);
    
    for iFrame = 1:nFrames
        fileName = fullfile(bwdistPath, bwdistFiles(iFrame).name);
        tmp = load(fileName);
        distToEdge(:,:,iFrame) = tmp.distToEdge * pixelSize;
    end
    
    % Read the TM MPM
    load(fullfile(movieData.fsmDirectory{1}, 'tack', 'mpm.mat'));
    
    % Compute lifetime and average speed of TM tracks within each bands
    % defined by dLims
    [trackInfos1, xMap1, yMap1] = mpm2trackInfos(MPM,distToEdge,dLims,minLifetime);
    
    % Read the Actin MPM
    load(fullfile(movieData.fsmDirectory{2}, 'tack', 'mpm.mat'));
    
    % Compute lifetime and average speed of TM tracks within each bands
    % defined by dLims
    [trackInfos2, xMap2, yMap2] = mpm2trackInfos(MPM,distToEdge,dLims,minLifetime);
    
    for iBand = 1:nBands
        row1 = trackInfos1{iBand}(:,1);
        row2 = trackInfos2{iBand}(:,1);
        
        firstFrame1 = trackInfos1{iBand}(:,2);
        firstFrame2 = trackInfos2{iBand}(:,2);
        
        lastFrame1 = trackInfos1{iBand}(:,3);
        lastFrame2 = trackInfos2{iBand}(:,3);
        
        lifetime1 = trackInfos1{iBand}(:,3) - trackInfos1{iBand}(:,2) + 1;
        lifetime2 = trackInfos2{iBand}(:,3) - trackInfos2{iBand}(:,2) + 1;
        
        indFirst = sub2ind(size(xMap1), row1, firstFrame1);
        indLast = sub2ind(size(xMap1), row1, lastFrame1);
        
        avgVelocity1 = sqrt((xMap1(indFirst) - xMap1(indLast)).^2 + ...
            ((yMap1(indFirst) - yMap1(indLast)).^2)) ./ (lifetime1 - 1);
        avgVelocity1 = avgVelocity1 * pixelSize * 60 / timeInterval;
                
        indFirst = sub2ind(size(xMap2), row2, firstFrame2);
        indLast = sub2ind(size(xMap2), row2, lastFrame2);
        avgVelocity2 = sqrt((xMap2(indFirst) - xMap2(indLast)).^2 + ...
            ((yMap2(indFirst) - yMap2(indLast)).^2)) ./ (lifetime2 - 1);
        avgVelocity2 = avgVelocity2 * pixelSize * 60 / timeInterval;

        % Get all pair of indices combinations
        n1 = numel(row1);
        n2 = numel(row2);
        
        pairTracks = pcombs2(1:n1, 1:n2);        
        
        % Compute the latest first frame
        maxFirstFrame = max(firstFrame1(pairTracks(:,1)),firstFrame2(pairTracks(:,2)));
        % Compute the earliest last frame
        minLastFrame = min(lastFrame1(pairTracks(:,1)),lastFrame2(pairTracks(:,2)));
        
        % Compute the number of overlapping frames
        overlap = max(minLastFrame - maxFirstFrame + 1, 0);

        % Sorted overlap
        [sortedOverlap perms] = sort(overlap);
        ind = find(sortedOverlap,'first',1);
        sortedOverlap = sortedOverlap(ind:end);
        perms = perms(ind:end);
        
        % Find the set of overlaps
        overlapSet = unique(sortedOverlap);
        
        for i = 1:numel(overlapSet)
            ind = perms(overlapInterval(i):overlapInterval(i+1));
            
            firstFrame = maxFirstFrame(ind);
            
            rowRange1 = row1(pairTracks(ind,1));
            rowRange2 = row2(pairTracks(ind,2));
            
            ind1 = arrayfun(@(j) rowRange1 + (firstFrame + j - 2) * n1, 1:overlapSet(i));
            ind1 = ind1(:);
            
            ind2 = arrayfun(@(j) rowRange2 + (firstFrame + j - 2) * n2, 1:overlapSet(i));
            ind2 = ind2(:);
            
            xData1 = reshape(xMap1(ind1), numel(ind), iOverlap);
            yData1 = reshape(yMap1(ind1), numel(ind), iOverlap);

            xData2 = reshape(xMap2(ind2), numel(ind), iOverlap);
            yData2 = reshape(yMap2(ind2), numel(ind), iOverlap);
            
            avgPWD{i} = sqrt((mean(xData1,2) - mean(xData2,2)).^2 + ...
                (mean(yData1,2) - mean(yData2,2)).^2);    
        end
        
        
        
        
        toc
        disp([' (' num2str(nFrames) ') ' num2str(iTM) '/' num2str(iBand)]);
        
        indPair = vertcat(indPair{:});
        avgPWD = vertcat(avgPWD{:});

        th = avgPWD * pixelSize <= maxDist;
        
        costMatrix = zeros(n1,n2);
        costMatrix(:) = -1; % nonlinkmarker
        ind = sub2ind(size(costMatrix), pairTracks(indPair(th),1), pairTracks(indPair(th),2));
        costMatrix(ind) = avgPWD(th);
        
        % Associate TM and Actin tracks using LAP.
        link12 = lap(costMatrix,-1,0,1);
        
        % Save Data
        ind = find(link12(1:n1) <= n2);
        
        hFig = figure('Visible', 'off');
        set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
        set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
        scatter(avgVelocity1(ind), avgVelocity2(link12(ind)), 8, [.7 .7 .7], 'Marker', '.');
        
        set(gca,'XLim',[0 2000]);
        set(gca,'YLim',[0 2000]);
        
        xlabel(['Average ' tmNames{iTM} ' Velocity (nm/min)']);
        ylabel('Average Actin Velocity (nm/min)');
        
        fileName = [outputDirectory filesep 'Fig5_' num2str(iTM) num2str(iBand) '.eps'];
        print(hFig, '-depsc', fileName);
        fixEpsFile(fileName);
        close(hFig);
    end
end
