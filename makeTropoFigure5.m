function makeTropoFigure5(paths, outputDirectory)

outputDirectory = fullfile(outputDirectory, 'Fig5');

if ~exist(outputDirectory,'dir')
    mkdir(outputDirectory);
end

nBands = 3;
dLims = [0 2 4 6] * 1000;
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
        maxFirstFrame = max(firstFrame1(pairTracks(:,1)),...
            firstFrame2(pairTracks(:,2)));
        % Compute the earliest last frame
        minLastFrame = min(lastFrame1(pairTracks(:,1)),...
            lastFrame2(pairTracks(:,2)));
        
        % Compute the number of overlapping frames
        overlap = max(minLastFrame - maxFirstFrame + 1, 0);
        
        clear minLastFrame;
        
        % Sorted overlap
        [sortedOverlap indPair] = sort(overlap);
        
        clear overlap;
        
        ind = find(sortedOverlap >= 3,1,'first');
        sortedOverlap = sortedOverlap(ind:end);
        indPair = indPair(ind:end);
        
        % Find the set of overlaps
        overlapSet = unique(sortedOverlap);
        
        % Find the indices in sortedOverlap where an overlap transition
        % occurs
        overlapInterval = [0; find(sortedOverlap(2:end) ~= ...
            sortedOverlap(1:end-1)); numel(sortedOverlap)] + 1;
        
        avgPWD = cell(numel(overlapSet), 1);
        
        for i = 1:numel(overlapSet)
            ind = indPair(overlapInterval(i):overlapInterval(i+1)-1);
            
            firstFrame = maxFirstFrame(ind);
            
            rowRange1 = row1(pairTracks(ind,1));
            rowRange2 = row2(pairTracks(ind,2));
            
            ind1 = arrayfun(@(j) rowRange1 + (firstFrame + j - 2) * size(xMap1,1),...
                1:overlapSet(i),'UniformOutput',false);
            ind1 = horzcat(ind1{:});
            
            ind2 = arrayfun(@(j) rowRange2 + (firstFrame + j - 2) * size(xMap2,1),...
                1:overlapSet(i),'UniformOutput',false);
            ind2 = horzcat(ind2{:});
            
            xData1 = reshape(xMap1(ind1), numel(ind), overlapSet(i));
            yData1 = reshape(yMap1(ind1), numel(ind), overlapSet(i));

            xData2 = reshape(xMap2(ind2), numel(ind), overlapSet(i));
            yData2 = reshape(yMap2(ind2), numel(ind), overlapSet(i));
            
            avgPWD{i} = sqrt((mean(xData1,2) - mean(xData2,2)).^2 + ...
                (mean(yData1,2) - mean(yData2,2)).^2);    
        end
        
        avgPWD = vertcat(avgPWD{:}) * pixelSize;

        th = avgPWD <= maxDist;
        
        indPair = indPair(th);
        avgPWD = avgPWD(th);
        
        % Associate TM and Actin tracks using LAP.
        avgPWD(avgPWD == 0) = eps;
        costMatrix2 = sparse(pairTracks(indPair,1), pairTracks(indPair,2), avgPWD, n1, n2);
        link12 = lap(costMatrix2,[],0,1);
        
        % Save Data
        ind = find(link12(1:n1) <= n2);
        X = avgVelocity1(ind);
        Y = avgVelocity2(link12(ind));
        covXY = cov(X,Y);
        r = covXY(2) / (std(X) * std(Y));
        hFig = figure('Visible', 'off');
        set(gca, 'FontName', 'Helvetica', 'FontSize', 18);
        set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
        scatter(X, Y, 10, [.7 .7 .7], 'Marker', '.');
        text(1000,1800,['r = ' num2str(r)],'HorizontalAlignment','center', 'FontSize', 20);
                
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
