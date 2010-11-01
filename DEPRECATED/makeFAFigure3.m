function makeFAFigure3(paths, outputDirectory)
% Compute the average Actin speed under focal adhesion footprints.

% TODO: we still need to remove 1-frame long tracks

avgSpeeds = cell(2, 1);
stdSpeeds = cell(2, 1);

for iMovie = 1:2
    % Load movieData
    fileName = fullfile(paths{iMovie},'movieData.mat');
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    % Get image size
    imSize = movieData.imSize;
    pixelSize = movieData.pixelSize_nm;
    timeInterval = movieData.timeInterval_s;
    nFrames = movieData.nImages;

    % Check that detection has been performed
    if ~checkMovieSegmentDetection(movieData) 
        error('Need detections to make Figure 3.');
    end

    % Load detections
    load(fullfile(movieData.segmentDetection.directory, movieData.segmentDetection.filename));

    % Check that Actin tracking has been done
    filename = fullfile(movieData.imageDirectory, movieData.channelDirectory{2}, '..', 'analysis', 'tack', 'mpm.mat');
    
    lengthMap = zeros(imSize(1), imSize(2), nFrames);

    % Generate FA's footprint
    for iFrame = 1:nFrames
        nMap = zeros(imSize(1), imSize(2));
        
        clusters = segmentParams{iFrame}; %#ok<USENS>
        segments = vertcat(clusters(:).segments);
        
        nSegments = size(segments, 1);

        % TODO: discard segment that are outside the first 3um and are
        % smaller than 0.5um.
        
        for iSegment = 1:nSegments
            xC = segments(iSegment,1);
            yC = segments(iSegment,2);
            sigmaPSF = segments(iSegment,4);
            l = segments(iSegment,5);
            theta = segments(iSegment,6);
        
            [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigmaPSF,l,theta,imSize);
        
            S = zeros(numel(yRange),numel(xRange));
            S(nzIdx) = l;
            
            lengthMap(yRange,xRange,iFrame) = lengthMap(yRange,xRange,iFrame) + S;
            nMap(yRange,xRange) = nMap(yRange,xRange) + (S ~= 0);
        end
        
        ind = find(nMap);
        offset = (iFrame - 1) * imSize(1) * imSize(2);
        lengthMap(ind + offset) = lengthMap(ind + offset) ./ nMap(ind);       
    end

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
    
    if ~exist(filename, 'file')
        error('Unable to find mpm.mat');
    end
    
    % Load MPM.mat
    load(filename);    
    
    % Compute Actin speed

    % TODO: One should vary minimum lifetime to 1 or 2.
    [trackInfos, xMap, yMap] = mpm2trackInfos(MPM,distToEdge,[0 +Inf],3);
    
    row = trackInfos{1}(:,1);
    
    firstFrame = trackInfos{1}(:,2);
    
    lastFrame = trackInfos{1}(:,3);
                
    lifetime = trackInfos{1}(:,3) - trackInfos{1}(:,2) + 1;
    
    indFirst = sub2ind(size(xMap), row, firstFrame);
    indLast = sub2ind(size(xMap), row, lastFrame);
    
    speeds = sqrt((xMap(indFirst) - xMap(indLast)).^2 + ...
        ((yMap(indFirst) - yMap(indLast)).^2)) ./ (lifetime - 1);
    speeds = speeds * pixelSize * 60 / timeInterval;
    
    % Compute the length
    lengths = zeros(size(speeds));
    
    sortedLifetime = sort(lifetime);
    sortedLifetime = unique(sortedLifetime);
   
   for i = 1:numel(sortedLifetime)
       
       l = sortedLifetime(i);
       
       ind = find(lifetime == l);
       
       frames = cell2mat(arrayfun(@(j) firstFrame(ind) + j - 1, 1:l, 'UniformOutput', false));
       rowMap = repmat(row(ind),1,l);
       
       I = yMap(sub2ind(size(yMap),rowMap(:), frames(:)));
       J = xMap(sub2ind(size(xMap),rowMap(:), frames(:)));
       
       tmp = reshape(lengthMap(sub2ind(size(lengthMap), ...
           I(:), J(:), frames(:))), numel(ind), l);
       
       nnzFrames = sum(tmp ~= 0, 2);
       
       ind2 = find(nnzFrames);
       
       if ~isempty(ind2)
           lengths(ind(ind2)) = mean(tmp(ind2,:), 2) * l ./ nnzFrames(ind2);
       end
   end
   
   % Remove zeros or NaN elements
   ind = lengths ~= 0 & ~isnan(lengths);
   lengths = lengths(ind);
   speeds = speeds(ind);
   
   lengths = lengths * pixelSize;
   
   lengthRange = [0 250 500 1000:1000:8000];
   
   avgSpeeds{iMovie} = zeros(numel(lengthRange)-1,1);
   stdSpeeds{iMovie} = zeros(numel(lengthRange)-1,1);
   
   for iLength = 1:numel(lengthRange)-1
       ind = lengths >= lengthRange(iLength) & lengths < lengthRange(iLength+1);
       
       avgSpeeds{iMovie}(iLength) = mean(speeds(ind));
       stdSpeeds{iMovie}(iLength) = std(speeds(ind));
   end
end

n1 = numel(avgSpeeds{1});
n2 = numel(avgSpeeds{2});

if n1 ~= n2
    [r, i] = max([n1,n2]);
    avgSpeeds{2-i+1}(end+1:r) = NaN;
    stdSpeeds{2-i+1}(end+1:r) = NaN;
end

avgSpeeds = horzcat(avgSpeeds{:});
stdSpeeds = horzcat(stdSpeeds{:});

colors = [
   0.290000000000000   0.617000000000000   0.547000000000000;
   0.360000000000000   0.630000000000000   0.900000000000000;
   0.060000000000000   0.330000000000000   0.600000000000000;
   0.000000000000000   0.000000000000000   0.000000000000000
   0.000000000000000   0.000000000000000   0.000000000000000
   0.000000000000000   0.000000000000000   0.000000000000000];

hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
h = bar(gca, avgSpeeds, 'group'); hold on;
XTicks = get(h(1), 'XData');
XTicks = XTicks - .5;
set(h(1), 'XData', XTicks);
set(h(2), 'XData', XTicks);
set(get(h(1),'Children'), 'FaceColor',  colors(1,:));
set(get(h(2),'Children'), 'FaceColor',  colors(2,:));

set(gca,'XTickLabel',[]);
arrayfun(@(x) text(x, -.8/800,num2str(lengthRange(x+1)),...
    'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right', 'Rotation', 45,...
    'FontSize', 12,'Parent',gca), 1:numel(lengthRange)-1);

errorbar([XTicks-.15 XTicks+.15],avgSpeeds(:),zeros(2*numel(XTicks),1),stdSpeeds(:),'xk'); hold off;
ylabel('Actin Speed (nm/min)');
h = legend('con', 'cre'); legend('boxoff');

hC = get(h,'Children');
hC = hC(1:2:end);
for i=1:numel(hC)
    hCC = get(hC(i),'Children');
    set(hCC,'FaceColor', colors(numel(hC) - i + 1,:));
end

fileName = [outputDirectory filesep 'FA length - Actin speed.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);text
close(hFig);


