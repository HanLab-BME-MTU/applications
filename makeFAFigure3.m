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
    sigmaPSF = 1.6255;

    % Check that detection has been performed
    if ~checkMovieSegmentDetection(movieData) 
        error('Need detections to make Figure 3.');
    end

    % Load detections
    load(fullfile(movieData.segmentDetection.directory, movieData.segmentDetection.filename));

    % Check that Actin tracking has been done
    filename = fullfile(movieData.imageDirectory, movieData.channelDirectory{2}, '..', 'analysis', 'tack', 'mpm.mat');
    
    if ~exist(filename, 'file')
        error('Unable to find mpm.mat');
    end
    
    % Load MPM.mat
    load(filename);

    % Read mask file list
    %maskPath = fullfile(movieData.masks.directory, movieData.masks.channelDirectory{1});
    %maskFiles = dir([maskPath filesep '*.tif']);
    
    nFrames = size(M,3); %#ok<NODEF>
    
    % Get the indices where there are vectors in M
    nnzInd = all(M(:,:,:),2);
    
    length = zeros(size(nnzInd));
    
    % Generate FA's footprint
    for iFrame = 1:nFrames
        lengthMap = zeros(imSize(1), imSize(2));
        nMap = zeros(imSize(1), imSize(2));
        
        segments = segmentParams{iFrame}; %#ok<USENS>
        
        nSegments = size(segments, 1);
        
        for iSegment = 1:nSegments
            xC = segments(iSegment,1);
            yC = segments(iSegment,2);
            l = segments(iSegment,4);
            theta = segments(iSegment,5);
        
            [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigmaPSF,l,theta,imSize);
        
            S = zeros(numel(yRange),numel(xRange));
            S(nzIdx) = l;
            
            lengthMap(yRange,xRange) = lengthMap(yRange,xRange) + S;
            nMap(yRange,xRange) = nMap(yRange,xRange) + (S ~= 0);
        end
        
        ind = find(nMap);
        lengthMap(ind) = lengthMap(ind) ./ nMap(ind);
        
        %BW = imread(fullfile(maskPath, maskFiles(iFrame).name));
        %dist = bwdist(1 - double(BW)) * pixelSize;
        %lengthMap(dist > 5000) = 0;
        
        ind = nnzInd(:,:,iFrame);
        
        I = M(ind,1,iFrame);
        J = M(ind,2,iFrame);
        
        locMax = sub2ind(imSize, I, J);
        
        length(ind,:,iFrame) = lengthMap(locMax);        
    end

    % Compute the speed
    D = cell(nFrames,1);
    L = cell(nFrames,1);
    for iFrame = 1:nFrames
        ind = length(:,:,iFrame) ~= 0;
        
        D{iFrame} = sqrt((M(ind,1,iFrame) - M(ind,3,iFrame)).^2 + ...
            (M(ind,2,iFrame) - M(ind,4,iFrame)).^2);
        
        L{iFrame} = length(ind,:,iFrame);
    end
    
    speeds =  vertcat(D{:}) * pixelSize / timeInterval;
    lengths = vertcat(L{:}) * pixelSize;
    
    lengthRange = min(lengths):1000:max(lengths);
  
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

hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
h = bar(gca, avgSpeeds, 'group'); hold on;
XTicks = get(h(1), 'XData');
XTicks = XTicks - .5;
set(h(1), 'XData', XTicks);
set(h(2), 'XData', XTicks);
errorbar([XTicks-.15 XTicks+.15],avgSpeeds(:),zeros(2*numel(XTicks),1),stdSpeeds(:),'xk'); hold off;
xlabel('Adhesion length (nm)');
ylabel('Actin Speed (nm/s)');
legend('fl/fl', 'WT');

fileName = [outputDirectory filesep 'Fig3-1.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);text
close(hFig);


