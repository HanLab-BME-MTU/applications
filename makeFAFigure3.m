function makeFAFigure3(paths, outputDirectory)
% Compute the average Actin speed under focal adhesion footprints.

% TODO: we still need to remove 1-frame long tracks

speeds = cell(2,1);

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
    if ~checkMovieDetection(movieData) 
        error('Need detections to make Figure 1.');
    end

    % Load detections
    load(fullfile(movieData.detection.directory, movieData.detection.filename));

    % Check that Actin tracking has been done
    filename = fullfile(movieData.channels(2).analysisDirectory, 'tack', 'mpm.mat');
    
    if ~exist(filename, 'file')
        error('Unable to find mpm.mat');
    end
    
    % Load MPM.mat
    load(filename);

    nFrames = size(M,3); %#ok<NODEF>
    
    % Get the indices where there are vectors in M
    nnzInd = all(M(:,:,:),2);
    
    % Get the indices where speckles are inside FA footprints.
    insideInd = false(size(nnzInd));
    
    % Generate FA's footprint
    for iFrame = 1:nFrames
        BW = false(imSize(1), imSize(2));
        
        segments = segmentParams{iFrame}; %#ok<USENS>
        
        nSegments = size(segments, 1);
        
        for iSegment = 1:nSegments
            xC = segments(iSegment,1);
            yC = segments(iSegment,2);
            l = segments(iSegment,4);
            theta = segments(iSegment,5);
        
            [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigmaPSF,l,theta,imSize);
        
            S = false(numel(yRange),numel(xRange));
            S(nzIdx) = true;
            
            BW(yRange,xRange) = BW(yRange,xRange) | S;
        end
        
        ind = nnzInd(:,:,iFrame);
        
        I = M(ind,1,iFrame);
        J = M(ind,2,iFrame);
        
        locMax = sub2ind(imSize, I, J);
        
        insideInd(ind,:,iFrame) = BW(locMax);        
    end

    % Compute the speed
    D = cell(nFrames,1);
    for iFrame = 1:nFrames
        ind = insideInd(:,:,iFrame);
        
        D{iFrame} = sqrt((M(ind,1,iFrame) - M(ind,3,iFrame)).^2 + ...
            (M(ind,2,iFrame) - M(ind,4,iFrame)).^2);
    end
    
    speeds{iMovie} =  vertcat(D{:}) * pixelSize / timeInterval;
end

avgSpeeds = cellfun(@mean,speeds);
stdSpeeds = cellfun(@std,speeds);

hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
h = bar(gca, avgSpeeds, 'group'); hold on;
% Get the central position over the bars
XTicks = get(h, 'XData');
errorbar(XTicks,avgSpeeds(:),zeros(size(XTicks)),stdSpeeds(:),'xk'); hold off;
set(gca, 'XTickLabel', {'Vcl fl/fl', 'Vcl -/-'});
ylabel('Actin Speed (nm/s)');

fileName = [outputDirectory filesep 'Fig3-1.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);text
close(hFig);

