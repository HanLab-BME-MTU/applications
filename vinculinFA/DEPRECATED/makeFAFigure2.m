function makeFAFigure2(paths, outputDirectory)
% Compute the mean speed of tracks for 2 different conditions. First path
% needs to be the control, 2nd path is the -/- condition.

speeds = cell(2,1);

for iMovie = 1:2
    % Load movieData
    fileName = fullfile(paths{iMovie},'movieData.mat');
    if ~exist(fileName, 'file')
        error(['Unable to locate ' fileName]);
    end
    load(fileName);

    % Check that detection has been performed
    if ~checkMovieDetection(movieData) 
        error('Need detections to make Figure 1.');
    end

    % Check that tracking has been performed
    if ~checkMovieTracking(movieData)
        error('Need tracks to make Figure 1.');
    end

    % Load detections
    load(fullfile(movieData.detection.directory, movieData.detection.filename));

    nFrames = numel(segmentParams); 

    % Load tracks
    load(fullfile(movieData.tracking.directory, movieData.tracking.filename));

    % Remove tracks which either
    % - lifetime < 2,
    % - starts at first frame
    % - ends at last frame

    trackSEL = getTrackSEL(tracksFinal); 
    ind = trackSEL(:,3) >= 2 | trackSEL(:,1) ~= 1 | trackSEL(:,2) ~= nFrames;
    tracksFinal = tracksFinal(ind);
    trackSEL = trackSEL(ind,:);

    % Compute translocation (distance between FA position at birth and death)
    nTracks = numel(tracksFinal);
    
    %get number of segments making each track
    numSegments = zeros(nTracks,1);
    for i = 1 : nTracks
        numSegments(i) = size(tracksFinal(i).tracksCoordAmpCG,1);
    end

    %locate the row of the first segment of each compound track in the
    %big matrices of all tracks (to be constructed in the next step)
    trackStartRow = ones(nTracks,1);
    for iTrack = 2:nTracks
        trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
    end

    %find total number of segments in all tracks (i.e. number of rows in big
    %matrices)
    numSegmentsTracks = trackStartRow(end)+numSegments(end)-1;

    xCoord = zeros(numSegmentsTracks,nFrames);
    yCoord = xCoord;
    
    for iTrack = 1:nTracks
        %get track start and end times
        startTime = trackSEL(iTrack,1);
        endTime   = trackSEL(iTrack,2);
            
        %store x-coordinates
        xCoord(trackStartRow(iTrack):trackStartRow(iTrack)+...
            numSegments(iTrack)-1,startTime:endTime) = ...
            tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
    
        %store y-coordinates
        yCoord(trackStartRow(iTrack):trackStartRow(iTrack)+...
            numSegments(iTrack)-1,startTime:endTime) = ...
            tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
    end
    
    xCoord = arrayfun(@(iSegment) nonzeros(xCoord(iSegment,:)), ...
        (1:numSegmentsTracks)', 'UniformOutput', false);

    yCoord = arrayfun(@(iSegment) nonzeros(yCoord(iSegment,:)), ...
        (1:numSegmentsTracks)', 'UniformOutput', false);
    
    d = arrayfun(@(iSegment) sqrt((xCoord{iSegment}(1:end-1) - ...
        xCoord{iSegment}(2:end)).^2 + (yCoord{iSegment}(1:end-1) - ...
        yCoord{iSegment}(2:end)).^2), (1:numSegmentsTracks)', ...
        'UniformOutput', false);
    
    % Convert to nanometers
    speeds{iMovie} = vertcat(d{:}) * movieData.pixelSize_nm / movieData.timeInterval_s;
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
ylabel('Speed (nm/s)');

fileName = [outputDirectory filesep 'Fig2-1.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);text
close(hFig);

