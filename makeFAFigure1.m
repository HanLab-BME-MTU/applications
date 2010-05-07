function makeFAFigure1(paths, outputDirectory)
% Compute the mean translocation of tracks for 2 different conditions.
% First path needs to be the control, 2nd path is the -/- condition.

transLoc = cell(2,1);

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

    nFrames = numel(segmentParams); %#ok<USENS>

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
    
    transLoc{iMovie} = zeros(nTracks,1);
    
    for iTrack = 1:nTracks
        %get track start and end times
        startTime = trackSEL(iTrack,1);
        endTime   = trackSEL(iTrack,2);
        
        x1 = segmentParams{startTime}(tracksFinal(iTrack).tracksFeatIndxCG(1),1);
        x2 = segmentParams{endTime}(tracksFinal(iTrack).tracksFeatIndxCG(end),1);
        
        y1 = segmentParams{startTime}(tracksFinal(iTrack).tracksFeatIndxCG(1),2);
        y2 = segmentParams{endTime}(tracksFinal(iTrack).tracksFeatIndxCG(end),2);
        
        transLoc{iMovie}(iTrack) = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    end
    
    % Convert to nanometers
    transLoc{iMovie} = transLoc{iMovie} * movieData.pixelSize_nm;
end

avgTransLoc = cellfun(@mean,transLoc);
stdTransLoc = cellfun(@std,transLoc);

hFig = figure('Visible', 'off');
set(gca, 'FontName', 'Helvetica', 'FontSize', 20);
set(gcf, 'Position', [680 678 560 400], 'PaperPositionMode', 'auto');
h = bar(gca, avgTransLoc, 'group'); hold on;
% Get the central position over the bars
XTicks = get(h, 'XData');
errorbar(XTicks,avgTransLoc(:),zeros(size(XTicks)),stdTransLoc(:),'xk'); hold off;
set(gca, 'XTickLabel', {'Vcl fl/fl', 'Vcl -/-'});
ylabel('Translocation (nm)');

fileName = [outputDirectory filesep 'Fig1-2.eps'];
print(hFig, '-depsc', fileName);
fixEpsFile(fileName);text
close(hFig);

