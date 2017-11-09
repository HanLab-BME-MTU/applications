function plotCompTrackReceptors(receptorInfoLabeled,sample,track,varargin)
% function plotCompTrack(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure,...
%     plotAS,timeStep,markMS)
%PLOTCOMPTRACKRECEPTORS plots the x-coordinates, y-coordinates and/or intensities along a compound track
%
%SYNOPSIS plotCompTrackAmp(trackedFeatureInfo,plotX,plotY,plotA,inOneFigure,...
%    plotAS,timeStep,markMS)
%
%INPUT  trackedFeatureInfo: Output of trackCloseGapsKalman for one track:
%                           Contains the fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%           .aggregationState: This field results from running the function
%                              aggregStateFromCompTracks. Only needed if
%                              aggregation state is to be plotted.
%       plotX             : 1 if x-coordinate is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       plotY             : 1 if y-coordinate is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       plotA             : 1 if amplitude is to be plotted, zero
%                           otherwise. Optional. Default: 1.
%       inOneFigure       : 1 if all plots appear in one figure window (one
%                           above the other), 0 if each figure is in its
%                           own window. Optional. Default: 1.
%       plotAS   : 1 to plot particle aggregation state (if
%                           supplied), 0 otherwise. Optional. Default: 1.
%       timeStep          : Time step for x-axis. Can be whatever units are
%                           relevant. Optional. Default: 1.
%                           NOTE: If timeStep is supplie, x-axis starts at
%                                 0. If not supplied, x-axis starts at 1.
%       markMS            : 1 to mark merges and splits, 0 otherwise.
%                           Optional. Default: 0.
%
%OUTPUT The plot(s).
%
%This is a modification of plotCompTrack.m from Khuloud Jaqaman, May 2007.
%This program plots the trajectors of the individual receptors in a
%compound track instead of the individual clusters, which is only different
%when some error has been introduced into the receptorTraj matrix and the
%traces of each individual receptor are desired.

%% Input

ip = inputParserRetrofit;
ip.addRequired('trackedFeatureInfo', ...
    @(s) isstruct(s) || isnumeric(s));
ip.addArgument('plotX',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('plotY',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('plotA',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('inOneFigure',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('plotAS',1, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.addArgument('timeStep',[], ...
    @(x) isnumeric(x) || isempty(x));
ip.addArgument('markMS',0, ...
    @(x) islogical(x) || ismember(x,[0 1]));
ip.parse(receptorInfoLabeled, varargin{:});

assignFieldsHere(ip.Results);

if(isempty(ip.Results.timeStep))
    timeStep = 1;
    axisOffset = 0;
else
    axisOffset = 1;
end

%extract information from input
seqOfEvents = receptorInfoLabeled(sample).compTracks(track).seqOfEvents;
tracksCoordAmpCG = receptorInfoLabeled(sample).compTracks(track).tracksCoordAmpCG;
receptTraj = receptorInfoLabeled(sample).receptorTraj;
tracksFeatIndxCG = receptorInfoLabeled(sample).compTracks(track).tracksFeatIndxCG;
clust2recept = receptorInfoLabeled(sample).clust2receptAssign;
[~,~,numFrames] = size(receptTraj);

if isfield(receptorInfoLabeled,'aggregState')
    aggregState = receptorInfoLabeled.aggregState;
else
    plotAS = 0;
end

%convert sparse to full if necessary
if issparse(tracksFeatIndxCG)
    tracksFeatIndxCG = full(tracksFeatIndxCG);
end

%determine whether track is sampled regularly or with doubled frequency
doubleFreq = mod(seqOfEvents(1,1)*2,2)==1;

%% Plotting

%extract receptors for the given compound track
    clustInCompTrack = tracksFeatIndxCG(tracksFeatIndxCG(:,1)~=0,1);
    receptInCompTrack = unique(clust2recept(clustInCompTrack,:,1));
    for iFrame = 2 : numFrames
        clustInCompTrack = tracksFeatIndxCG(tracksFeatIndxCG(:,iFrame)~=0,iFrame);
        receptToAdd = clust2recept(clustInCompTrack,:,iFrame);
        receptInCompTrack = unique([receptInCompTrack; receptToAdd(:)]);
    end
    receptInCompTrack = unique(receptInCompTrack(:));
    receptInCompTrack = receptInCompTrack(receptInCompTrack~=0);

%x-coordinates
if plotX

    %open new figure window and hold on to it
    figure, hold on

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAS,1,1)
        hold on
    end

    %get and plot the x-coordinates
    xCoordSequence = squeeze(receptTraj(receptInCompTrack,1,:));
    plotCompTrackCore(xCoordSequence,seqOfEvents,doubleFreq,timeStep,axisOffset);
    
    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('X-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotY && ~plotA && ~plotAS)
        hold off
    end

end %(if plotX)

%y-coordinates
if plotY

    %open new figure window if needed and hold on to it
    if ~inOneFigure || ~plotX
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAS,1,plotX+1)
        hold on
    end
    
    %get and plot the y-coordinates
    yCoordSequence = squeeze(receptTraj(receptInCompTrack,2,:));
    plotCompTrackCore(yCoordSequence,seqOfEvents,doubleFreq,timeStep,axisOffset);

    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('Y-coordinate (pixels)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotA && ~plotAS)
        hold off
    end

end %(if plotY)

%amplitudes
if plotA

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY)
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAS,1,plotX+plotY+1)
        hold on
    end

    %extract amplitudes from input
    ampSequence = tracksCoordAmpCG(:,4:8:end);

    %mark merges and splits
    if markMS
        yaxisMin = min(ampSequence(:));
        yaxisMax = max(ampSequence(:));
        plot((splitTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(splitTimes)),'k+');
        plot((mergeTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(mergeTimes)),'rx');
        legend({'Splits','Merges'})
    end
    
    %plot the amplitudes, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(ampSequence,seqOfEvents,doubleFreq,timeStep,axisOffset);

    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('Amplitude (a.u.)');

    %hold off of figure if each plot is in a separate figure window or if
    %this is the last plot
    if ~inOneFigure || (~plotAS)
        hold off
    end

end %(if plotA)

%aggregation state
if plotAS

    %open new figure window if needed and hold on to it
    if ~inOneFigure || (~plotX && ~plotY && ~plotA)
        figure, hold on
    end

    %if all figures will be in the same window, specify subplot
    if inOneFigure
        subplot(plotX+plotY+plotA+plotAS,1,plotX+plotY+plotA+1)
        hold on
    end

    %replace zeros with NaNs in aggregation state matrix
    aggregState(aggregState==0) = NaN;
    
    %mark merges and splits
    if markMS
        yaxisMin = min(aggregState(:));
        yaxisMax = max(aggregState(:));
        plot((splitTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(splitTimes)),'k+');
        plot((mergeTimes-axisOffset)*timeStep,(yaxisMin-0.1*(yaxisMax-yaxisMin))*ones(size(mergeTimes)),'rx');
        legend({'Splits','Merges'})
    end
    
    %plot the aggregation states, taking into account closed gaps, merges and
    %splits and sampling frequency doubling
    plotCompTrackCore(aggregState,seqOfEvents,doubleFreq,timeStep,axisOffset);

    %put axes labels
    if ~inOneFigure
        if axisOffset
            xlabel('Time');
        else
            xlabel('Frame number');
        end
    end
    ylabel('Aggregation state');

    %hold off of figure
    hold off
    
end %(if plotAS)

%put x-axis label at bottom in case of everything in one figure
if inOneFigure
    if axisOffset
        xlabel('Time');
    else
        xlabel('Frame number');
    end
end


%% Subfunction

function plotCompTrackCore(valuesMatrix,seqOfEvents,doubleFreq,timeStep,axisOffset)

if nargin < 4 || isempty(timeStep)
    timeStep = 1;
end
if nargin < 5 || isempty(axisOffset)
    axisOffset = 0;
end

%get first frame, last frame and number of frames
firstFrame = 1;
lastFrame = length(valuesMatrix);

%get sampling frequency
samplingFreq = 1 / (1+doubleFreq);

%get number of segments making compound track
numSegments = size(valuesMatrix,1);

%plot values as dotted black lines, closing gaps
for i = 1 : numSegments
    indx = find(~isnan(valuesMatrix(i,:)));
    plot(((indx-1)*samplingFreq+firstFrame-axisOffset)*timeStep,valuesMatrix(i,indx),'k:');
end

%plot values in color, leaving gaps as blank (so that they appear as
%dotted lines in the final figure)
plot(((firstFrame:samplingFreq:lastFrame)'-axisOffset)*timeStep,valuesMatrix','marker','.');


%% ~~~ the end ~~~

