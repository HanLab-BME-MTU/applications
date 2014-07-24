function [trackFront, speed] = filterUtrackSegmentedMask(pathMD)
% filterUtrackSegmentedMask filters tracks with a mask that is
% modified from segmented mask (e.g. showing only anterior part of a cell)
% 
% input:    pathMD:    path to the movieData file
%               erodeDis:   distance of erosion from the cell edge.
%                               This is usually in the case of flows from
%                               non-speclized cell where flows at the very
%                               edge is not accurate.
% output:   flowMagAllCh: speed of all tracks in the mask
% the function also stores filtered utracks mat file to analysis folder after
% backing up the original one
% speed is in the unit of nm/min.

% Load movieData
% movieDataPath = [pathMD '/movieData.mat'];
movieData = MovieData.load(pathMD);%movieDataPath);
nFrames = movieData.nFrames_;
nChannels = length(movieData.channels_);

% Load flow vectors
iTracking = movieData.getProcessIndex('TrackingProcess');
if isempty(iTracking)
    error('uTrack has to be run.')
end
trackProcess = movieData.getProcess(iTracking);
iChan = 1; % assuming only one channel is used
minLifetime = 10;
tracksOrg = trackProcess.loadChannelOutput(iChan);
% if ii>minLifetime
SEL = getTrackSEL(tracksOrg); %SEL: StartEndLifetime
% Remove any less than 3-frame long track.
isValid = SEL(:,3) >= minLifetime;
tracksOrg = tracksOrg(isValid);

detectedPointProc = movieData.getProcess(movieData.getProcessIndex('PointSourceDetectionProcess'));
detectedPoints = detectedPointProc.loadChannelOutput(iChan);

disp('reformating NA tracks...')
tracks = formatTracks(tracksOrg,detectedPoints,nFrames); 
    
% % Load segmented masks
% iMask = movieData.getProcessIndex('MaskRefinementProcess');
% maskProcess = movieData.getProcess(iMask);

% Backup the original vectors to backup folder
display('Backing up the original data')
[pathstr,~,~] = fileparts(trackProcess.outFilePaths_{1});
[pathstr1,analysis,~] = fileparts(pathstr);
backupFolder = fullfile(pathstr1, [analysis ' Backup']); % name]);
if ~exist(backupFolder,'dir')
    mkdir(backupFolder);
    copyfile(pathstr, backupFolder)
end

channelName = @(x)movieData.getChannelPaths{x}(max(regexp(movieData.getChannelPaths{x},filesep))+1:end);   
for k = 1:nChannels    
    %Create string for current directory
    [pathstr,~,~] = fileparts(trackProcess.outFilePaths_{k});
    outputDir{k} = fullfile(pathstr,channelName(k));
    mkClrDir(outputDir{k});
    outFileTif_frames=@(chan) [outputDir{chan} filesep 'tiffs'];
    outFileFig_frames=@(chan) [outputDir{chan} filesep 'figs'];
    mkClrDir(outFileTif_frames(k));
    mkClrDir(outFileFig_frames(k));
end
%Format string for zero-padding file names
outFileTracks=@(chan) [outputDir{chan} filesep 'tracksFront.mat'];
outFileFlow=@(chan) [outputDir{chan} filesep 'flow.mat'];
outFileSpeed=@(chan) [outputDir{chan} filesep 'speed.mat'];
outFileSpeedXls=@(chan) [outputDir{chan} filesep 'speed.xls'];
outFileFig_alltracks=@(chan) [outputDir{chan} filesep 'all_tracks.fig'];
outFileFig_FrontTracks=@(chan) [outputDir{chan} filesep 'frontTracks.fig'];
outFileTif_alltracks=@(chan) [outputDir{chan} filesep 'all_tracks.tif'];
outFileTif_FrontTracks=@(chan) [outputDir{chan} filesep 'frontTracks.tif'];
% mask = maskProcess.loadChannelOutput(1, 1);
% [nrow,ncol] = size(mask);
% roi = false(nrow,ncol);
iiformat = ['%.' '3' 'd'];
trackIdx = false(numel(tracks),1);
mag = 30; % magnification factor for quiver plot
for iChan=1:nChannels
    % show static tracks and cell
    h=figure;
    pax1 = movieData.getChannel(iChan).loadImage(1);
    imshow(pax1,[]);
    hold on
    for jj=1:numel(tracks)
        plot(tracks(jj).xCoord,tracks(jj).yCoord)
    end
    % Get user input by impoly that will be excluded
    display('Please draw a polygon where you want to show the vector field. Finish with right click after a final point')
    hpoly = impoly;
    polyMask = createMask(hpoly);
    delete(hpoly);

    % Apply the mask to the tracks
    for k=1:numel(tracks)
        for ii=1:nFrames
            if ~isnan(tracks(k).yCoord(ii)) && polyMask(round(tracks(k).yCoord(ii)),round(tracks(k).xCoord(ii)))
                trackIdx(k) = true;
            end
        end
    end
    tracksFront = tracks(trackIdx);
    for jj=1:numel(tracksFront)
        plot(tracksFront(jj).xCoord,tracksFront(jj).yCoord,'k','Linewidth',0.5)
    end
    % For selected tracks, perform linear regression to find average
    % distances and velocities
    speed = zeros(numel(tracksFront),1);
    timeFrame = zeros(numel(tracksFront),1);
    flow = zeros(numel(tracksFront),4);

    for k=1:numel(tracksFront)
        iStart = tracksFront(k).startingFrame;
        iEnd = tracksFront(k).endingFrame;
        % For tracks that are long in y-direction, I'll flip the
        % regression between x and y
        % Length comparison
%             xlength = nanmax(tracksFront(k).xCoord(iStart:iEnd))-nanmin(tracksFront(k).xCoord(iStart:iEnd));
%             ylength = nanmax(tracksFront(k).yCoord(iStart:iEnd))-nanmin(tracksFront(k).yCoord(iStart:iEnd));
        xlength = abs(tracksFront(k).xCoord(iEnd)-tracksFront(k).xCoord(iStart));
        ylength = abs(tracksFront(k).yCoord(iEnd)-tracksFront(k).yCoord(iStart));
        if ylength>xlength
            beta = fit(tracksFront(k).yCoord(iStart:iEnd)',tracksFront(k).xCoord(iStart:iEnd)','poly1');
            yStart = tracksFront(k).yCoord(iStart);
            xStart = feval(beta,yStart);
            yEnd = tracksFront(k).yCoord(iEnd);
            xEnd = feval(beta,yEnd);
        else
            beta = fit(tracksFront(k).xCoord(iStart:iEnd)',tracksFront(k).yCoord(iStart:iEnd)','poly1');
            xStart = tracksFront(k).xCoord(iStart);
            yStart = feval(beta,xStart);
            xEnd = tracksFront(k).xCoord(iEnd);
            yEnd = feval(beta,xEnd);
        end
        % Get the first and last point of the fit
        speed(k) = norm([xEnd-xStart, yEnd-yStart])/(iEnd-iStart)*movieData.pixelSize_/(movieData.timeInterval_/60); % nm/min
        flow(k,:)=[xStart yStart xEnd yEnd];
        timeFrame(k)=iEnd-iStart;
        trackFront(k).fit_flow = flow(k,:);
        trackFront(k).fit_speed = speed(k); % nm/min
        trackFront(k).fit_time = timeFrame(k); % nm/min
        % Show the average distance vectors
        quiver(xStart,yStart, xEnd-xStart, yEnd-yStart,0,'m', 'Linewidth',0.5);
    end
    % scale arrow (distance)
    line([5 5+5000/movieData.pixelSize_],[5 5], 'LineWidth',2,'Color','k'); % 5000 nm
    % save velocity vectors: % more detailed orientation analys can be
    % done... but later...
    save(outFileTracks(iChan),'trackFront');
    save(outFileFlow(iChan),'flow');
    save(outFileSpeed(iChan),'speed');
    xlswrite(outFileSpeedXls(iChan),speed);
%     print('-depsc2', '-r150', outFileFig_alltracks(iChan),'.eps'));
    print('-dtiff', '-r300', outFileTif_alltracks(iChan));
    hgsave(h,outFileFig_alltracks(iChan),'-v7.3')
    close(h)
    % Show only fitered tracks
    h2=figure;
    imshow(pax1,[]);
    hold on
    for jj=1:numel(tracksFront)
        plot(tracksFront(jj).xCoord,tracksFront(jj).yCoord,'k','Linewidth',0.5)
    end
    quiver(flow(:,1),flow(:,2), mag*(flow(:,3)-flow(:,1))./timeFrame, mag*(flow(:,4)-flow(:,2))./timeFrame,0,'r', 'Linewidth',1);
    % scale arrow (distance)
    line([5 5+5000/movieData.pixelSize_],[5 5], 'LineWidth',2,'Color','k'); % 5000 nm
    % scale arrow (velocity)
    quiver(movieData.imSize_(2)-5- mag*1000/movieData.pixelSize_/60*movieData.timeInterval_ ,5, mag*1000/movieData.pixelSize_/60*movieData.timeInterval_, 0,0,'k', 'Linewidth',0.5); % 1000 nm/min
    print('-dtiff', '-r300', outFileTif_FrontTracks(iChan));
    hgsave(h2,outFileFig_FrontTracks(iChan),'-v7.3')
    close(h2)
    % Creating movie with tracks
    for ii=1:nFrames
        h3=figure;
        pax = movieData.getChannel(iChan).loadImage(ii);
        imshow(pax,[]);
        hold on
        for jj=1:numel(tracksFront)
            if ii<=tracksFront(jj).endingFrame
                plot(tracksFront(jj).xCoord(1:ii),tracksFront(jj).yCoord(1:ii),'k')
            end
            if ii == tracksFront(jj).endingFrame
                quiver(trackFront(jj).fit_flow(1),trackFront(jj).fit_flow(2), mag*(trackFront(jj).fit_flow(3)-trackFront(jj).fit_flow(1))/trackFront(jj).fit_time, mag*(trackFront(jj).fit_flow(4)-trackFront(jj).fit_flow(2))/trackFront(jj).fit_time,0,'r', 'Linewidth',1);
            end
        end
        % scale arrow (distance)
        line([5 5+5000/movieData.pixelSize_],[5 5], 'LineWidth',2,'Color','k'); % 5000 nm
        % scale arrow (velocity)
        quiver(movieData.imSize_(2)-5- mag*1000/movieData.pixelSize_/60*movieData.timeInterval_ ,5, mag*1000/movieData.pixelSize_/60*movieData.timeInterval_, 0,0,'k', 'Linewidth',0.5); % 1000 nm/min
        print('-dtiff', '-r300', strcat(outFileTif_frames(iChan),'/frame',num2str(ii,iiformat),'.tif'));
        hgsave(h3,strcat(outFileFig_frames(iChan),'/frame',num2str(ii,iiformat),'.fig'),'-v7.3')
        close(h3)
    end
end
