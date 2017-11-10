function plotAdaptiveWindowsProtEvent(eventWindows,winPositions,frameOnset,frameStart,frameEnd)

%% Collect windows

%get number of slices
numSlice = length(eventWindows);

%static windows
winsOnsetFrame = winPositions(frameOnset,:);
for iSlice = 1 : numSlice
    indxSlice = eventWindows(iSlice).onset{1,1}(1,2);
    winsOnset{iSlice} = winsOnsetFrame{indxSlice}(eventWindows(iSlice).onset{1,1}(:,1));
    winsBand2{iSlice} = winsOnsetFrame{indxSlice}(eventWindows(iSlice).onset{1,2}(:,1));
    winsBand3{iSlice} = winsOnsetFrame{indxSlice}(eventWindows(iSlice).onset{1,3}(:,1));
    tmp = [];
    for iBigBand = 4 : 10
        if ~isempty(eventWindows(iSlice).onset{1,iBigBand})
            tmp = [tmp; eventWindows(iSlice).onset{1,iBigBand}(:,1)];
        end
    end
    winsLamella{iSlice} = winsOnsetFrame{indxSlice}(tmp);
end

%dynamic windows before
for iFrame = frameStart : frameOnset-1
    winsCurrentFrame = winPositions(iFrame,:);
    frameDiff = frameOnset - iFrame;
    for iSlice = 1 : numSlice
        indxSlice = eventWindows(iSlice).onset{1,1}(1,2);
        winsBefDynamic{iSlice,frameDiff} = winsCurrentFrame{indxSlice}(eventWindows(iSlice).befDynamic{frameDiff,1}(:,1));
    end
end

%dynamic windows after
for iFrame = frameOnset + 1 : frameEnd
    winsCurrentFrame = winPositions(iFrame,:);
    frameDiff = iFrame - frameOnset;
    for iSlice = 1 : numSlice
        indxSlice = eventWindows(iSlice).onset{1,1}(1,2);
        winsAftDynamic{iSlice,frameDiff} = winsCurrentFrame{indxSlice}(eventWindows(iSlice).aftDynamic{frameDiff,frameDiff}(:,1));
    end    
end
 
%% Plot windows

%before protrusion onset
for iFrame = frameStart : frameOnset-1
    frameDiff = frameOnset - iFrame;
    %right behind edge
    plotWindows([winsBefDynamic{:,frameDiff}],{'c','FaceAlpha',0.3,'EdgeColor','c'});
    hold on
    %onset
    plotWindows(winsOnset,{'y','FaceAlpha',0.3,'EdgeColor','y'});
    %band 2
    plotWindows(winsBand2,{'b','FaceAlpha',0.3,'EdgeColor','b'});
    %band 3
    plotWindows(winsBand3,{'g','FaceAlpha',0.3,'EdgeColor','g'});
    %lamella
    plotWindows(winsLamella,{'m','FaceAlpha',0.3,'EdgeColor','m'});
end

%at protrusion onset
%right behind edge
plotWindows(winsOnset,{'c','FaceAlpha',0.3,'EdgeColor','c'});
hold on
%onset
plotWindows(winsOnset,{'y','FaceAlpha',0.3,'EdgeColor','y'});
%band 2
plotWindows(winsBand2,{'b','FaceAlpha',0.3,'EdgeColor','b'});
%band 3
plotWindows(winsBand3,{'g','FaceAlpha',0.3,'EdgeColor','g'});
%lamella
plotWindows(winsLamella,{'m','FaceAlpha',0.3,'EdgeColor','m'});

%during protrusion
for iFrame = frameOnset+1 : frameEnd
    frameDiff = iFrame - frameOnset;
    %new protrusion area
    plotWindows([winsAftDynamic{:,1:frameDiff}],{'r','FaceAlpha',0.3,'EdgeColor','r'});    
    hold on
    %right behind edge
    plotWindows([winsAftDynamic{:,frameDiff}],{'c','FaceAlpha',0.3,'EdgeColor','c'});
    %onset
    plotWindows(winsOnset,{'y','FaceAlpha',0.3,'EdgeColor','y'});
    %band 2
    plotWindows(winsBand2,{'b','FaceAlpha',0.3,'EdgeColor','b'});
    %band 3
    plotWindows(winsBand3,{'g','FaceAlpha',0.3,'EdgeColor','g'});
    %lamella
    plotWindows(winsLamella,{'m','FaceAlpha',0.3,'EdgeColor','m'});
end


