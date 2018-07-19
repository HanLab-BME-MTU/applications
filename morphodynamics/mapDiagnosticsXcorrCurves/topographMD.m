function topoFig = topographMD(MD, iFrame, iChan, topoMap, title0, figVisibleFlag)
% topographMD Visualize the spatial distribution of activities averaged
% over the whole time period on the specified frame of the cell image using
% plotWindows.m function.
%
% Usage:
%       topoFig = topographMD(MD, iFrame, iChan, topoMap, title0, figVisibleFlag)
%        
% Input:
%       MD          - a movieData object
%       iFrame      - an index of frame to plot
%       iChan       - a channel index
%       topoMap     - 2-dim'l array (windows, layers) of activities to plot
%       title0      - characters for title
%       figVisibleFlag - whether plotted or not
%
% Output:
%       topoFig     - a matlab figure object
%
% Jungsik Noh, 2016/06/21 


topoFig = figure('Visible', figVisibleFlag);
hold on

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
%iChan = 1;
windows = MD.processes_{iWinProc}.loadChannelOutput(iFrame, iChan);  %% win{iStrip}{iBand}
%preBandmax = 5;
%nBandMax_ = MD.processes_{iWinProc}.nBandMax_;

nColor = 1024;
winColors = colormap(jet(nColor));
%minval = min(topoMap(:));
%maxval = max(topoMap(:));

minval = quantile(topoMap(:), 0.01);
maxval = quantile(topoMap(:), 0.99);
topoMap(topoMap < minval) = minval;
topoMap(topoMap > maxval) = maxval;


%
    nStrip = numel(windows);
    for w = 1:nStrip
        nBandCur = numel(windows{w});
        for l = 1:nBandCur

            if ~isnan(topoMap(w,l))
                if (maxval > minval)
                    iCol = round(interp1(linspace(minval,maxval,1024),1:nColor,topoMap(w,l)));%This is probably slow - do it manually?                    
                    plotString = {winColors(iCol,:)};
                else
                    iCol = 1;
                    plotString = {winColors(iCol,:)};
                end
            else
                plotString = {'k','FaceAlpha',0};
            end
            if ~isempty(windows{w}{l})
                plotWindows(windows{w}{l},plotString, 50)
            end                
        end
    end
    
    axis off
    if (maxval > minval)
        caxis([minval, maxval])
    end
    colorbar
    title(title0)

refWin = 1:30:nStrip;
for k = 1:numel(refWin)
    if ~isempty(windows{refWin(k)}) 
    if ~isempty(windows{refWin(k)}{1}) 
        refwindows = windows{refWin(k)}{1};
        xx = refwindows{1}(1, 1);
        yy = refwindows{1}(2, 1);
        text(xx, yy, ['w', num2str(refWin(k))])
    end
    end
end
    
end
