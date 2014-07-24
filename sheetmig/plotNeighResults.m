function plotNeighResults(neighborhood,maxDist,dFrames)
marker=['ro','bs','m*','c+','gd','yo','ks'];
[maxBin,numFrames]=size(neighborhood);
ppr=2;        
for frame=1:dFrames:numFrames
    figure('Name',['Nearest Neighbor Distance, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    for binID=1:maxBin
        if~isempty(neighborhood(binID,frame).c2cDis)
        subplot(ceil(maxBin/ppr),ppr,binID)
        errorbar(1:size(neighborhood(binID,frame).c2cDis,2),nanmean(neighborhood(binID,frame).c2cDis,1),nanstd(neighborhood(binID,frame).c2cDis,[],1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        title(['for cell to edge distance: ',num2str(neighborhood(binID,frame).c2edMean,'%.1f'),'+-',num2str(neighborhood(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([1,60])
        ylim([1,200])
        end
    end
end

marker=['ro','bs','m*','c+','gd','yo','ks'];
[maxBin,numFrames]=size(neighborhood);
ppr=2;        
for frame=1:dFrames:numFrames
    figure('Name',['Change of Nearest Neighbor Distance , frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    for binID=1:maxBin
        if~isempty(neighborhood(binID,frame).c2cDD)
        subplot(ceil(maxBin/ppr),ppr,binID)
        errorbar(1:size(neighborhood(binID,frame).c2cDis,2),nanmean(neighborhood(binID,frame).c2cDD,1),nanstd(neighborhood(binID,frame).c2cDD,[],1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        title(['for cell to edge distance: ',num2str(neighborhood(binID,frame).c2edMean,'%.1f'),'+-',num2str(neighborhood(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([1,60])
        ylim([1,200])
        end
    end
end

marker=['ro','bs','m*','c+','gd','yo','ks'];
[maxBin,numFrames]=size(neighborhood);
ppr=2;        
for frame=1:dFrames:numFrames
    figure('Name',['Nearest Neighbor Distance Normalized by traveled distance, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    for binID=1:maxBin
        if~isempty(neighborhood(binID,frame).c2cDD)
        subplot(ceil(maxBin/ppr),ppr,binID)
        %display('What to do with 0/0???')
        errorbar(1:size(neighborhood(binID,frame).c2cDis,2),nanmean(neighborhood(binID,frame).c2cDD./neighborhood(binID,frame).trvDis,1),nanstd(neighborhood(binID,frame).c2cDD./neighborhood(binID,frame).trvDis,[],1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        title(['for cell to edge distance: ',num2str(neighborhood(binID,frame).c2edMean,'%.1f'),'+-',num2str(neighborhood(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([1,60])
        ylim([0,1])
        end
    end
end

marker=['ro','bs','m*','c+','gd','yo','ks'];
[maxBin,numFrames]=size(neighborhood);
ppr=2;        
for frame=1:dFrames:numFrames
    figure('Name',['Nearest Neighbor Distance Normalized by direct distance, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    for binID=1:maxBin
        if~isempty(neighborhood(binID,frame).c2cDD)
        subplot(ceil(maxBin/ppr),ppr,binID)
        %display('What to do with 0/0???')
        errorbar(1:size(neighborhood(binID,frame).c2cDis,2),nanmean(neighborhood(binID,frame).c2cDD./neighborhood(binID,frame).dirDis,1),nanstd(neighborhood(binID,frame).c2cDD./neighborhood(binID,frame).dirDis,[],1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        title(['for cell to edge distance: ',num2str(neighborhood(binID,frame).c2edMean,'%.1f'),'+-',num2str(neighborhood(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([1,60])
        ylim([0,1])
        end
    end
end