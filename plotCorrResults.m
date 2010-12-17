function plotCorrResults(corr,maxDist,dFrames)
marker=['ro','bs','m*','c+','gd','yo','ks'];
[maxBin,numFrames]=size(corr);
ppr=2;        
for frame=1:dFrames:numFrames
    figure('Name',['cos(a)-correlation, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    title('This is a test')
    for binID=1:maxBin
        subplot(ceil(maxBin/ppr),ppr,binID)
        errorbarxy(corr(binID,frame).RMean4CosaMean,corr(binID,frame).cosaMean,corr(binID,frame).RMean4CosaSTD,corr(binID,frame).cosaSTD,[],[],marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)),marker(2*(mod(frame,7)+1)-1));
        title(['for cell to edge distance: ',num2str(corr(binID,frame).c2edMean,'%.1f'),'+-',num2str(corr(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([0,maxDist])
        ylim([-1,1])
    end
end

for frame=1:dFrames:numFrames
    figure('Name',['Full velocity correlation, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    title('This is a test')
    for binID=1:maxBin
        subplot(ceil(maxBin/ppr),ppr,binID)
        if ~isempty(corr(binID,frame).funcVel)
            plot(corr(binID,frame).RMean,corr(binID,frame).funcVel/corr(binID,frame).funcVel(1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        end
        title(['for cell to edge distance: ',num2str(corr(binID,frame).c2edMean,'%.1f'),'+-',num2str(corr(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([0,maxDist])
        ylim([-1,1])
    end
end

for frame=1:dFrames:numFrames
    figure('Name',['x-velocity correlation, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    title('This is a test')
    for binID=1:maxBin
        subplot(ceil(maxBin/ppr),ppr,binID)
        if ~isempty(corr(binID,frame).funcVelx)
            plot(corr(binID,frame).RMean,corr(binID,frame).funcVelx/corr(binID,frame).funcVelx(1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        end
        title(['for cell to edge distance: ',num2str(corr(binID,frame).c2edMean,'%.1f'),'+-',num2str(corr(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([0,maxDist])
        ylim([-1,1])
    end
end

for frame=1:dFrames:numFrames
    figure('Name',['y-velocity correlation, frame: ',num2str(frame,'%03.f')],'NumberTitle','off')
    title('This is a test')
    for binID=1:maxBin
        subplot(ceil(maxBin/ppr),ppr,binID)
        if ~isempty(corr(binID,frame).funcVely)
            plot(corr(binID,frame).RMean,corr(binID,frame).funcVely/corr(binID,frame).funcVely(1),marker(2*(mod(frame,7)+1)-1:2*(mod(frame,7)+1)));
        end
        title(['for cell to edge distance: ',num2str(corr(binID,frame).c2edMean,'%.1f'),'+-',num2str(corr(binID,frame).c2edSTD,'%.1f'),'[Pix]'])
        %xlim([0,maxR+binPixR])
        xlim([0,maxDist])
        ylim([-1,1])
    end
end