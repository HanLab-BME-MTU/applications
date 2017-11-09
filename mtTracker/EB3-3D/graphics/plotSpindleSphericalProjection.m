function plotSpindleSphericalProjection(handles,cumulAzi,cumulElev,cumulPoleId,cumulTimePt,cumulAziKin,cumulElevKin,cumulTimePtKin,cumulTrackIdKin,timeFrame,EB3MarkerSize,KinMarkerSize,cmapEB3,cmapKin,catchingEB3Id)

EB3InTimeFrame=(cumulTimePt<=timeFrame(2))&(cumulTimePt>=timeFrame(1));
KinInTimeFrame=(cumulTimePtKin==timeFrame(2));

plotTipIdx=(cumulElev<0)&(cumulPoleId==1)&EB3InTimeFrame;
plotSphericalProjection(handles(1),cumulAzi(plotTipIdx),abs(-pi/2-cumulElev(plotTipIdx)),ceil(cumulTimePt(plotTipIdx)),EB3MarkerSize,cmapEB3,'filled','Marker','o');
title(handles(1),[' astral Pole 1']);

% Pole 1

if(~isempty(cumulAziKin))
P1cumulElevKin=cumulElevKin(:,1);
P1cumulAziKin=cumulAziKin(:,1);
plotTipIdx=(P1cumulElevKin>0)&KinInTimeFrame;
plotSphericalProjection(handles(2),P1cumulAziKin(plotTipIdx,1),abs(pi/2-P1cumulElevKin(plotTipIdx,1)),cumulTrackIdKin(plotTipIdx),KinMarkerSize,cmapKin,'Marker','o'); hold(handles(2),'on');
end
plotTipIdx=(cumulElev>0)&(cumulPoleId==1)&EB3InTimeFrame;
plotSphericalProjection(handles(2),cumulAzi(plotTipIdx),abs(pi/2-cumulElev(plotTipIdx)),1,EB3MarkerSize,cmapEB3,'Marker','+'); hold(handles(2),'on')
if(~isempty(catchingEB3Id))
    CatchingEB3BufferTime=(cumulTimePt<=timeFrame(2))&(cumulTimePt>=max(timeFrame(2)-20,1));
    catchingEB3Mask=zeros(size(CatchingEB3BufferTime)); catchingEB3Mask(catchingEB3Id)=1;
    plotTipIdx=(cumulElev>0)&(cumulPoleId==1)&catchingEB3Mask&CatchingEB3BufferTime;
    localTime=ceil(cumulTimePt(plotTipIdx));
    m=min(localTime(:));M=max(localTime(:));
    colormapIdx=1;
    if(m<M)
        colormapIdx=max(min(ceil((localTime-m)*length(cmapEB3)/(M-m)),length(cmapEB3)),1);
    end
    plotSphericalProjection(handles(2),cumulAzi(plotTipIdx),abs(pi/2-cumulElev(plotTipIdx)),colormapIdx,EB3MarkerSize,cmapEB3,'filled','Marker','o'); hold(handles(2),'off')
end
hold(handles(2),'off')
title(handles(2),[' interpolar Pole 1']);

% Pole 2
if(~isempty(cumulAziKin))
P2cumulElevKin=cumulElevKin(:,2);
P2cumulAziKin=cumulAziKin(:,2);
plotTipIdx=(P2cumulElevKin>0)&KinInTimeFrame;
plotSphericalProjection(handles(3),P2cumulAziKin(plotTipIdx),abs(pi/2-P2cumulElevKin(plotTipIdx)),cumulTrackIdKin(plotTipIdx),KinMarkerSize,cmapKin,'Marker','o'); hold(handles(3),'on');
end
plotTipIdx=(cumulElev>0)&(cumulPoleId==2)&EB3InTimeFrame;
plotSphericalProjection(handles(3),cumulAzi(plotTipIdx),abs(pi/2-cumulElev(plotTipIdx)),1,EB3MarkerSize,cmapEB3,'Marker','+'); hold(handles(3),'on');
if(~isempty(catchingEB3Id))
    CatchingEB3BufferTime=(cumulTimePt<=timeFrame(2))&(cumulTimePt>=max(timeFrame(2)-20,1));
    catchingEB3Mask=zeros(size(CatchingEB3BufferTime)); catchingEB3Mask(catchingEB3Id)=1;
    plotTipIdx=(cumulElev>0)&(cumulPoleId==2)&catchingEB3Mask&CatchingEB3BufferTime;
    localTime=ceil(cumulTimePt(plotTipIdx));
    m=min(localTime(:));M=max(localTime(:));
    colormapIdx=1;
    if(m<M)
        colormapIdx=max(min(ceil((localTime-m)*length(cmapEB3)/(M-m)),length(cmapEB3)),1);
    end
    plotSphericalProjection(handles(3),cumulAzi(plotTipIdx),abs(pi/2-cumulElev(plotTipIdx)),colormapIdx,EB3MarkerSize,cmapEB3,'filled','Marker','o'); 
end
hold(handles(3),'off')
title(handles(3),[' interpolar Pole 2']);

plotTipIdx=(cumulElev<0)&(cumulPoleId==2)&EB3InTimeFrame;
plotSphericalProjection(handles(4),cumulAzi(plotTipIdx),abs(-pi/2-cumulElev(plotTipIdx)),ceil(cumulTimePt(plotTipIdx)),EB3MarkerSize,cmapEB3,'filled','Marker','o')
title(handles(4),[' astral Pole 2']);



