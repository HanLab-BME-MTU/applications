function h = plotIntensityForceSingle(curTrack,curID)
    h=figure;
    set(h,'Position',[200,400,400,200])
    
    ax8=subplot(1,2,1);
    curStartFrameEE=curTrack.startingFrameExtraExtra;
    curEndFrameEE=curTrack.endingFrameExtraExtra;
    chosenStartFrame = curTrack.startingFrameExtra;
    chosenEndFrame = curTrack.endingFrameExtra;
    plot((curStartFrameEE:curEndFrameEE),curTrack.ampTotal(curStartFrameEE:curEndFrameEE),'-k'), hold on
    plot((curStartFrameEE:curEndFrameEE),curTrack.bkgAmp(curStartFrameEE:curEndFrameEE),'Color',[0.5 0.5 0.5])
    plot((chosenStartFrame:chosenEndFrame),curTrack.ampTotal(chosenStartFrame:chosenEndFrame),'-b','Marker','o','MarkerFaceColor','blue')
    xlabel('Time (s)'); ylabel('Fluorescence intensity (a.u.)')
    set(ax8,'FontSize',8)
    % force time series
    ax9=subplot(1,2,2);
    plot((curStartFrameEE:curEndFrameEE),curTrack.forceMag(curStartFrameEE:curEndFrameEE),'-k'), hold on
    plot((chosenStartFrame:chosenEndFrame),curTrack.forceMag(chosenStartFrame:chosenEndFrame),'-r','Marker','o','MarkerFaceColor','red')
    xlabel('Time (s)'); ylabel('Traction (Pa)')
    title(num2str(curID))
    set(ax9,'FontSize',8)
    
    splineParam = 0.1;
    d = curTrack.ampTotal(curStartFrameEE:curEndFrameEE);
    tRange = curTrack.iFrame(curStartFrameEE:curEndFrameEE);
%     numNan = find(isnan(d),1,'last');
%     tRange(isnan(d)) = [];
%     d(isnan(d)) = [];
    axes(ax8)
    sd_spline= csaps(tRange,d,splineParam);
    sd=ppval(sd_spline,tRange);
    plot((curStartFrameEE:curEndFrameEE),sd,'Color',[229/255 84/255 0],'Linewidth',2)

    axes(ax9)
    curForce = curTrack.forceMag(curStartFrameEE:curEndFrameEE);
    sCurForce_spline= csaps(tRange,curForce,splineParam);
    sCurForce=ppval(sCurForce_spline,tRange);
    plot((curStartFrameEE:curEndFrameEE),sCurForce,'Color',[229/255 84/255 0],'Linewidth',2)
end