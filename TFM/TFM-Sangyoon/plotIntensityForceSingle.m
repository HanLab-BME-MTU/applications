function h = plotIntensityForceSingle(curTrack,curID,intenOrForce,splineParamPeak,splineParamInit,tInterval)
    if nargin<6
        tInterval = 1;
    end
    if nargin<5
        splineParamInit = 0.5;
    end
    if nargin<4
        splineParamPeak = 0.1;
        splineParamInit = 0.5;
    end
    if nargin<3
        intenOrForce = 0;
        h=figure;
        set(h,'Position',[200,400,400,200])
    end
    curStartFrameEE=curTrack.startingFrameExtraExtra;
    curEndFrameEE=curTrack.endingFrameExtraExtra;
    chosenStartFrame = curTrack.startingFrameExtra;
    chosenEndFrame = curTrack.endingFrameExtra;
    if intenOrForce==0
        subplot(1,2,1);
    end
    
    d = curTrack.ampTotal(curStartFrameEE:curEndFrameEE);
    tRange = curTrack.iFrame(curStartFrameEE:curEndFrameEE);
    if intenOrForce==0 || intenOrForce==1
%         if intenOrForce==0
        plot((curStartFrameEE:curEndFrameEE)*tInterval,curTrack.ampTotal(curStartFrameEE:curEndFrameEE),'-k'), hold on
        plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.ampTotal(chosenStartFrame:chosenEndFrame),'-b')
%         else
%             plot((chosenStartFrame:chosenEndFrame),curTrack.ampTotal(chosenStartFrame:chosenEndFrame),'-k'), hold on
%             tRange=chosenStartFrame:chosenEndFrame;
%         end
        plot((curStartFrameEE:curEndFrameEE)*tInterval,curTrack.bkgAmp(curStartFrameEE:curEndFrameEE),'Color',[0.5 0.5 0.5])
        xlabel('Time (s)'); ylabel('Fluorescence intensity (a.u.)')
        set(gca,'FontSize',7)
        sd_splinePeak= csaps(curStartFrameEE:curEndFrameEE,d,splineParamPeak);
        sd=ppval(sd_splinePeak,tRange);
        plot(tRange*tInterval,sd,'Color',[229/255 84/255 0],'Linewidth',2)
        sd_splineInit= csaps(curStartFrameEE:curEndFrameEE,d,splineParamInit);
        sdinit=ppval(sd_splineInit,tRange);
        plot(tRange*tInterval,sdinit,'Color',[58/255 84/255 212/255],'Linewidth',1)
        %background level
        if ~isempty(curTrack.bkgMaxInt)
            line([0 (curEndFrameEE)*tInterval],[curTrack.bkgMaxInt curTrack.bkgMaxInt],'linestyle',':','Color','k');
        end
        if ~isempty(curTrack.forceTransmitting) && curTrack.forceTransmitting
            plot((curTrack.firstIncreseTimeInt/tInterval)*tInterval,...
                sdinit(round(curTrack.firstIncreseTimeInt/tInterval-curStartFrameEE+1)),'o','MarkerFaceColor','b','MarkerEdgeColor','w')
        end
        if ~isempty(curTrack.intenPeakness) && curTrack.intenPeakness
            plot((curTrack.intenPeakFrame)*tInterval,sd(curTrack.intenPeakFrame-curStartFrameEE+1),'o',...
                'MarkerFaceColor','w','MarkerEdgeColor',[84/255 84/255 255/255])
        end
    end
    xlim([tRange(1)*tInterval tRange(end)*tInterval])
    
    if intenOrForce==0
        subplot(1,2,2);
    end
    if intenOrForce==0 || intenOrForce==2
        % force time series
        plot((curStartFrameEE:curEndFrameEE)*tInterval,curTrack.forceMag(curStartFrameEE:curEndFrameEE),'-k'), hold on
        if intenOrForce==0
            plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.forceMag(chosenStartFrame:chosenEndFrame),'-r')
            title(num2str(curID))
        else
            plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.forceMag(chosenStartFrame:chosenEndFrame),'-b')
        end
%         else
%             plot((chosenStartFrame:chosenEndFrame),curTrack.forceMag(chosenStartFrame:chosenEndFrame),'-k'), hold on
%             tRange=chosenStartFrame:chosenEndFrame;
%         end
        xlabel('Time (s)'); ylabel('Traction (Pa)')
        set(gca,'FontSize',7)

    %     numNan = find(isnan(d),1,'last');
    %     tRange(isnan(d)) = [];
    %     d(isnan(d)) = [];
%         axes(ax9)
        curForce = curTrack.forceMag(curStartFrameEE:curEndFrameEE);
        sCurForce_spline= csaps(curStartFrameEE:curEndFrameEE,curForce,splineParamPeak);
        sCurForce=ppval(sCurForce_spline,tRange);
        plot(tRange*tInterval,sCurForce,'Color',[229/255 84/255 0],'Linewidth',2)
        sCurForce_splineInit= csaps(curStartFrameEE:curEndFrameEE,curForce,splineParamInit);
        sCurForceInit=ppval(sCurForce_splineInit,tRange);
        plot(tRange*tInterval,sCurForceInit,'Color',[58/255 84/255 212/255],'Linewidth',1)
        if ~isempty(curTrack.forceTransmitting) && curTrack.forceTransmitting
            plot(curTrack.firstIncreseTimeForce,...
                sCurForceInit(round(curTrack.firstIncreseTimeForce/tInterval)-curStartFrameEE+1),'o',...
                'MarkerFaceColor','b','MarkerEdgeColor','w')
        end
        if ~isempty(curTrack.forcePeakness) && curTrack.forcePeakness
            plot((curTrack.forcePeakFrame)*tInterval,sCurForce(curTrack.forcePeakFrame-curStartFrameEE+1),'o',...
                'MarkerFaceColor','w','MarkerEdgeColor',[229/255 84/255 84/255])
        end
        %background level
        if ~isempty(curTrack.bkgMaxForce)
            line([0 (curEndFrameEE)*tInterval],[curTrack.bkgMaxForce curTrack.bkgMaxForce],'linestyle',':','Color','k')
        end
    end
    if intenOrForce==3
        % edge advance time series
%         plot((curStartFrameEE:curEndFrameEE),curTrack.edgeAdvanceDist(curStartFrameEE:curEndFrameEE),'-k'), hold on
        plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.edgeAdvanceDist(chosenStartFrame:chosenEndFrame),'-k'), hold on
        xlabel('Time (s)'); ylabel('Edge protrusion (pixel)')
        set(gca,'FontSize',7)
        curEdge = curTrack.edgeAdvanceDist(chosenStartFrame:chosenEndFrame);
        tRange=chosenStartFrame:chosenEndFrame;
        sCurEdge_spline= csaps(tRange,curEdge,splineParamPeak);
        sCurEdge=ppval(sCurEdge_spline,tRange);
        plot((tRange)*tInterval,sCurEdge,'Color',[229/255 84/255 0],'Linewidth',2)
    end
    if intenOrForce==4
        % edge advance time series
%         plot((curStartFrameEE:curEndFrameEE),curTrack.edgeAdvanceDist(curStartFrameEE:curEndFrameEE),'-k'), hold on
        plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.advanceDist(chosenStartFrame:chosenEndFrame),'-k'), hold on
        xlabel('Time (s)'); ylabel('Adhesion advance (pixel)')
        set(gca,'FontSize',7)
        curEdge = curTrack.advanceDist(chosenStartFrame:chosenEndFrame);
        tRange=chosenStartFrame:chosenEndFrame;
        sCurEdge_spline= csaps(tRange,curEdge,splineParamPeak);
        sCurEdge=ppval(sCurEdge_spline,tRange);
        plot((tRange)*tInterval,sCurEdge,'Color',[229/255 84/255 0],'Linewidth',2)
    end
    if intenOrForce==5
        % edge advance time series
%         plot((curStartFrameEE:curEndFrameEE),curTrack.edgeAdvanceDist(curStartFrameEE:curEndFrameEE),'-k'), hold on
        plot((chosenStartFrame:chosenEndFrame)*tInterval,curTrack.distToEdge(chosenStartFrame:chosenEndFrame),'-k'), hold on
        xlabel('Time (s)'); ylabel('Distance to edge (pixel)')
        set(gca,'FontSize',7)
        curEdge = curTrack.distToEdge(chosenStartFrame:chosenEndFrame);
        tRange=chosenStartFrame:chosenEndFrame;
        sCurEdge_spline= csaps(tRange,curEdge,splineParamPeak);
        sCurEdge=ppval(sCurEdge_spline,tRange);
        plot((tRange)*tInterval,sCurEdge,'Color',[229/255 84/255 0],'Linewidth',2)
    end
    xlim([tRange(1)*tInterval tRange(end)*tInterval])
end