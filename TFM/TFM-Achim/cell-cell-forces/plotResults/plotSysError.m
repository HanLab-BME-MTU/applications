function []=plotSysError(groupedClusters)

goodCellSet=findCells(groupedClusters);
goodCluster=unique(horzcat(goodCellSet.clusterId));
   
tooShort=0;
for iCluster=goodCluster
    close all
    currErrVecList=[];
    for iFrame=1:length(groupedClusters.cluster{iCluster}.trackedNet)
        if ~isempty(groupedClusters.cluster{iCluster}.trackedNet{iFrame})
            currErrVecList=vertcat(currErrVecList,groupedClusters.cluster{iCluster}.trackedNet{iFrame}.stats.errorSumForce.vec);
        end
    end
    figure(1)
    plot(currErrVecList(:,1),currErrVecList(:,2),'*')
    hold on;
    plot([-1.1*max(abs(currErrVecList(:,1))) 1.1*max(abs(currErrVecList(:,1)))], [0 0],'k')
    plot([0 0], [-1.1*max(abs(currErrVecList(:,2))) 1.1*max(abs(currErrVecList(:,2)))],'k')
    xlim([-1.1*max(abs(currErrVecList(:,1))) 1.1*max(abs(currErrVecList(:,1)))])
    ylim([-1.1*max(abs(currErrVecList(:,2))) 1.1*max(abs(currErrVecList(:,2)))])
    axis equal
    hold off;
    
   
    %corrected:
    if length(currErrVecList(:,1))>10
        currErrVecListNew=[currErrVecList(:,1)-nanmean(currErrVecList(:,1)) currErrVecList(:,2)-nanmean(currErrVecList(:,2))];
    else
        currErrVecListNew=currErrVecList;
        tooShort=tooShort+1;
    end
    
    figure(2)
    plot(currErrVecListNew(:,1),currErrVecListNew(:,2),'*')
    hold on;
    plot([-1.1*max(abs(currErrVecListNew(:,1))) 1.1*max(abs(currErrVecListNew(:,1)))], [0 0],'k')
    plot([0 0], [-1.1*max(abs(currErrVecListNew(:,2))) 1.1*max(abs(currErrVecListNew(:,2)))],'k')
    xlim([-1.1*max(abs(currErrVecListNew(:,1))) 1.1*max(abs(currErrVecListNew(:,1)))])
    ylim([-1.1*max(abs(currErrVecListNew(:,2))) 1.1*max(abs(currErrVecListNew(:,2)))])
    axis equal
    hold off
    
    avgErr    = nanmean(sqrt(nansum(currErrVecList.^2   ,2)),1);
    avgErrNew = nanmean(sqrt(nansum(currErrVecListNew.^2,2)),1);
    relImprov = avgErrNew/avgErr;
    
    errAnalysis(iCluster).errVecList    = currErrVecList;
    errAnalysis(iCluster).errVecListNew = currErrVecListNew;
    errAnalysis(iCluster).avgErr        = avgErr;
    errAnalysis(iCluster).avgErrNew     = avgErrNew;
    errAnalysis(iCluster).relImprov     = relImprov;
    if relImprov<0.95
        display(['relImprov: ',num2str(relImprov),' for cluster: ',num2str(iCluster),' fileName: ',groupedClusters.clusterList{iCluster}]);
        %input('')
    end
end
hist(horzcat(errAnalysis.relImprov),10)
mean(horzcat(errAnalysis.relImprov))
tooShort