function [maxSpeed,maxAcceleration,densityCell]=estimateTrackability(detections,aquisitionFreq,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addParameter('debugMode',[false]);
ip.parse(varargin{:});
p=ip.Results;

maxSpeed=cell(1,length(detections));
maxAcceleration=cell(1,length(detections));
densityCell=cell(1,length(detections));


for dIdx=1:length(detections)
    d=(detections(dIdx));
    
    M=[d.xCoord(:,1) d.yCoord(:,1) d.zCoord(:,1)]*0.1;
%     D = createSparseDistanceMatrix(M,M,Radius);
%     [sortedD,indx]=sort(full(D));
%     densityMinDist=1./((4/3)*pi*(sortedD(2,:)/2).^3);
%     densityOrder2MinDist=2./((4/3)*pi*(sortedD(2,:)+(sortedD(3,:)-sortedD(2,:))/2).^3);

    density=estimateDensity(d,10);
    densityCell{dIdx}=density{1};
    maxSpeed{dIdx}=aquisitionFreq*(3./(4*pi*density{1})).^(1/3);
    maxAcceleration{dIdx}=aquisitionFreq^2*(3/(4*pi*density{1})).^(1/3);

    if(strcmp(p.debugMode,'AmbiguityHeuristic'))
        dp1=detections(dIdx+1);
    	ambiguityHeuristicEstimation(d,dp1,6,0.5)
    end 
end
end

function [hVal,ambiguityCost]=ambiguityHeuristicEstimation(detectionAtT,detectionsAtTp1,searchRadius,percentile)
    % TODO: Complete, debug and ... probably revamped
	%H0 is trackable
	M=[detectionAtT.xCoord(:,1) detectionAtT.yCoord(:,1) detectionAtT.zCoord(:,1)]*0.1;
	N=[detectionsAtTp1.xCoord(:,1) detectionsAtTp1.yCoord(:,1) detectionsAtTp1.zCoord(:,1)]*0.1;
	D = createSparseDistanceMatrix(M,N,searchRadius);
	[sortedD,indx]=sort(full(D),2);
    closeIdx=arrayfun(@(r) find(sortedD(r,:)>0),1:size(sortedD,1),'unif',0);
    closeIdx=cellfun(@(c) c(1:min(2,end)),closeIdx,'unif',0);
    closeIdx=vertcat(closeIdx{:});
    
    %%
	closestZIdxForward=indx(sub2ind(size(indx),(1:size(indx,1)),closeIdx(:,1)'))';
	nextClosestZIdxForward=indx(sub2ind(size(indx),(1:size(indx,1)),closeIdx(:,2)'))';
	[bestIdxForWard, bestIdxBackward] = lap(D, [], [], 1);
	competingPredictions=bestIdxBackward(nextClosestZIdxForward);

	%% Best total cost
    unlinked=bestIdxForWard(1:size(D,2))>length(closestZIdxForward);
    linkedCost=bestIdxForWard(1:size(D,2));
    linkedCost=linkedCost(~unlinked);
    bestCost=zeros(size(closestZIdxForward));
    alternativeCost=zeros(size(closestZIdxForward));
    % debug and rather rederive with simple density for now 
	bestCost(unlinked)=abs(D(unlinked,linkedCost))+abs(D(unlinked,nextClosestZIdxForward));
    bestIsClosest=(bestIdxForWard==closestZIdxForward);
    alternativeCost(bestIsClosest)=abs(D(bestIsClosest,nextClosestZIdxForward))+abs(D(competingPredictions(bestIsClosest),closestZIdxForward));
    alternativeCost(~bestIsClosest)=abs(D(~bestIsClosest,bestIdxForWard))+abs(D(competingPredictions(~bestIsClosest),nextClosestZIdxForward));
	
	ambiguityCost=bestCost./alternativeCost;
	% Sum of possible association

	hVal=(ambiguityCost<percentile);
end

