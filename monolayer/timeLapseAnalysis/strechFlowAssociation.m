function [rData, rSim, pVal,rsDelay] = strechFlowAssociation(cumsumStrech,cumsumFlow,maxTimeDelay)

[istart, iend] = findRelevantIndices(cumsumFlow);
cumsumStrech = cumsumStrech(istart:iend,:);
cumsumFlow = cumsumFlow(istart:iend,:);

[rData, pvalData] = (corr(cumsumStrech(:),cumsumFlow(:)));

rsDelay = getTemporalDelay(cumsumStrech,cumsumFlow,maxTimeDelay);

nIter = 5000;
rs = zeros(1,nIter);
pvals = zeros(1,nIter);
for i = 1 : nIter
    rndFlow = (cumsumFlow(randperm(size(cumsumFlow,1)),:));
    [rs(i), pvals(i)] = corr(cumsumStrech(:),rndFlow(:));
end

rSim = mean(rs);
pVal = sum(rData < rs) / nIter;
end

function [istart, iend] = findRelevantIndices(cumsumFlow)
istart = find(sum(cumsumFlow,2),1,'first');
iend = find(sum(cumsumFlow,2),1,'last');
end

function rsDelay = getTemporalDelay(cumsumStrech,cumsumFlow,maxTimeDelay)
rsDelay = nan(1,2*maxTimeDelay+1);
for lag = -maxTimeDelay : maxTimeDelay
    curStretch = cumsumStrech(:,max(1,1-lag):min(end,end-lag));
    curFlow = cumsumFlow(:,max(1,1+lag):min(end,end+lag));
    
    %     if lag > 0
    %         curStretch = cumsumStrech(:,1:end-lag);
    %         curFlow = cumsumFlow(:,);
    %     else if lag < 0
    %             curStretch = cumsumStrech
    %             curFlow = cumsumFlow
    %         else
    %             curStretch = cumsumStrech;
    %             curFlow = cumsumFlow;
    %         end
    %     end
    rsDelay(lag+maxTimeDelay+1) = corr(curStretch(:),curFlow(:));
end
end

% function [stretchFlow,protNoStretchFlow,noStretchFlow] = getFlowForStretchingAndNonStretchingEvents(flowAccumulationPatch,strechEvents,protrudingCellsNoStretch)
% stretchFlow = flowAccumulationPatch(strechEvents);
% protNoStretchFlow = flowAccumulationPatch(protrudingCellsNoStretch);
% noStretchFlow = flowAccumulationPatch(~strechEvents);
% end