function [ distToBranch] = GCAAnalysisExtract_distanceToBranch( analInfo,filoFilterSet )
%
%mkPlot = 1;
distToBranchCell = cell(length(analInfo)-1,1);

for iFrame = 1:length(analInfo) -1
    filoInfo = analInfo(iFrame).filoInfo;
    filterFrameC= filoFilterSet{iFrame};
    filoInfoFilt  = filoInfo(filterFrameC);
 
    conXYD =  vertcat(filoInfoFilt(:).conXYCoords); % this 
    % format gives the xy coords and the measured distance from veil to branch point
    distToBranch = conXYD(:,3) ; 
    distToBranchCell{iFrame} = distToBranchToVeil.*0.216;
    clear distToBranch
       
end

end





