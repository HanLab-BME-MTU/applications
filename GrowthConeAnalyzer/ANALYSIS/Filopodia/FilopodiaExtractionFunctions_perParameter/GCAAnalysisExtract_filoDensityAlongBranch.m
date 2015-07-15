function [ branchDensitiesCell] = GCAAnalysisExtract_filoDensityAlongBranch( analInfo,filoFilterSet )
%%% GCAAnalysisExtract_filoDensityAlongBranch
% Collects filopodia density along a branch for an entire movie for a
% filopodia filter set 
%mkPlot = 1;
branchDensitiesCell = cell(length(analInfo)-1,1);

for iFrame = 1:length(analInfo) -1
   
    filoInfo = analInfo(iFrame).filoInfo;
    filterFrameC= filoFilterSet{iFrame};
    filoInfoFilt  = filoInfo(filterFrameC); % after this filter should have only N order branch and its corresponding N-1 branchstem
 %% get the branch stems 
   NTypes = unique(vertcat(filoInfoFilt(:).type)); 
   NTypes(NTypes==0)= 1; % for now just switch the 0 order (no filo to 1) 
   typeStem = min(NTypes);
   
   if length(NTypes)==2 % filter ok
       % get the stem lengths 
       idxStem= vertcat(filoInfoFilt(:).type)==typeStem; 
       lengthsStem = vertcat(filoInfoFilt(idxStem).Ext_length); 
       
       filoInfoStem = filoInfoFilt(idxStem); 
       
       
       % get the number of branches per stem - tricky part here is this need to
       % likewise be filtered by fit etc which it will not be in the length of the .conIdx. 
       % get the number of filo 
       IDsCurrentSet = find(filterFrameC); 
       % for each stem get the conIdx and filter by filtInfo 
       numFiloBeforeFilt = arrayfun(@(x) length(filoInfoStem(x).conIdx),1:sum(idxStem)); 
       numFilo = arrayfun(@(x) length(intersect(filoInfoStem(x).conIdx,IDsCurrentSet)),1:sum(idxStem)) ;
       
       % test if that filo is in the current filter set 
       
       % density = number/length 
    
   % note maybe should include 0 order in this set?
    
    %% Get Filopodia Number 
    filterFrameC= filoFilterSet{iFrame};
    numFilo = sum(filterFrameC);
       
    %% Calculate Density 
    
  
    
    % for now just convert

    densitiesCell{iFrame,1} = numFilo/branchLength*10; % output 
   else 
       display('Check Branch Filter: N~=2'); 
end

%  if mkPlot ==1
%      scatter((1:length(densities))*5,densities,50,'k','filled');
%      ylabel('Filopodia Density OverTime','FontName','Arial','FontSize',14);
%      xlabel('Time (s)')
%      saveas(gcf,[saveDir filesep '001.fig']);
%
%  end
% if ~isempty(saveDir)
%
%  save([saveDir filesep 'filopodiaDensityCell'],'densitiesCell');
%  save([saveDir filesep 'toPlotMovie'],'densities');
end





