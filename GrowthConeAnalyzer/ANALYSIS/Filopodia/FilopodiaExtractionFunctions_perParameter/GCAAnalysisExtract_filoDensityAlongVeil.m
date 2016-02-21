function [ densitiesCell] = GCAAnalysisExtract_filoDensityAlongVeil( analInfo,filoFilterSet, veilStem )
%
%mkPlot = 1;
densitiesCell = cell(length(analInfo)-1,1);

for iFrame = 1:length(analInfo) -1
    % get filopodia 
   
    
 
    
    % get masks 
    neuriteMask = veilStem(iFrame).finalMask;
    
    % sort pixels
    [ny,nx] = size(neuriteMask);
    roiYX = bwboundaries(neuriteMask);
    edgeMask = zeros([ny,nx]);
    
    pixEdge =  sub2ind([ny,nx],roiYX{1}(:,1),roiYX{1}(:,2));
    edgeMask(pixEdge) = 1;
    
    % add a thinning step
    edgeMask = bwmorph(edgeMask,'thin','inf');
    
    
    % take out border pixels
    
    boundaryMask = zeros([ny,nx]);
    boundaryMask(1:ny,1) =1;
    boundaryMask(1:ny,nx)=1;
    boundaryMask(1,1:nx)= 1;
    boundaryMask(ny,1:nx) =1;
    
    edgeMask(boundaryMask==1) = 0;
    pixEdgeMask = find(edgeMask==1);
    EPs = getEndpoints(pixEdgeMask,[ny,nx],0);
    
    
    pixIdxBack = nan(length(pixEdgeMask),1); % overinitialize to make happy
    edgeMask = logical(edgeMask);
    
    if isempty(EPs)
        % try breaking
        EPs = getEndpoints(pixEdgeMask(2:end),[ny,nx],1);
    end
    % before 20141221 transform thin as will ultimately have two pixels
    % with
    
    transform = bwdistgeodesic(edgeMask,EPs(1,1),EPs(1,2));
    
    iPix = 0;
    while length(find(transform==iPix)) == 1
        pixIdxBack(iPix+1) = find(transform==iPix); % start at the endpoint
        iPix = iPix +1;
    end
    
    distBoundMicron =  calculateDistance(pixIdxBack,[ny,nx],0) ; % note currently calculate distance is
    
    
    %% Get Filopodia Number 
    filterFrameC= filoFilterSet{iFrame};
    if ~isempty(filterFrameC)
    numFilo = sum(filterFrameC(:,1));
    
    
    
    
    %% Calculate Density 
    
  
    
    % for now just convert
    
    
    
    
    densitiesCell{iFrame,1} = numFilo/distBoundMicron*10;
    else 
        densitiesCell{iFrame,1} = []; 
    
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





