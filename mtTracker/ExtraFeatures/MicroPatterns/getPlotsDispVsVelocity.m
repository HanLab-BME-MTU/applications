function [ forPlot ] = getPlotsDispVsVelocity(saveDir,groupList)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% takes the 
if nargin<1 || isempty(saveDir)
    saveDir=uigetdir(pwd,'Select output directory for plots.');
end


projGroupName=groupList(:,1);
projGroupDir= groupList(:,2);

% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);

groupNames = cell(length(btwGrpNames),1);
for iGroup = 1:length(btwGrpNames)
    groupNames(iGroup,1)= btwGrpNames(iGroup,1);
end 


for iGroup = 1:length(btwGrpNames)
  % get indices of projects in iGroup
    tempIdx= strmatch(btwGrpNames(iGroup),projGroupName,'exact');
    
    

   for iProj = 1:length(tempIdx)
     
       temp = load([projGroupDir{tempIdx(iProj)} filesep 'meta' filesep 'projData']);
    
       pathup1 = getFilenameBody(projGroupDir{tempIdx(iProj)});
       pathup2 = getFilenameBody(pathup1);
       
       roiMask = imread([pathup2 filesep 'roiMaskSeg.tif']);
       [imL,imW] = size(roiMask);
       dataMat = temp.projData.dataMat_FOR_STATS;
       xMat = temp.projData.xMatSubtracks;
       yMat = temp.projData.yMatSubtracks;
      
      
      
      avgX = zeros(length(xMat(:,1)),1);
      avgY = zeros(length(yMat(:,1)),1);
      for iSubtrack = 1:length(xMat(:,1))
          xCoords = xMat(iSubtrack,:);
          avgX(iSubtrack,1) = mean(xCoords(~isnan(xCoords)));
          yCoords = yMat(iSubtrack,:);
          avgY(iSubtrack,1) = mean(yCoords(~isnan(yCoords)));
      end
      
      
      pixIdx = sub2ind([imL,imW],ceil(avgY-.5),ceil(avgX-.5));
      
      weightedRoi=bwdist(swapMaskValues(roiMask)).*temp.projData.pixSizeNm/1000;
      
      if iProj == 1
          
      forPlot(iGroup).distValues = weightedRoi(pixIdx);
      forPlot(iGroup).vel = dataMat(dataMat(:,5) == 1,4); 
      
      else 
        forPlot(iGroup).distValues = [forPlot(iGroup).distValues ; weightedRoi(pixIdx) ];
        forPlot(iGroup).vel = [forPlot(iGroup).vel; dataMat(dataMat(:,5) == 1,4)];
      end 
        
   end % for iProj
   
end % for iGroup
           
for iGroup = 1:length(btwGrpNames)
    figure;
    scatter(forPlot(iGroup).distValues,forPlot(iGroup).vel,7,'filled');
    forTitle = strrep(char(btwGrpNames(iGroup)),'_',' ');
    title({forTitle ; 'Distance From the Cell Edge (um) Vs. Growth Velocity (um/min)'},'FontWeight','bold','FontSize',14);
    xlabel('Distance From Cell Edge (um)','Fontweight','bold','FontSize',12);
    ylabel('Growth Velocity of Subtrack (um/min)', 'Fontweight','bold','FontSize',12);
    saveas(gcf,[saveDir filesep 'Dis_Vs_Vel' char(btwGrpNames(iGroup)) '.fig']);
    saveas(gcf,[saveDir filesep 'Dis_Vs_Vel' char(btwGrpNames(iGroup)) '.eps']);
    save([saveDir filesep 'DisVsVel_ForPlot'],'forPlot');
end     
end

