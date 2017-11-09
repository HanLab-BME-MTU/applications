function [ output_args ] = GCAGroupAnalysisPlotHMMTransitionsCart(toPlot,varargin)
% GCAGroupAnalysisPlotHMMTransitions

%% Check input 
ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('InvertTimeAxis',false,@(x) islogical(x)); % it turned out the data 
% is easier to interpret if start at 0 and move downward to 90s 
% However want to keep the same positive y as for the MDS scatter

ip.addParameter('plotNetVector',true); 

ip.addParameter('plotMag',true); 

ip.addParameter('groupIDs',[9,10]); 


ip.parse(varargin{:});
%%

if ip.Results.InvertTimeAxis 
    mult = -1; 
else 
    mult = 1; 
end 

groupIDs = ip.Results.groupIDs; 
nGroups = length(groupIDs); 
setAxis('on')
hold on

% filter by the max natural magnitude of the control 
% for now just hardcode the number of projects
% firstY = NaN(3,nGroups);
% firstU = NaN(3,nGroups);
% firstV = NaN(3,nGroups);

for iGroup = 1:nGroups
    nProjs = size(toPlot.info.projList{groupIDs(iGroup)},1); 
    cmap = toPlot.info.colorShades{groupIDs(iGroup)}; 
    
    % for now just hardwire for three projects as I know I only have three
    % projects for each perturbation 
    cmapFinal{iGroup}{1}  = cmap(1,:); %
    cmapFinal{iGroup}{2} = cmap(5,:); 
    cmapFinal{iGroup}{3} = cmap(end,:); 
    
   
    for iProj = 1:nProjs 
       x =  toPlot.HMM.transTime{groupIDs(iGroup)}{iProj};
       xconvert = (x.*5-5)-300; 
       if ~isempty(xconvert)
           idx = find(xconvert > 0 & xconvert < 90); 
           u = toPlot.HMM.transCart{groupIDs(iGroup)}{iProj}(:,1);
           v = toPlot.HMM.transCart{groupIDs(iGroup)}{iProj}(:,2); 
  scaleFactor = 10000; % for MDS 
 % scaleFactor = .3; % for tsne
%            firstU{iProj,iGroup} = u(idx);
%            firstV{iProj,iGroup} = v(idx); 
%            firstY{iProj,iGroup} = x(idx); 
firstU{iGroup}{iProj} = u(idx).*scaleFactor;
firstV{iGroup}{iProj} = v(idx).*scaleFactor;
%firstY{iGroup}{iProj} = (x(idx).*5 - 5) -300;
firstY{iGroup}{iProj} = xconvert(idx);
           %scatter(firstY,firstX,50,cmapFinal(iProj,:),'filled');
       end 
           
 
    end 
end
% firstU = firstU.*100;
% firstV = firstV.*100; 
% firstY = firstY-61; 
% firstY = firstY.*5-5; 
hold on 
scaleForAxis = 90/4; 
plotPos = scaleForAxis*([1,3]);

for iGroup = 1:nGroups;
    
    dataU = firstU{iGroup}; 
    dataV = firstV{iGroup}; 
    dataY = firstY{iGroup}; 
    cmap = cmapFinal{iGroup};
    % remove empty cells 
    idxKeep = cellfun(@(x) ~isempty(x),dataU);  
    
    dataU = dataU(idxKeep); 
    dataV = dataV(idxKeep); 
    dataY = dataY(idxKeep); 
    cmap = cmap(idxKeep); 
    
     %plot the line vectors for each group in different individual color 
     % shades 
     cellfun(@(u,v,y,c) arrayfun(@(i) line([plotPos(iGroup),plotPos(iGroup) + u(i)],[mult.*y(i),mult.*y(i) + v(i)],'Linewidth',2,...
         'color',c),1:length(u)),dataU,dataV,dataY,cmap); 
     netVect = cellfun(@(u,v) [sum(u) sum(v)],dataU,dataV,'uniformoutput',0);
     
     netVectAll{iGroup} = netVect;
     
      if ip.Results.plotNetVector
              % for now just plot in the middle the net vectors can move them later 
     

          
          cellfun(@(x,c) arrayfun(@(i) line([scaleForAxis*2,scaleForAxis*2+x(i,1)],[20*iGroup,20*iGroup+x(i,2)],'color',c ,...
              'Linewidth',2 ),1:size(x,1)),netVect,cmap);
          
          %      cellfun(@(x,c) arrayfun(@(i) scatter(scaleForAxis*2,20*iGroup,500,c,'filled','MarkerEdgeColor','w'),cmap);
          scatter(scaleForAxis*2,20*iGroup,500,toPlot.info.color{groupIDs(iGroup)},'filled','MarkerEdgeColor','w')
      end

      %netVect{iGroup} = cellfun(@(u,v) [sum(u) sum(v)],dataU,dataV,'uniformoutput',0); 
    
    % plot the scatter per group ; color each project by an individual
    % color shade. 
    cellfun(@(y,c) scatter(repmat(plotPos(iGroup),length(y),1),mult.*y,400,c,'filled','MarkerEdgeColor','w'),dataY,cmap); 
    
    projID{1,1} = '1'; 
    projID{1,2} = '2'; 
    projID{1,3} = '3'; 
    projID= projID(idxKeep); 
     cellfun(@(y,n) text(repmat(plotPos(iGroup)-0.75,length(y),1),mult.*y,n,'color','w','FontSize',20) , dataY,projID);
     
     
    
    
    
    
    % plot lines at x,y with
%     arrayfun(@(i) line([iGroup,iGroup+firstU(i,iGroup)],[firstY(i,iGroup),firstY(i,iGroup)+firstV(i,iGroup)],...
%         'Linewidth', 2, 'color',toPlot.info.color{groupIDs(iGroup)}),1:size(firstY,1));
% plot all the transitions after 


    %notBoxPlotQuiver(firstX,[],[],[],firstU,firstV);
%     scatter(repmat(iGroup,size(firstY,1),1),firstY(:,iGroup),100,toPlot.info.color{groupIDs(iGroup)},'filled','MarkerEdgeColor','w');
end

ylabel('Time After Perturbation (s)');
h1 = get(gcf,'CurrentAxes');


if ip.Results.InvertTimeAxis
    
    axis([0,90,-90,0]);

else
    axis([0,90,0,90]);
end

set(gca,'XTick',plotPos)
% yLim = h1.YLim; % use whatever they used
% axis([0.5 2.5 yLim(1) yLim(2)]);
%groupNames = toPlot.info.names(groupIDs);
groupNames{1} = '+ DMSO';
groupNames{2} = '+ 25 uM CK666';
% set(gca,'XTick',1:numel(groupNames));
set(gca,'XTickLabel',groupNames,'FontSize',20);

% saveas(gcf,'20161105_TransVects90s_FirstMDS.fig');
% saveas(gcf,'20161105_TransVects90s_FirstMDS.eps','psc2'); % old 

% saveas(gcf,'20161108_TransVects90s_FirstTsne.fig');
% saveas(gcf,'20161108_TransVects90s_FirstTsne.eps','psc2'); % new

saveas(gcf,[date 'TransVects90s_MDSNoVel.fig']); 
saveas(gcf,[date 'TransVects90s_MDSNoVel.eps'],'psc2'); 


close gcf
%      %% plot net vector for that time
%      setAxis('on')
%      hold on 
%      scatter(0,0,50,'k','filled')
%      for iGroup = 1:nGroups
%       [t,r] = cellfun(@(x) cart2pol(x(1),x(2)),netVectAll{iGroup});
%       d = rad2deg(t); 
%       
%       scatter(rad2deg(t),r/10000,50,toPlot.info.color{groupIDs(iGroup)},'filled'); 
%      end 
%%
if ip.Results.plotMag
    setAxis('on')
    netVectAll{1}{end+1} = [0,0]; % only for the MDS
    
    dataCell = cellfun(@(x) vertcat(x{:}),netVectAll,'uniformoutput',0);
    [t,r]= cellfun(@(x) cart2pol(x(:,1),x(:,2)),dataCell,'uniformoutput',0);
    % dataMagCell = cellfun(@(x) x(:,2),dataCell,'uniformoutput',0);
    %
    dataMat = reformatDataCell(r);
    h = notBoxPlot(dataMat./scaleFactor);
    
    for iGroup= 1: length(groupIDs)
        colorShades{iGroup}=  [toPlot.info.colorShades{groupIDs(iGroup)}(end,:) ; toPlot.info.colorShades{groupIDs(iGroup)}(7,:)];
    end
    %set the colors
    arrayfun(@(i) set(h(i).data,'markerFaceColor', colorShades{i}(1,:)),1:2);
    arrayfun(@(i) set(h(i).data,'markerEdgeColor','w'),1:2);
    arrayfun(@(i) set(h(i).mu,'color',colorShades{i}(1,:)),1:2);
    arrayfun(@(i) set(h(i).semPtch,'faceColor',colorShades{i}(2,:)),1:2);
    arrayfun(@(i) set(h(i).sdPtch,'faceColor',colorShades{i}(2,:)),1:2);
    
    
    ylabel({'Net Transition Magnitude :';'1.5 Min After Treatment'});
    %ylabel(toPlot.(params{iParam}).yLabel);
    set(gca,'XTick',1:numel(groupNames));
    set(gca,'XTickLabel',groupNames,'FontSize',20);
    h1 = get(gcf,'CurrentAxes');
    yLim = h1.YLim; % use whatever they used
    axis([0.5 2.5 yLim(1) yLim(2)]);
    
    % saveas(gcf,'20161105_NetTransMag90s_FirstMDS.fig');
    % saveas(gcf,'20161105_NetTransMag90s_FirstMDS.eps','psc2'); old
    
    % saveas(gcf,'20161108_NetTransMag90s_Firsttsne.fig');
    % saveas(gcf,'20171108_NetTransMag90s_Firsttsne.eps','psc2');
    
    saveas(gcf, [date 'NetTransMag90s_MDSNoVel.fig']);
    saveas(gcf, [date 'NetTransMag90s_MDSNoVel.eps'],'psc2');
end % if ip.Results.plotMag

end    
     
     
     
            
            


