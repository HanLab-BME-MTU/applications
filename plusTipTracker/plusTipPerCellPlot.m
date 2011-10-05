function [ cMap] = perCellPlot(groupList,dataStruct,mean_stdOverAllCellsInGroup, fieldName,cMap)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


projGroupName=groupList(:,1);
projGroupDir=cellfun(@(x) formatPath(x),groupList(:,2),'uniformoutput',0);


% fix the names if there are spaces or hyphens and append prefix 'grp'
projGroupName=cellfun(@(x) strrep(x,'-','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) strrep(x,' ','_'),projGroupName,'uniformoutput',0);
projGroupName=cellfun(@(x) ['grp_' x],projGroupName,'uniformoutput',0);


% count unique groups and keep them in order of the original list
[btwGrpNames,m,projGroupIdx] = unique(projGroupName);
[b,idx]=sort(m);
btwGrpNames=btwGrpNames(idx);


nGroups = length(dataStruct);

 
for iGroup = 1:nGroups
 if iGroup==1
            % properties for plotting
            nRep=ceil(nGroups/3);
            mrkTpe=repmat(['o';'^';'s'],nRep,1); % we'll make this into a cell later
            prop_name(1) = {'Marker'};
            prop_name(2) = {'MarkerFaceColor'};
            prop_name(3) = {'MarkerEdgeColor'};

            % set the colormap to be random
            if isempty(cMap)
                cM=varycolor(nGroups);
                cM=cM(randsample(nGroups,nGroups),:);
                cMap=mat2cell(cM,ones(nGroups,1),3);
            end

            figure; hold on
 end
        
 nC=length(dataStruct(iGroup).(fieldName));
        prop_values(1:nC,1) = {mrkTpe(iGroup)};
        prop_values(1:nC,2) = cMap(iGroup,:);
        prop_values(1:nC,3) = cMap(iGroup,:);
        
   % use plot instead of scatter so more flexibility with properties. to
        % do this, make 2 x nPoints matrix where the second row is all NaNs and
        % then use the plot function
        if iGroup == 1
             xStart = 1; 
             xEnd = 0; 
        else 
            xStart = xCoord(1,end) + 1; 
            xEnd = xCoord(1,end);   
        end 
        
        x1 = (xStart:length(dataStruct(iGroup).(fieldName))+xEnd);
        x2 = nan(size(x1)); 
        
        xCoord = [x1 ; x2]; 
        xValues(iGroup) = mean(xCoord(1,:));
        y1 = abs(dataStruct(iGroup).(fieldName));
        y2 = nan(size(y1));
        yCoord= [y1 ; y2]; 
        
     %subplot(2,1,1) 
     %hold on
      h =  plot(xCoord,yCoord,'.'); 
      
       set(h,prop_name,prop_values,'MarkerSize',16)
        % store the handle for the first member of each group for the legend
       h1(iGroup)=h(1);

        % these have to be defined each time
        clear prop_values
        
        %FOR boxplots - this has a bug in it now 
      % if length(y1) < 30 
       %   yfix = [y1 ; nan(30-length(y1),1)']; 
       %else 
       %end 
       %if iGroup == 1
        %   dataStructAll = yfix; 
        %else 
         %   dataStructAll = [dataStructAll yfix]; 
      % end 
        
        
       
        % keep all the parameter values in a vector
        %gsAll(c:c+length(paramValues)-1,1)=paramValues;

        %c=c+length(paramValues);
end       
       
      hold on 
      
      for iGroup = 1:length(btwGrpNames) 
         
          yValue = cell2mat(mean_stdOverAllCellsInGroup.(fieldName)(1+iGroup,2)); 
          
          
          plot(xValues(iGroup), yValue,'Marker','+','MarkerSize',16,'MarkerEdgeColor',cMap{iGroup},'MarkerFaceColor',cMap{iGroup});
      end
      
      % this was for me
% btwGrpNames=cellfun(@(x)strrep(x,'grp_6_22_10','dsRNA-'),btwGrpNames,'uniformoutput',0); 
% btwGrpNames=cellfun(@(x) strrep(x,'Orbit','CLASP'),btwGrpNames,'uniformoutput',0);
   legend(h1,btwGrpNames,'location','BestOutside');
   
  figure
   
   % subplot(2,1,2)
    %boxplot(dataStructAll,btwGrpNames,'labelorientation','inline','widths',0.3);
    %boxplot(dataStructAll)
    
    
end

