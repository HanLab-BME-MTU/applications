function [ branchDensitiesCell] = GCAAnalysisExtract_filoDensityAlongBranch( analInfo,filoFilterSet )
%%% GCAAnalysisExtract_filoDensityAlongBranch
% Collects filopodia density along a branch for an entire movie for a
% filopodia filter set 
%mkPlot = 1;
branchDensitiesCell = cell(length(analInfo)-1,1);

for iFrame = 1:length(analInfo) -1
    
    filoInfo = analInfo(iFrame).filoInfo;
   
        filterFrameC= filoFilterSet{iFrame};
         if ~isempty(filterFrameC)
        filoInfoFilt  = filoInfo(filterFrameC(:,1)); % after this filter should have only N order branch and its corresponding N-1 branchstem
        %% get the branch stems
        types = vertcat(filoInfoFilt(:).type);
        % change type 0 (no filo) to type 1 as want to calculate density of
        % branching using ALL filo attached to veil.
        types(types==0) = 1;
        NTypes = unique(types);
        % NTypes(NTypes==0)= 1; % for now just switch the 0 order (no filo to 1)
        typeStem = min(NTypes);
        
        if length(NTypes)==2 % filter ok
            % get the stem lengths
            idxStem= vertcat(types)==typeStem;
            lengthStem = vertcat(filoInfoFilt(idxStem).Ext_length).*0.216;
            
            filoInfoStem = filoInfoFilt(idxStem);
            
            
            % get the number of branches per stem - tricky part here is this need to
            % likewise be filtered by fit etc which it will not be in the length of the .conIdx.
            % get the number of filo
            IDsCurrentSet = find(filterFrameC);
            % for each stem get the conIdx and filter by filtInfo
            numFiloBeforeFilt = arrayfun(@(x) length(filoInfoStem(x).conIdx),1:sum(idxStem));
            numFilo = arrayfun(@(x) length(intersect(filoInfoStem(x).conIdx,IDsCurrentSet)),1:sum(idxStem)) ;
            
            % test if that filo is in the current filter set
            
            density = numFilo'./lengthStem*10;
            
            % note maybe should include 0 order in this set?
            
            %% Calculate Density
            branchDensitiesCell{iFrame,1} = density; % output
        else
            if ~(isempty(filoInfoFilt) || length(filoInfoFilt)==1);
                display('Check Branch Filter: N~=2');
            end
        end
        
    else
        branchDensitiesCell{iFrame,1} = [];
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


