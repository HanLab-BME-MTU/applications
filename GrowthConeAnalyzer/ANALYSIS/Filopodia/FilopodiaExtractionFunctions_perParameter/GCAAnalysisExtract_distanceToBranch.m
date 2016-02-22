function [ distToBranchCell] = GCAAnalysisExtract_distanceToBranch( analInfo,filoFilterSet )
%
%mkPlot = 1;
distToBranchCell = cell(length(analInfo)-1,1);

% for iFrame = 1:length(analInfo) -1
for iFrame = 1:length(analInfo)-1
    
    filoInfo = analInfo(iFrame).filoInfo;
    filterFrameC= filoFilterSet{iFrame}(:,1); % take first column don't consider embedded
    filoInfoFilt  = filoInfo(filterFrameC);
    
    if ~isempty(filoInfoFilt)
        types = vertcat(filoInfoFilt(:).type);
        
        
        %take out type 0 - veil/stem filo with no branch (just in case included
        %in the filter set.
        if ~isempty(types ==0)
            filoInfoFilt(types==0) = [];
        end
        
        typesNoZ = vertcat(filoInfoFilt(:).type);
        NTypes = unique(typesNoZ);
        
        % get the stem filo - this will be the lowest order for the set
        % ie if filtered type 1 and type 2 , 1 will be the stem and 2 the
        % branch, if 2 and 3, 2 will be the stem, and 3 the branch ..etc
        typeStem = min(NTypes);
        idxStem= vertcat(typesNoZ)==typeStem;
        
        filoInfoStem = filoInfoFilt(idxStem);
        
        % do not extract distance measurements for branches not included in the filter set
        % (this removes for example small filo or filo that didn't pass the
        % filtering criteria and give the user some flexibility to this setting post processing)
        
        % find all the IDs of the current filo from the filter (each
        % filo/branch of the original set will have it's own ID)
        IDsCurrentSet = find(filterFrameC);
        
        % get individual filters for each list of filo branches connected to
        % the stem based on the the presence of the filo ID in the filter set
        % (Note another way of dealing with this would have simply been to
        % save the information as to the coordinate of connection with the
        % individual branch- I chose not to do this so as not to have redundant
        % information. The conIdx always indicates the ID of any potential
        % branches connected directly to that filo and the position and distance of each attachment is
        % saved in conXY.
        [~,idxKeep] = arrayfun(@(x) intersect(filoInfoStem(x).conIdx,IDsCurrentSet),1:sum(idxStem),'uniformoutput',0);
        
        %     % remove empty
            noIntersect = cellfun(@(x) isempty(x),idxKeep);
            filoInfoStem(noIntersect) = [];
            idxKeep(noIntersect) = [];
        %
        % get rid of empty cells 
        %idxKeep = idxKeep(cellfun(@(x) ~isempty(x),idxKeep)); 
        %if ~isempty(idxKeep)
            
            distValues = arrayfun(@(x) filoInfoStem(x).conXYCoords(idxKeep{x},3),1:length(idxKeep),'uniformoutput',0); % note if the length is zero 
            % it will simply result in an empty double which is what you
            % would like 
            IDs  = arrayfun(@(x) filoInfoStem(x).conIdx(idxKeep{x}),1:length(idxKeep),'uniformoutput',0);
            IDs = vertcat(IDs{:}); % note the IDs get messed up so the distValues are not in the same order as the filoBranch
            distValues = vertcat(distValues{:});
            
            
            distToBranchCell{iFrame} = distValues.*0.216;
            clear distValues
%         else
%             distToBranchCell{iFrame} = [];
%        % end
    else
        distToBranchCell{iFrame} = [];
    end
end

end





