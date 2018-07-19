function [ filoBranchComplexityCell] = GCAAnalysisExtract_filoBranchComplexity( analInfo,filoFilterSet,varargin )
%%% GCAAnalysisExtract_filoBranchComplexity 


%% Check Parameters 

ip = inputParser;
ip.KeepUnmatched = true;

ip.CaseSensitive = false;

ip.addParameter('umPerPixel',.216);

ip.parse(varargin{:});

filoBranchComplexityCell = cell(length(analInfo)-1,1);

for iFrame = 1:length(analInfo) -1
    
    filoInfo = analInfo(iFrame).filoInfo;
    
    filterFrameC= filoFilterSet{iFrame};
    if ~isempty(filterFrameC)
        filoInfoFilt  = filoInfo(filterFrameC(:,1)); % after this filter should have only N order branch and its corresponding N-1 branchstem
        %% get the branch stems
        types = vertcat(filoInfoFilt(:).type);
        NBranches = sum(types>1); % get the total number of branches
        
        lengths = vertcat(filoInfoFilt(:).Ext_length);
        totalLength = sum(lengths(~isnan(lengths)));
        totalLength = totalLength.*ip.Results.umPerPixel; 
        
        filoBranchComplexityCell{iFrame,1} = NBranches./totalLength*10; %  output per 10 um 
        
    else
        filoBranchComplexityCell{iFrame,1} = [];
    end
     
end


