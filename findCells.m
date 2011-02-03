function cellList=findCells(groupedClusters,varargin)
% constraints that can be included in the search through the cluster list:
% 'deg'
% 'kPa'

optargin = size(varargin,2);

% strcmp(fnameFirstBeadImg,groupedClusters.clusterList)

degPos=find(strcmp('deg',varargin));
if ~isempty(degPos)
    degCheck = 1;
    % it is the next entry which contains the numeric value:
    degVal   = varargin{degPos+1};
else
    degCheck = 0;
end
    

kPaPos=find(strcmp('kPa',varargin));
if ~isempty(kPaPos)
    kPaCheck = 1;
    % it is the next entry which contains the numeric value:
    kPaVal   = varargin{kPaPos+1};
else
    kPaCheck = 0;
end

for clusterId=1:groupedClusters.numClusters
    
end