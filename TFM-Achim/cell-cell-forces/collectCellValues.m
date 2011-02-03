function varargout=collectCellValues(groupedClusters,goodSet,varargin)
% [deg_vals,elE_vals,area_vals,resF_vals]=collectCellValues(groupedClusters,goodSet,'deg','elE','area','resF')
% Runs through all cells and collects the data from fields:
% 'deg'
% 'elE'
% 'area'
% 'resF'

degPos=find(strcmp('deg',varargin));
if ~isempty(degPos)
    degCheck = 1;
    deg_vals = [];
else
    degCheck = 0;
end
    

elEPos=find(strcmp('elE',varargin));
if ~isempty(elEPos)
    elECheck = 1;
    % it is the next entry which contains the numeric value:
    elE_vals   = [];
else
    elECheck = 0;
end

areaPos=find(strcmp('area',varargin));
if ~isempty(areaPos)
    areaCheck = 1;
    % it is the next entry which contains the numeric value:
    area_vals   = [];
else
    areaCheck = 0;
end

resFPos=find(strcmp('resF',varargin));
if ~isempty(resFPos)
    resFCheck = 1;
    % it is the next entry which contains the numeric value:
    resF_vals   = [];
else
    resFCheck = 0;
end


for idx=1:length(goodSet)
    clusterId=goodSet(idx).clusterId;
    cellId   =goodSet(idx).cellId;
    toDoList =goodSet(idx).frames;
    
    for frame=toDoList
        % Now go through all checks:
        if  degCheck
            deg_vals=vertcat(deg_vals,groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.deg);
        end
        
        if  elECheck
            if isempty(groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.elE)
                elE_vals=vertcat(elE_vals ,NaN);
            else
                elE_vals=vertcat(elE_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.elE);
            end
        end
        
        if  areaCheck
            area_vals=vertcat(area_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.area);
        end
        
        if  resFCheck
            resF_vals=vertcat(resF_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.vec);
        end
    end
end

% number of optional input = number of output arguments:
optargin = size(varargin,2);


if degCheck
    varargout(degPos)  = {deg_vals};
end

if elECheck
    varargout(elEPos)  = {elE_vals};
end

if areaCheck
    varargout(areaPos) = {area_vals};
end

if  resFCheck
    varargout(resFPos) = {resF_vals};
end

if length(varargout)~=optargin
    error('Something went wrong');
end

for k=2:optargin
    if size(varargout{1},1)~=size(varargout{k},1)
        error('Something went wrong');
    end
end