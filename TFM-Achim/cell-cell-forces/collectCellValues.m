function varargout=collectCellValues(groupedClusters,goodCellSet,varargin)
% [deg_vals,elE_vals,area_vals,resF_vals]=collectCellValues(groupedClusters,goodCellSet,'deg','elE','area','resF')
% Runs through the groupedClusters and collects all the data from fields
% specified in the input arguments. Potential fields are:
% 'deg'  : The degree of connectivity of a cell.
% 'elE'  : The contraction/elastic energy of the cell.
% 'area' : The area of the cell.
% 'resF' : The residual force of the cell.
% 'intF' : The interfacial forces exerted on this cell.

degPos=find(strcmp('deg',varargin));
if ~isempty(degPos)
    degCheck = 1;
    % initialize:
    deg_vals = [];
else
    degCheck = 0;
end
    

elEPos=find(strcmp('elE',varargin));
if ~isempty(elEPos)
    elECheck = 1;
    % initialize:
    elE_vals   = [];
else
    elECheck = 0;
end

areaPos=find(strcmp('area',varargin));
if ~isempty(areaPos)
    areaCheck = 1;
    % initialize:
    area_vals   = [];
else
    areaCheck = 0;
end

resFPos=find(strcmp('resF',varargin));
if ~isempty(resFPos)
    resFCheck = 1;
    % initialize:
    resF_vals   = [];
else
    resFCheck = 0;
end

intFPos=find(strcmp('intF',varargin));
if ~isempty(intFPos)
    intFCheck = 1;
    % initialize:
    maxIdx =length(goodCellSet);
    intF_vals(maxIdx).edge = [];
    
%     maxEdgeNum=100;  % This can be read out of the data set!
%     % maxFrame  =200;  % This can be read out of the data set!
%     intF_vals(maxIdx).edge(maxEdgeNum).f1    = []; % Make sure to take the right direction!!!
%     intF_vals(maxIdx).edge(maxEdgeNum).f2    = []; % Make sure to take the right direction!!!
%     intF_vals(maxIdx).edge(maxEdgeNum).fn    = []; % Is the mean of f1 and f2. Make sure to take the right direction!!!
%     intF_vals(maxIdx).edge(maxEdgeNum).fc    = []; % Make sure to take the right direction!!!
%     intF_vals(maxIdx).edge(maxEdgeNum).frame = [];
else
    intFCheck = 0;
end


for idx=1:length(goodCellSet)
    clusterId=goodCellSet(idx).clusterId;
    cellId   =goodCellSet(idx).cellId;
    toDoList =goodCellSet(idx).frames;
    
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
        
        tic;
        if  intFCheck
            edges=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.edges;
            for edgeId=edges
                if edgeId>length(intF_vals(idx).edge)
                   intF_vals(idx).edge(edgeId).f1=[]; 
                   intF_vals(idx).edge(edgeId).f2=[];
                   intF_vals(idx).edge(edgeId).fn=[];
                   intF_vals(idx).edge(edgeId).fc=[];
                   intF_vals(idx).edge(edgeId).frames=[];
                end
                % pick the right direction:
                display('To do: pick the right direction!')
                f1=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f1;
                f2=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f2;
                fc=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc;

                intF_vals(idx).edge(edgeId).f1=vertcat(intF_vals(idx).edge(edgeId).f1,f1);
                intF_vals(idx).edge(edgeId).f2=vertcat(intF_vals(idx).edge(edgeId).f2,f2);
                intF_vals(idx).edge(edgeId).fn=vertcat(intF_vals(idx).edge(edgeId).fn,0.5*(f1-f2));
                intF_vals(idx).edge(edgeId).fc=vertcat(intF_vals(idx).edge(edgeId).fc,fc);
                intF_vals(idx).edge(edgeId).frames=vertcat(intF_vals(idx).edge(edgeId).frames,frame);
            end
            
        end
        toc;
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

if  intFCheck
    varargout(intFPos) = {intF_vals};
end

% Do some simple output checks:
if length(varargout)~=optargin
    error('Something went wrong');
end

for k=2:optargin
    if size(varargout{1},1)~=size(varargout{k},1)
        error('Something went wrong');
    end
end