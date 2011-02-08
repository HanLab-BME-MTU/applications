function varargout=collectCellValues(groupedClusters,goodCellSet,varargin)
% [deg_out,elE_out,area_out,resF_out,intF_out]=collectCellValues(groupedClusters,goodCellSet,'deg','elE','area','resF','corr')
% [corr_vals]=collectCellValues(groupedClusters,goodCellSet,'corr')
% Runs through the groupedClusters and collects all the data from fields
% specified in the input arguments. Potential fields are:
% 'deg'  : The degree of connectivity of a cell.
% 'elE'  : The contraction/elastic energy of the cell.
% 'area' : The area of the cell.
% 'resF' : The residual force of the cell.
% 'corr' : The interfacial forces exerted on this cell.

degPos=find(strcmp('deg',varargin));
if ~isempty(degPos)
    degCheck = 1;
else
    degCheck = 0;
end
% initialize:
deg_out = [];
    

elEPos=find(strcmp('elE',varargin));
if ~isempty(elEPos)
    elECheck = 1;
else
    elECheck = 0;
end
% initialize:
elE_out  = [];

areaPos=find(strcmp('area',varargin));
if ~isempty(areaPos)
    areaCheck = 1;
else
    areaCheck = 0;
end
% initialize:
area_out   = [];

resFPos=find(strcmp('resF',varargin));
if ~isempty(resFPos)
    resFCheck = 1;
else
    resFCheck = 0;
end
% initialize:
resF_out   = [];

corrPos=find(strcmp('corr',varargin));
if ~isempty(corrPos)
    corrCheck = 1;
    % initialize:
    maxIdx =length(goodCellSet);
    corr_out(maxIdx).edge = [];
    
%     maxEdgeNum=100;  % This can be read out of the data set!
%     % maxFrame  =200;  % This can be read out of the data set!
%     corr_out(maxIdx).edge(maxEdgeNum).f1    = []; % Make sure to take the right direction!!!
%     corr_out(maxIdx).edge(maxEdgeNum).f2    = []; % Make sure to take the right direction!!!
%     corr_out(maxIdx).edge(maxEdgeNum).fn    = []; % Is the mean of f1 and f2. Make sure to take the right direction!!!
%     corr_out(maxIdx).edge(maxEdgeNum).fc    = []; % Make sure to take the right direction!!!
%     corr_out(maxIdx).edge(maxEdgeNum).frame = [];
else
    corrCheck = 0;
end


for idx=1:length(goodCellSet)
    clusterId=goodCellSet(idx).clusterId;
    cellId   =goodCellSet(idx).cellId;
    toDoList =goodCellSet(idx).frames;
    
    % initialize all values:
    deg_vals =NaN+zeros(toDoList(end),1);
    elE_vals =NaN+zeros(toDoList(end),1);
    area_vals=NaN+zeros(toDoList(end),1);
    resF_vals=NaN+zeros(toDoList(end),2);
    
    for frame=toDoList
        % Now go through all checks:
        if  degCheck
            deg_vals(frame)=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.deg;
        end
        
        if  elECheck
            if isempty(groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.elE)
                elE_vals(frame)=NaN;
            else
                elE_vals(frame)=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.elE;
            end
        end
        
        if  areaCheck
            area_vals(frame)=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.area;
        end
        
        if  resFCheck
            resF_vals(frame,:)=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.vec;
        end
        
        if  corrCheck
            edges=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.edges;
            for edgeId=edges
                if edgeId>length(corr_out(idx).edge)
                   % set all initial force values to NaNs:
                   corr_out(idx).edge(edgeId).f1=NaN+zeros(toDoList(end),2); 
                   corr_out(idx).edge(edgeId).f2=NaN+zeros(toDoList(end),2);
                   corr_out(idx).edge(edgeId).fn=NaN+zeros(toDoList(end),2);
                   corr_out(idx).edge(edgeId).fc=NaN+zeros(toDoList(end),2);
                   % The frame list runs from 1 to the maximum of the
                   % toDoList:
                   corr_out(idx).edge(edgeId).frames=1:toDoList(end);
                   % just to keep track, input the original edgeId:
                   corr_out(idx).edge(edgeId).edgeId=edgeId;
                end
                % display('To do: pick the right direction!')
                f1 =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f1;
                f2 =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f2;
                fc1=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc1;
                
                % find the position of the current cell/node and pick the
                % right direction! fn and fc are the force exerted by the
                % neighboring cells on to the current cell:
                if find(groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.nodes==cellId)==1
                    fh = f1;
                    f1 = f2;
                    f2 = fh;
                    fc1=-fc1;
                elseif isempty(find(groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.nodes==cellId, 1))
                    error('Something went wrong')
                end
                corr_out(idx).edge(edgeId).f1(frame,:)=f1;
                corr_out(idx).edge(edgeId).f2(frame,:)=f2;
                corr_out(idx).edge(edgeId).fn(frame,:)=0.5*(f1-f2);
                corr_out(idx).edge(edgeId).fc(frame,:)=fc1;
                % corr_out(idx).edge(edgeId).frames=vertcat(corr_out(idx).edge(edgeId).frames,frame);
            end
            
        end
    end
    % append the found values:
    deg_out = vertcat( deg_out, deg_vals);
    elE_out = vertcat( elE_out, elE_vals);
    area_out= vertcat(area_out,area_vals);
    resF_out= vertcat(resF_out,resF_vals);
    
    %remove empty edges
    if  corrCheck
        corr_out(idx).resF=resF_vals;
        edgeId=1;
        while edgeId<=length(corr_out(idx).edge)
            checkVal=(isempty(corr_out(idx).edge(edgeId)) || isempty(corr_out(idx).edge(edgeId).edgeId));
            if checkVal
                corr_out(idx).edge(edgeId)=[];
            else
                edgeId=edgeId+1;
            end
        end
    end
end

% number of optional input = number of output arguments:
optargin = size(varargin,2);


if degCheck
    varargout(degPos)  = {deg_out};
end

if elECheck
    varargout(elEPos)  = {elE_out};
end

if areaCheck
    varargout(areaPos) = {area_out};
end

if  resFCheck
    varargout(resFPos) = {resF_out};
end

if  corrCheck
    varargout(corrPos) = {corr_out};
end

% Do some simple output checks:
if length(varargout)~=optargin
    error('Something went wrong');
end

for k=2:optargin
    if size(varargout{1},1)~=size(varargout{k},1)
        display('Something might have went wrong');
    end
end