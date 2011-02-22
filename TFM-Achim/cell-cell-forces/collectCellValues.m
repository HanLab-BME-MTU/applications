function varargout=collectCellValues(groupedClusters,goodCellSet,varargin)
% [deg_out,elE_out,area_out,sumFmag_out,resF_out,intF_out]=collectCellValues(groupedClusters,goodCellSet,'deg','elE','area','sumFmag','resF','corr')
% [corr_vals]=collectCellValues(groupedClusters,goodCellSet,'corr')
% Runs through the groupedClusters and collects all the data from fields
% specified in the input arguments. Potential fields are:
% 'deg'    : The degree of connectivity of a cell.
% 'elE'    : The contraction/elastic energy of the cell.
% 'area'   : The area of the cell.
% 'sumFmag': The sum of the force magnitudes over the footprint of a cell.
% 'resF'   : The residual force of the cell.
% 'sumFi'  : The sum of interfacial forces (magnitude) exerted on the cell.
% 'sumLi'  : The total length of the cell's interfaces.
% 'corr'   : The interfacial forces exerted on this cell.

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

sumFmagPos=find(strcmp('sumFmag',varargin));
if ~isempty(sumFmagPos)
    sumFmagCheck = 1;
else
    sumFmagCheck = 0;
end
% initialize:
sumFmag_out  = [];

resFPos=find(strcmp('resF',varargin));
if ~isempty(resFPos)
    resFCheck = 1;
else
    resFCheck = 0;
end
% initialize:
resF_out   = [];

sumFiPos=find(strcmp('sumFi',varargin));
if ~isempty(sumFiPos)
    sumFiCheck = 1;
else
    sumFiCheck = 0;
end
% initialize:
sumFi_out   = [];

sumLiPos=find(strcmp('sumLi',varargin));
if ~isempty(sumLiPos)
    sumLiCheck = 1;
else
    sumLiCheck = 0;
end
% initialize:
sumLi_out   = [];

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
    deg_vals     = NaN+zeros(toDoList(end),1);
    elE_vals     = NaN+zeros(toDoList(end),1);
    area_vals    = NaN+zeros(toDoList(end),1);
    sumFmag_vals = NaN+zeros(toDoList(end),1);
    resF_vals    = NaN+zeros(toDoList(end),2);
    sumFi_vals   = NaN+zeros(toDoList(end),1);
    sumLi_vals   = NaN+zeros(toDoList(end),1);
    
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
        
        if  sumFmagCheck
            sumFmag_vals(frame)=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.sumFmag;
        end
        
        if  resFCheck || corrCheck
            resF_vals(frame,:)=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.vec;
        end
        
        if  sumFiCheck
            edges=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.edges;
            % Since we are only interested in the magnitude here, it
            % doesn't matter which direction we pick:
            sumFi=0;
            for edgeId=edges
                f1 =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f1;
                f2 =groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f2;
                fc1=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc1;
                fn_mag=sqrt(sum((0.5*(f1-f2)).^2,2));
                fc_mag=sqrt(sum(fc1.^2,2));
                
                % Take the cluster force only if networkforce doesn't
                % exist?!
                if ~isnan(fn_mag)
                    sumFi=sumFi+fn_mag;
                elseif ~isnan(fc_mag)
                    sumFi=sumFi+fc_mag;
                else
                    error(['Couldn''t find a force value for: cluster: ',num2str(clusterId),' edge: ',num2str(edgeId),' frame: ',num2str(frame)])
                    % This might happen if the edge is empty! The program
                    % should treat this case correctly. The error message
                    % can be removed when checked that everything is done
                    % correctly!
                end
            end
            sumFi_vals(frame,:)=sumFi;
        end
        
        if  sumLiCheck
            edges=groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.edges;
            sumLi=0;
            for edgeId=edges
                if ~isempty(groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.intf_internal_L)
                    sumLi = sumLi+groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.intf_internal_L;
                else
                    error(['Couldn''t find a length value for: cluster: ',num2str(clusterId),' edge: ',num2str(edgeId),' frame: ',num2str(frame)])
                    % This might happen if the edge is empty! The program
                    % should treat this case correctly. The error message
                    % can be removed when checked that everything is done
                    % correctly!
                end
            end
            sumLi_vals(frame,:)=sumLi;
        end
        
        if  corrCheck
            edges=sort(groupedClusters.cluster{clusterId}.trackedNet{frame}.node{cellId}.edges);
            for edgeId=edges
                if edgeId>length(corr_out(idx).edge)
                   % set all initial force values to NaNs:
                   corr_out(idx).edge(edgeId).fc   = NaN+zeros(toDoList(end),2);
                   corr_out(idx).edge(edgeId).fm   = NaN+zeros(toDoList(end),2);
                   corr_out(idx).edge(edgeId).flag = ones(toDoList(end),1);
                   % corr_out(idx).edge(edgeId).fn=NaN+zeros(toDoList(end),2);
                   % As I have set it up right now, fn cannot be used for
                   % the correlation, since an empty edge (created by
                   % trackNetwork) can not be distinguished from an edge
                   % with undetermined value of the network force. Since fc
                   % is determined for every edge, .fc=[NaN NaN]; means
                   % that this edge is actually empty. Since fi_tot is
                   % calculated as nansum, this value won't contribute as
                   % intended. This doesn't work with .fn=[NaN NaN]; which
                   % is ambiguous.
                   % To prevent misuse don't store these values in the
                   % corr_out structure!:
                   % corr_out(idx).edge(edgeId).f1=NaN+zeros(toDoList(end),2); 
                   % corr_out(idx).edge(edgeId).f2=NaN+zeros(toDoList(end),2);
                   
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
                % It is not save to use these values for reasons described
                % above:
                % corr_out(idx).edge(edgeId).f1(frame,:)=f1;
                % corr_out(idx).edge(edgeId).f2(frame,:)=f2;
                % corr_out(idx).edge(edgeId).fn(frame,:)=0.5*(f1-f2);
                fn = 0.5*(f1-f2);
                
                corr_out(idx).edge(edgeId).fc(frame,:)=fc1;
                % create the mixed force vector which takes the network
                % force whenever possible.
                if ~isnan(sum(fn))
                    corr_out(idx).edge(edgeId).fm(frame,:)=fn;
                    % to keep track which value we have chosen:
                    corr_out(idx).edge(edgeId).flag(frame,:)=1;
                else
                    corr_out(idx).edge(edgeId).fm(frame,:)=corr_out(idx).edge(edgeId).fc(frame,:);
                    % to keep track which value we have chosen:
                    if ~isnan(sum(corr_out(idx).edge(edgeId).fc(frame,:)))
                        corr_out(idx).edge(edgeId).flag(frame,:)=0;
                    else
                        % then, the edge was empty and it is OK to use the
                        % network result
                        corr_out(idx).edge(edgeId).flag(frame,:)=1;
                    end
                end
            end
        end
    end
    % append the found values:
    deg_out     = vertcat( deg_out, deg_vals);
    elE_out     = vertcat( elE_out, elE_vals);
    area_out    = vertcat(area_out,area_vals);
    sumFmag_out = vertcat(sumFmag_out,sumFmag_vals);
    resF_out    = vertcat(resF_out,resF_vals);
    sumFi_out   = vertcat(sumFi_out,sumFi_vals);
    sumLi_out   = vertcat(sumLi_out,sumLi_vals);
    
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

if sumFmagCheck
    varargout(sumFmagPos) = {sumFmag_out};
end

if  resFCheck
    varargout(resFPos) = {resF_out};
end

if  sumFiCheck
    varargout(sumFiPos) = {sumFi_out};
end

if  sumLiCheck
    varargout(sumLiPos) = {sumLi_out};
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