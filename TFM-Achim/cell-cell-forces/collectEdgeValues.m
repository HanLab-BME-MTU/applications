function varargout=collectEdgeValues(groupedClusters,goodSet,varargin)
% [deg_vals,lgth_vals,f1_vals,f2_vals,fc1_vals,fc2_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'deg','lgth','f1','f2','fc1','fc2')
% Runs through the groupedClusters and collects all the data from fields
% specified in the input arguments. Potential fields are:
% 'deg'  : The degree of connectivity of a cell.
% 'lgth' : The interface length.
% 'f1'   : The force f1 determined by the network analysis
% 'f2'   : The force f2 determined by the network analysis
% 'fc1'  : The force f1 determined by the network analysis
% 'fc2'  : The force f2 determined by the network analysis
% 'Itot' : The total Ecad intensity along the interface
% 'Iavg' : The total Ecad intensity along the interface

degPos=find(strcmp('deg',varargin));
if ~isempty(degPos)
    degCheck = 1;
    deg_vals = [];
else
    degCheck = 0;
end
    

lgthPos=find(strcmp('lgth',varargin));
if ~isempty(lgthPos)
    lgthCheck = 1;
    % it is the next entry which contains the numeric value:
    lgth_vals   = [];
else
    lgthCheck = 0;
end

f1Pos=find(strcmp('f1',varargin));
if ~isempty(f1Pos)
    f1Check = 1;
    % it is the next entry which contains the numeric value:
    f1_vals   = [];
else
    f1Check = 0;
end

f2Pos=find(strcmp('f2',varargin));
if ~isempty(f2Pos)
    f2Check = 1;
    % it is the next entry which contains the numeric value:
    f2_vals   = [];
else
    f2Check = 0;
end

fc1Pos=find(strcmp('fc1',varargin));
if ~isempty(fc1Pos)
    fc1Check = 1;
    % it is the next entry which contains the numeric value:
    fc1_vals   = [];
else
    fc1Check = 0;
end

fc2Pos=find(strcmp('fc2',varargin));
if ~isempty(fc2Pos)
    fc2Check = 1;
    % it is the next entry which contains the numeric value:
    fc2_vals   = [];
else
    fc2Check = 0;
end

ItotPos=find(strcmp('Itot',varargin));
if ~isempty(ItotPos)
    ItotCheck = 1;
    % it is the next entry which contains the numeric value:
    Itot_vals   = [];
else
    ItotCheck = 0;
end

IavgPos=find(strcmp('Iavg',varargin));
if ~isempty(IavgPos)
    IavgCheck = 1;
    % it is the next entry which contains the numeric value:
    Iavg_vals   = [];
else
    IavgCheck = 0;
end

for idx=1:length(goodSet)
    clusterId=goodSet(idx).clusterId;
    edgeId   =goodSet(idx).edgeId;
    toDoList =goodSet(idx).frames;
    
    for frame=toDoList
        % Now go through all checks:
        if  degCheck
            nodes=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.nodes;
            deg_vals=vertcat(deg_vals,[groupedClusters.cluster{clusterId}.trackedNet{frame}.node{nodes(1)}.deg groupedClusters.cluster{clusterId}.trackedNet{frame}.node{nodes(2)}.deg]);
        end
        
        if  lgthCheck
            lgth_vals=vertcat(lgth_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.intf_internal_L);
        end
        
        if  f1Check
            f1_vals=vertcat(f1_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f1);
        end
        
        if  f2Check
            f2_vals=vertcat(f2_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.f2);
        end
        
        if  fc1Check
            fc1_vals=vertcat(fc1_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc1);
        end
        
        if  fc2Check
            fc2_vals=vertcat(fc2_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.fc2);
        end
        
        if  ItotCheck
            Itot_vals=vertcat(Itot_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.tot);
        end
        
        if  IavgCheck
            Iavg_vals=vertcat(Iavg_vals ,groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.avg);
        end
    end
end

% number of optional input = number of output arguments:
optargin = size(varargin,2);


if degCheck
    varargout(degPos)  = {deg_vals};
end

if lgthCheck
    varargout(lgthPos)  = {lgth_vals};
end

if f1Check
    varargout(f1Pos) = {f1_vals};
end

if f2Check
    varargout(f2Pos) = {f2_vals};
end

if fc1Check
    varargout(fc1Pos) = {fc1_vals};
end

if fc2Check
    varargout(fc2Pos) = {fc2_vals};
end

if ItotCheck
    varargout(ItotPos) = {Itot_vals};
end

if IavgCheck
    varargout(IavgPos) = {Iavg_vals};
end


if length(varargout)~=optargin
    error('Something went wrong');
end

for k=2:optargin
    if size(varargout{1},1)~=size(varargout{k},1)
        display('Something might have went wrong, vectors are of different length!');
    end
end