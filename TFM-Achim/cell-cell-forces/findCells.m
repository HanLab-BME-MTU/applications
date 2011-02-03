function goodSet=findCells(groupedClusters,varargin)
% goodSet=findCells(groupedClusters,'deg',2,'kPa',10,'myo',[0 1])
% constraints that can be included in the search through the cluster list:
% 'deg'
% 'kPa'
% 'myo'

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

myoPos=find(strcmp('myo',varargin));
if ~isempty(myoPos)
    myoCheck = 1;
    % it is the next entry which contains the numeric value:
    myoVal   = varargin{myoPos+1};
else
    myoCheck = 0;
end


idx=1;
goodSet(idx).clusterId=[];
goodSet(idx).cellId   =[];
goodSet(idx).frames   =[];
for clusterId=1:groupedClusters.numClusters
    toDoList=[];
    trackedNet=groupedClusters.cluster{clusterId}.trackedNet;
    maxCell=0;
    for frame=1:length(trackedNet)
        if ~isempty(trackedNet{frame})
            toDoList=horzcat(toDoList,frame);
            maxCell=max(maxCell,length(trackedNet{frame}.node));
        end
    end
    
    
    for cellId=1:maxCell
        for frame=toDoList
            % Now go through all checks:
            if     (cellId<=length(trackedNet{frame}.node))...
                    && (~isempty(trackedNet{frame}.node{cellId}))...
                    && (~kPaCheck || ismember(trackedNet{frame}.par.yModu_Pa/1000,kPaVal))...
                    && (~degCheck || ismember(trackedNet{frame}.node{cellId}.deg ,degVal))...
                    && (~myoCheck || ismember(trackedNet{frame}.node{cellId}.spec,myoVal))
                
                % If all of this is true we have found a good entry
                display(['cluster: ',num2str(clusterId),' cell: ',num2str(cellId),' frame: ',num2str(frame)]);
                
                goodSet(idx).clusterId=clusterId;
                goodSet(idx).cellId   =cellId;
                goodSet(idx).frames   =[goodSet(idx).frames,frame];
            end
            
        end
        idx=idx+1;
    end
end

% sort out the empty ones:
setId=1;
while setId<length(goodSet)
    if isempty(goodSet(setId).frames)
        goodSet(setId)    =[];
    else
        setId=setId+1;
    end
end