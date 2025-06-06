function goodCellSet=findCells(groupedClusters,varargin)
% goodCellSet=findCells(groupedClusters,'deg',2,'kPa',10,'myo',[0 1],'errs',0)
% Optional constraints that can be included in the search through the
% cluster list. Only cells that match ALL search patterns will be regarded
% as a good data and given in the output.
% 'deg'    : List of deg-values. Cells with a degree that is not in the
%            given list will be dismissed in the search.
% 'kPa'    : List of stiffness values in [kPa]. Experiments on substrates of
%            different stiffness will be dismissed.
% 'myo'    : 0/1: 0 = control cell; 1 = myosin cell. Cells that don't match
%            the search pattern will be dismissed.
% 'type'   : Select for specifc myosinII hairpins (e.g. myoIIA_hp93,
%            myoIIA_hp94, myoIIB_hp103,tln1,ctrl). 'type' can be a list of patterns,
%            e.g. all myoIIA hair pins.
% 'myoGlb' : 0/1/-1: 0 = cluster with only control cells; 1 = cluster with only
%            myosin cells. -1 = mixed clusters. Clusters that don't match
%            the search pattern will be dismissed.
% 'divGlb' : 0/1/-1/NaN:
%            This will fail if one cell leaves and another one joins the
%            cluster in the same frame!
%            0 = cluster with constant cluster size: no divisions
%            or cell deaths. 
%            1 = cluster sizes increase at least once over time (but never
%            shrinks) by: 
%            a) a cell division or by 
%            b) an additional cell joining the cluster.
%            -1 = cluster sizes decrease at least once over time (but the
%            cluster never grows) by: 
%            a) a cell died
%            b) a cell left the cluster
%            [1, -1]: look for all clusters that either shrinks or grows
%            over time. In other words a cluster that at least once
%            changes size (in any direction).
% 'errF'   : x: Only clusters with an error in the force measurement of <x
%            will be considered (the magnitude of the error in [nN]).
% 'relErrF': x: Only clusters with a rel error in the force measurement of
%            <x will be considered. The rel error is calculated as 
%            the (error force)/(maximal cell-cell force found in the cluster).
%            thus [x]=1.
% 'errs'   : 0/x: If 0, only clusters with no missing forces will be
%            considered. If x>0 is given then, clusters with an errors<x will 
%            be considered.

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

myoGlbPos=find(strcmp('myoGlb',varargin));
if ~isempty(myoGlbPos)
    myoGlbCheck = 1;
    % it is the next entry which contains the numeric value:
    myoGlbVal   = varargin{myoGlbPos+1};
else
    myoGlbCheck = 0;
end

divGlbPos=find(strcmp('divGlb',varargin));
if ~isempty(divGlbPos)
    divGlbCheck = 1;
    % it is the next entry which contains the numeric value:
    divGlbVal   = varargin{divGlbPos+1};
else
    divGlbCheck = 0;
end

typePos=find(strcmp('type',varargin));
if ~isempty(typePos)
    typeCheck = 1;
    % it is the next entry which contains the numeric value:
    typeVal   = varargin{typePos+1};
else
    typeCheck = 0;
end

errsPos=find(strcmp('errs',varargin));
if ~isempty(errsPos)
    errsCheck = 1;
    % it is the next entry which contains the numeric value:
    errsVal   = varargin{errsPos+1};
else
    errsCheck = 0;
end


errFPos=find(strcmp('errF',varargin));
if ~isempty(errFPos)
    errFCheck = 1;
    % it is the next entry which contains the numeric value:
    errFVal   = varargin{errFPos+1};
else
    errFCheck = 0;
end

relErrFPos=find(strcmp('relErrF',varargin));
if ~isempty(relErrFPos)
    relErrFCheck = 1;
    % it is the next entry which contains the numeric value:
    relErrFVal   = varargin{relErrFPos+1};
else
    relErrFCheck = 0;
end

idx=1;
goodCellSet(idx).clusterId=[];
goodCellSet(idx).cellId   =[];
goodCellSet(idx).frames   =[];
for clusterId=1:groupedClusters.numClusters
    toDoList=[];
    numCellsVec=[];
    trackedNet=groupedClusters.cluster{clusterId}.trackedNet;
    maxCell=0;
    for frame=1:length(trackedNet)
        if ~isempty(trackedNet{frame})
            toDoList=horzcat(toDoList,frame);
            maxCell=max(maxCell,length(trackedNet{frame}.node));
            numCellsVec=vertcat(numCellsVec, trackedNet{frame}.stats.numCells);
        end
    end
    diffVec=numCellsVec(2:end)-numCellsVec(1:(end-1));
    
    checkVecInc=(diffVec== 1);
    checkSumInc= sum(checkVecInc);
    checkVecDec=(diffVec==-1);
    checkSumDec= sum(checkVecDec);
    
    % find subsequent frames with fewer than the maximum number of cells.
    if checkSumInc>0 && checkSumDec==0
        % then the cluster is at least growing once but never shrinks.
        currDivGlb=1;
    elseif checkSumInc==0 && checkSumDec>0
        % Then the cluster is shrinking at least once but never grows.
        currDivGlb=-1;
    elseif checkSumInc>0 && checkSumDec>0
        % then the cluster is shrinking and growing at least once
        currDivGlb=[-1,1];
    else
        % teh cluster is not changing size:
        currDivGlb=0;
    end
    
    for cellId=1:maxCell
        for frame=toDoList
            % check the composition of the network:
            if trackedNet{frame}.stats.numMyo==0;
                currMyoGlb=0;
            elseif trackedNet{frame}.stats.numMyo==trackedNet{frame}.stats.numCells
                currMyoGlb=1;
            else
                currMyoGlb=-1;
            end
            
            % find maximal cell-cell force in the cluster:
            maxCCF=0;
            for iEdge=1:length(trackedNet{frame}.edge)
                if ~isempty(trackedNet{frame}.edge{iEdge})
                    maxCCF=max([maxCCF, norm(trackedNet{frame}.edge{iEdge}.fc)]) ;
                end
            end
            % Now go through all checks:
            if     (cellId<=length(trackedNet{frame}.node))...
                    && (~isempty(trackedNet{frame}.node{cellId}))...
                    && (~kPaCheck     || ismember(trackedNet{frame}.par.yModu_Pa/1000,kPaVal))...
                    && (~degCheck     || ismember(trackedNet{frame}.node{cellId}.deg ,degVal))...
                    && (~myoCheck     || ismember(trackedNet{frame}.node{cellId}.spec,myoVal))...
                    && (~myoGlbCheck  || ismember(currMyoGlb,myoGlbVal))...
                    && (~typeCheck    || sum(strcmp(trackedNet{frame}.node{cellId}.type,typeVal)>0))...
                    && (~errsCheck    || trackedNet{frame}.stats.errs<=errsVal)...
                    && (~errFCheck    || trackedNet{frame}.stats.errorSumForce.mag<=errFVal)...
                    && (~relErrFCheck || trackedNet{frame}.stats.errorSumForce.mag/maxCCF<=relErrFVal)...
                    && (~divGlbCheck  || numel(setdiff(currDivGlb,divGlbVal))==0)
                
                
                % If all of this is true we have found a good entry
                % display(['cluster: ',num2str(clusterId),' cell: ',num2str(cellId),' frame: ',num2str(frame)]);
                
                goodCellSet(idx).clusterId=clusterId;
                goodCellSet(idx).cellId   =cellId;
                goodCellSet(idx).frames   =[goodCellSet(idx).frames,frame];
            end
            
        end
        idx=idx+1;
    end
end

% sort out the empty ones:
setId=1;
while setId<=length(goodCellSet)
    if isempty(goodCellSet(setId).frames)
        goodCellSet(setId)    =[];
    else
        setId=setId+1;
    end
end
