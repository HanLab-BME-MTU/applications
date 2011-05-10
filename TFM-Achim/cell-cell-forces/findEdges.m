function goodEdgeSet=findEdges(groupedClusters,varargin)
% goodEdgeSet=findEdges(groupedClusters,'deg',[1 2],'kPa',10,'myo',[1 1],'type',{'myoIIA_hp93';'myoIIA_hp94'},'errs',0)
% Optional constraints that can be included in the search through the
% cluster list. Only cells that match ALL search patterns will be regarded
% as a good data and given in the output.
% 'deg'   : List of deg-values. Since edges connect two cells, 'deg' should
%           be at least a 2D-vector. If only one value x is given, it is
%           assumed that deg=[x x]. Edges with a degree that is not in the
%           given list will be dismissed in the search.
% 'kPa'   : List of stiffness values in [kPa]. Experiments on substrates of
%           different stiffness will be dismissed.
% 'myo'   : 0/1/-1: The edge connects: 
%                   0 = two control cells; 
%                   1 = two myosin cells;
%                  -1 = a mixture.
%           'myo' can be a list of values. Edges that don't match the
%           search pattern will be dismissed.
% 'type'  : Select for specifc myosinII hairpins e.g. {'myoIIA_hp93';
%           'myoIIA_hp94'; 'myoIIB_hp103'}. 'type' can be a list of patterns,
%           e.g. all myoIIA hair pins. This makes only sense if it is
%           required that the edge connects to at least one myosin cell.
%           Therefore this check will be performed only on the myosin cells
%           not on the control cells. Thus, if an edge that connects to a
%           myosin cell that is not of a type found in the list, then the
%           edge will be dismissed. An edge connecting two control cells
%           will always pass!!!
% 'myoGlb': 0/1/-1: 0 = cluster with only control cells; 1 = cluster with
%           only myosin cells. -1 = mixed clusters. Clusters that don't
%           match the search pattern will be dismissed.
% 'asmbly': Find edges that assemble or disassemble over time:
%            0: edge is there over the whole length of the movie
%            1: edge is formed during the movie
%           -1: edge that disassemble during the movie (cell is breaking
%               apart)
% 'dItotRel': Threshold that the relative change in intensity (Itot) has to
%           exceed to pass the test.
% 'errF'  : x: Only clusters with an error in the force measurement of <x
%           will be considered (the magnitude of the error in [nN]).
% 'errs'  : 0/x: If 0, only clusters with no missing forces will be 
%           considered. If x>0 is given then, clusters with an errors < x
%           will be considered.

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

asmblyPos=find(strcmp('asmbly',varargin));
if ~isempty(asmblyPos)
    asmblyCheck = 1;
    % it is the next entry which contains the numeric value:
    asmblyVal   = varargin{asmblyPos+1};
else
    asmblyCheck = 0;
end

dItotRelPos=find(strcmp('dItotRel',varargin));
if ~isempty(dItotRelPos)
    dItotRelCheck = 1;
    % it is the next entry which contains the numeric value:
    dItotRelVal   = varargin{dItotRelPos+1};
else
    dItotRelCheck = 0;
end

typePos=find(strcmp('type',varargin));
if ~isempty(typePos)
    typeCheck = 1;
    % it is the next entry which contains the numeric value:
    typeVal   = varargin{typePos+1};
    % make sure that a [ctrl ctrl] will always pass the test:
    typeVal   = vertcat(typeVal, 'ctrl');
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


idx=1;
goodEdgeSet(idx).clusterId=[];
goodEdgeSet(idx).edgeId   =[];
goodEdgeSet(idx).frames   =[];
for clusterId=1:groupedClusters.numClusters
    toDoList=[];
    trackedNet=groupedClusters.cluster{clusterId}.trackedNet;
    maxEdge=0;
    for frame=1:length(trackedNet)
        if ~isempty(trackedNet{frame})
            toDoList=horzcat(toDoList,frame);
            maxEdge=max(maxEdge,length(trackedNet{frame}.edge));
        end
    end
    
    
    for edgeId=1:maxEdge
        if asmblyCheck
            k=1;
            for frame=toDoList
                % Check if the edge exists in the frame:
                if length(trackedNet{frame}.edge)<edgeId
                    checkVec(k,1)=0;
                else
                    checkVec(k,1)=~isempty(trackedNet{frame}.edge{edgeId});
                end
                k=k+1;
            end
            % if this vector is 0, the edge is there, if it is 1 the edge
            % was formed, if it -1 the edge is lost:
            asmblyVec=checkVec(2:end)-checkVec(1:(end-1));
            
            checkVecAss=(asmblyVec== 1);
            checkSumAss= sum(checkVecAss);
            checkVecDis=(asmblyVec==-1);
            checkSumDis= sum(checkVecDis);
            
            if checkSumAss>0 && checkSumDis==0
                % then the edge was formed over time.
                currAssVal=1;
            elseif checkSumAss==0 && checkSumDis>0
                % then the edge was lost over time.
                currAssVal=-1;
            elseif checkSumAss>0 && checkSumDis>0
                % then the edge was formed and lost during the movie
                currAssVal=[-1,1];
            else
                % teh cluster is not changing size:
                currAssVal=0;
            end
            clear checkVec
        end
        
        for frame=toDoList
            % check the composition of the network:
            if trackedNet{frame}.stats.numMyo==0;
                currMyoGlb=0;
            elseif trackedNet{frame}.stats.numMyo==trackedNet{frame}.stats.numCells
                currMyoGlb=1;
            else
                currMyoGlb=-1;
            end
            
            % get the degree pair of this edge if needed:
            if degCheck && edgeId<=length(trackedNet{frame}.edge) && ~isempty(trackedNet{frame}.edge{edgeId})
                nodes=trackedNet{frame}.edge{edgeId}.nodes;                
                degPair=[trackedNet{frame}.node{nodes(1)}.deg trackedNet{frame}.node{nodes(2)}.deg];                
            end
            
            % get the myo-mixture value for this edge if needed
            if myoCheck && edgeId<=length(trackedNet{frame}.edge) && ~isempty(trackedNet{frame}.edge{edgeId})
                nodes=trackedNet{frame}.edge{edgeId}.nodes;
                if trackedNet{frame}.node{nodes(1)}.spec==1 && trackedNet{frame}.node{nodes(2)}.spec==1
                    myoEdgeVal =1;
                elseif trackedNet{frame}.node{nodes(1)}.spec==0 && trackedNet{frame}.node{nodes(2)}.spec==0
                    myoEdgeVal =0;
                else
                    myoEdgeVal=-1;
                end
            end
            
            % get the myo-mixture value for this edge if needed
            if typeCheck && edgeId<=length(trackedNet{frame}.edge) && ~isempty(trackedNet{frame}.edge{edgeId})
                myoEdgeType='ctrl';
                nodes=trackedNet{frame}.edge{edgeId}.nodes;
                if trackedNet{frame}.node{nodes(1)}.spec==1
                    myoEdgeType =trackedNet{frame}.node{nodes(1)}.type;
                end
                % overwrite (they have to be the same anyways) with the
                % next value if any: 
                if trackedNet{frame}.node{nodes(2)}.spec==1
                    myoEdgeType =[trackedNet{frame}.node{nodes(2)}.type];
                end
            end
            
            % Now go through all checks:
            if     (edgeId<=length(trackedNet{frame}.edge))...
                    && (~isempty(trackedNet{frame}.edge{edgeId}))...
                    && (~kPaCheck    || ismember(trackedNet{frame}.par.yModu_Pa/1000,kPaVal))...
                    && (~degCheck    || isempty(setdiff(degPair,degVal)))...
                    && (~myoCheck    || ismember(myoEdgeVal,myoVal))...
                    && (~myoGlbCheck || ismember(currMyoGlb,myoGlbVal))...
                    && (~typeCheck   || sum(strcmp(myoEdgeType,typeVal))>0)...
                    && (~errsCheck   || trackedNet{frame}.stats.errs<=errsVal)...
                    && (~errFCheck   || trackedNet{frame}.stats.errorSumForce.mag<=errFVal)...
                    && (~asmblyCheck || numel(setdiff(currAssVal,asmblyVal))==0)                               
                
                
                % If all of this is true we have found a good entry
                % display(['cluster: ',num2str(clusterId),' cell: ',num2str(edgeId),' frame: ',num2str(frame)]);
                
                goodEdgeSet(idx).clusterId=clusterId;
                goodEdgeSet(idx).edgeId   =edgeId;
                goodEdgeSet(idx).frames   =[goodEdgeSet(idx).frames,frame];
            end            
        end
        idx=idx+1;
    end
end

% sort out the empty ones:
setId=1;
while setId<length(goodEdgeSet)
    if isempty(goodEdgeSet(setId).frames)
        goodEdgeSet(setId)    =[];
    else
        setId=setId+1;
    end
end

% sort out the ones that do not change significantly in intensity:
if dItotRelCheck
    setId=1;
    while setId<=length(goodEdgeSet)
        edgeId    = goodEdgeSet(setId).edgeId;
        toDoList  = goodEdgeSet(setId).frames;
        clusterId = goodEdgeSet(setId).clusterId;
        k=1;
        for frame=toDoList
            ItotList(k)=groupedClusters.cluster{clusterId}.trackedNet{frame}.edge{edgeId}.int.tot;
            k=k+1;
        end
        [maxItot idMax]=max(ItotList);
        [minItot idMin]=min(ItotList);
        
        
        curr_dItotRel=maxItot/minItot;
        if (sign(dItotRelVal)~=sign(idMax-idMin)) || (sign(dItotRelVal)==sign(idMax-idMin) && curr_dItotRel<abs(dItotRelVal))
            goodEdgeSet(setId)    =[];
        else
            setId=setId+1;
        end
        clear ItotList
    end
    
end
