function [corrResultsXY, corrResultsProjected]=TFM_part_8_corrAnalysis(constrForceField)

load('fileAndFolderNames.mat')

if ~strcmp(pwd,path_ProjFolder)
    display('Before running this script browse to the FSM project folder')
    return
end

if nargin < 1 || isempty(constrForceField)
    display('Loading cellCellForces.mat. This file is large and thus takes some while:...')
    tic;
    filestruct=load(path_cellCellForces);
    constrForceField=filestruct.constrForceField;
    toc;    
end

%**************************************************************************
% 1) Track cells and interfaces in the network and create network_tracked *
%**************************************************************************
[constrForceField]=trackNetwork(constrForceField);


%**************************************************************************
% 2) Pick the cell of interest:                                           *
%**************************************************************************
goodFrame=0;
while goodFrame==0
    frameToPlot=input('Enter the frame number you are interested: ');
    if isfield(constrForceField{frameToPlot},'network_tracked')
        plotCellNetwork(constrForceField{frameToPlot}.network_tracked)
        goodFrame=input('Good frame 1/0, [1]: ');
    else
        display('This was not a good frame!');
        goodFrame=0;
    end
    
    if isempty(goodFrame)
        goodFrame=1;
    end
end
cellNodeID=input('Enter the cell/node ID number you are interested in: ');

initDeg=constrForceField{frameToPlot}.network_tracked.node{cellNodeID}.deg;
if initDeg<2
    display('This is a bad node, deg < 2!')
    return;
end

%**************************************************************************
% 3) Find all frames for this cell where the network has not changed      *
%**************************************************************************

% get the list of associated edges:
cellEdgeIDs=constrForceField{frameToPlot}.network_tracked.node{cellNodeID}.edges;
cellEdgeIDs=sort(cellEdgeIDs);

toDoList=[];
for frame=1:length(constrForceField)
    if isfield(constrForceField{frame},'network_tracked') && ... % the frame has been analyzed.
       length(constrForceField{frame}.network_tracked.node)>=cellNodeID && ... % the cell has not been lost.
       ~isempty(constrForceField{frame}.network_tracked.node{cellNodeID}) && ... % the node still exists.
       constrForceField{frame}.network_tracked.node{cellNodeID}.deg==initDeg && ... % the degree = the number of edges is the same.
       sum(ismember(cellEdgeIDs,constrForceField{frame}.network_tracked.node{cellNodeID}.edges))==initDeg; % all edges are identical.
    % only then, it is a good frame:
        toDoList=horzcat(toDoList,frame);
    end
end

%**************************************************************************
% 3) Extract the time course of the forces                                *
%**************************************************************************

X=NaN*zeros(max(toDoList),2+2*length(cellEdgeIDs)); % list of force vectors.
nVec=NaN*zeros(max(toDoList),2+2*length(cellEdgeIDs));
cVec=NaN*zeros(max(toDoList),2+2*length(cellEdgeIDs));
t=NaN*zeros(max(toDoList),1); % list of time points.
col=NaN*zeros(1,2+2*length(cellEdgeIDs)); % description of the columns.
% cols: cell IDs 

% run through all frames that have been analyzed.
for frame=toDoList
    % fill in the residual force for this cell:
    X(frame,1:2)=constrForceField{frame}.network_tracked.node{cellNodeID}.vec;
    col(1:2)=cellNodeID;
    
    edge=constrForceField{frame}.network_tracked.edge;
    node=constrForceField{frame}.network_tracked.node;
    % run through all edges of that particular node:    
    for k=1:length(cellEdgeIDs)
        edgeNum=cellEdgeIDs(k);
        % Check if the force could be determined by the network analysis:
        if ~isempty(~edge{edgeNum}.f1) && ~isempty(edge{edgeNum}.f2) 
            % f1 and f2 are the corresponding forces to .nodes=[n1 n2];
            % calculate the mean vector to get a little bit better value:
            % (In the ideal case, they are the same)
            fnet =1/2*(edge{edgeNum}.f1-edge{edgeNum}.f2);
            ferr =1/2*(edge{edgeNum}.f1+edge{edgeNum}.f2);
            % alternatively, one could take the cluster force here. But
            % then it is difficult to get fc into the right direction!
            
            % Make sure that fnet points in the right direction (outward!).
            % Find the force that belongs to the centered cell:
            [~, loc] = ismember(cellNodeID, edge{edgeNum}.nodes);
            
            if loc==1
                % Then f1 belongs to the center cell, thus fnet has to be
                % flipped (since above fnet=f1-f2):
                fnet=-fnet;
            end
            
            % fill fnet into X. Note that the forces have two components:
            X(frame,2+2*k-1:2+2*k)=fnet;
            X_err(frame,2+2*k-1:2+2*k)=ferr;
            col(2+2*k-1:2+2*k)= setdiff(edge{edgeNum}.nodes, cellNodeID);
            
            % read out the average normal vector of the interface of each
            % force in nVec:
            
            if ~isempty(edge{edgeNum}.nVec_internal)
                currNVec=edge{edgeNum}.nVec_internal;
                nVec(frame,2+2*k-1:2+2*k)=currNVec;
            else
                display(['No normal vector has been calculated for interface: ', num2str(edgeNum),' of frame: ',num2str(frame)]);
            end
            
            % read out the connection vector of the two
            % cells forming the interface, store that in cVec.            
            cNodes=edge{edgeNum}.nodes;
            pos1=node{cNodes(1)}.pos;
            pos2=node{cNodes(2)}.pos;
            vec2to1=(pos2-pos1);
            currCVec=vec2to1/sqrt(vec2to1*vec2to1');
            cVec(frame,2+2*k-1:2+2*k)=currCVec;            
        end
    end    
end
%!!! to be done: check if t is ~ evenly spaced

% determine maximal possible lag:
numPts=sum(~isnan(X(:,1)));
maxLag=floor(numPts/4);

% calculate the correlation coefficients:
[~,cols]=size(X);
c=NaN*zeros(cols,cols,2*maxLag+1);
for i=1:cols
    for j=i:cols %min(i+2,cols):cols
        out1=crossCorr(X(:,i),X(:,j),maxLag); % in the ideal case these entries are all -1!
        c(i,j,:)=out1(:,1);
        c_std(i,j,:)=out1(:,2);
    end
end
figure(1)
plotmatrix(X)

figure(2)
title('The cross correlation for cXY')
for i=1:cols
    for j=i:cols
        subplot(cols,cols,(i-1)*cols+j)
        plot(-maxLag:1:maxLag,reshape(c(i,j,:),[],1))
        ylim([-1 1])
        xlim([-maxLag maxLag])
    end
end

% create c1 and c0:
c1=c;
c0=c;
for i=1:cols
    for j=i:cols
        if mod(i,2)==mod(j,2)
            c0(i,j,:)=NaN;
        else
            c1(i,j,:)=NaN;
        end
    end
end

% figure(3)
% plot(reshape(    c1(5,5,:),[],1))

% create the more condensed matrices c1_sum and c0_sum:
c1_sum=NaN*zeros(cols/2,cols/2,2*maxLag+1);
c0_sum=NaN*zeros(cols/2,cols/2,2*maxLag+1);
for i=1:2:cols-1
    for j=i:2:cols-1
        c1_sum((i+1)/2,(j+1)/2,:)=1/2*(c1(i,j,:)+c1(i+1,j+1,:));        
    end
end

for i=1:2:cols-1
    for j=i:2:cols-1
        c0_sum((i+1)/2,(j+1)/2,:)=1/2*(c0(i+1,j,:)+c0(i,j+1,:));        
    end
end

corrResultsXY.X=X;
corrResultsXY.X_err=X_err;
corrResultsXY.numVal=sum(~isnan(X(:,1)));
corrResultsXY.c1=c1;
corrResultsXY.c0=c0;
corrResultsXY.c1_sum=c1_sum;
corrResultsXY.c0_sum=c0_sum;

%**************************************************************************
% calculate X_I_ for assuming each neighboring cell to be reference:      *
%**************************************************************************
for ref=2:cols/2
    for frame=toDoList
        refnVec=nVec(frame,2*ref-1:2*ref);
        refcVec=cVec(frame,2*ref-1:2*ref);
        for i=1:cols/2
            % Calculate the parallel component:
            corrResultsProjected{ref}.XnVecII(frame,i)=dot(X(frame,2*i-1:2*i),refnVec);
            corrResultsProjected{ref}.XcVecII(frame,i)=dot(X(frame,2*i-1:2*i),refcVec);
            
            % Calculate the perpendicular component:
            corrResultsProjected{ref}.XnVec_I_(frame,i)=dot(X(frame,2*i-1:2*i),[-refnVec(2),refnVec(1)]);
            corrResultsProjected{ref}.XcVec_I_(frame,i)=dot(X(frame,2*i-1:2*i),[-refcVec(2),refcVec(1)]);
        end
    end    
end

% calculate the correlation coefficients:
for ref=2:cols/2
    corrResultsProjected{ref}.cnVecII=NaN*zeros(cols/2,cols/2);
    corrResultsProjected{ref}.ccVecII=NaN*zeros(cols/2,cols/2);    
    corrResultsProjected{ref}.cnVec_I_=NaN*zeros(cols/2,cols/2);
    corrResultsProjected{ref}.ccVec_I_=NaN*zeros(cols/2,cols/2);
end

for ref=2:cols/2
    for i=1:cols/2
        for j=i:cols/2 %min(i+2,cols):cols
            out1=crossCorr(corrResultsProjected{ref}.XnVecII(:,i),corrResultsProjected{ref}.XnVecII(:,j),0); % in the ideal case these entries are all -1!
            corrResultsProjected{ref}.cnVecII(i,j,:)=out1(:,1);
            corrResultsProjected{ref}.cnVecII_std(i,j,:)=out1(:,2);
            
            out1=crossCorr(corrResultsProjected{ref}.XcVecII(:,i),corrResultsProjected{ref}.XcVecII(:,j),0); % in the ideal case these entries are all -1!
            corrResultsProjected{ref}.ccVecII(i,j,:)=out1(:,1);
            corrResultsProjected{ref}.ccVecII_std(i,j,:)=out1(:,2);
            
            out1=crossCorr(corrResultsProjected{ref}.XnVec_I_(:,i),corrResultsProjected{ref}.XnVec_I_(:,j),0); % in the ideal case these entries are all -1!
            corrResultsProjected{ref}.cnVec_I_(i,j,:)=out1(:,1);
            corrResultsProjected{ref}.cnVec_I__std(i,j,:)=out1(:,2);
            
            out1=crossCorr(corrResultsProjected{ref}.XcVec_I_(:,i),corrResultsProjected{ref}.XcVec_I_(:,j),0); % in the ideal case these entries are all -1!
            corrResultsProjected{ref}.ccVec_I_(i,j,:)=out1(:,1);
            corrResultsProjected{ref}.ccVec_I__std(i,j,:)=out1(:,2);
        end
    end
end
% cn and cc seem to be very consistent... as was to be expected!