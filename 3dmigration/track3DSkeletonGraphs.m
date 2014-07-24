function tracks = track3DSkeletonGraphs(skelIn,varargin)




showPlots = false;%SDhow plots for debug/dev


ip = inputParser;
ip.FunctionName = mfilename;
ip.addRequired('skelIn',@(x)(isstruct(x) && numel(x) > 1));
ip.addParamValue('MaxDisp',Inf,@(x)(numel(x)==1 && x > 0));
ip.parse(skelIn,varargin{:});
p = ip.Results;

nFrames = numel(skelIn);

%First we go through and get some skeleton properties to help initialize
%our tracking
iTips = cell(nFrames,1);
nVerts = zeros(nFrames,1);
nTips = zeros(nFrames,1);

for j = 1:nFrames
    
    nVerts(j) = size(skelIn(j).vertices,1);
    
    %Find just the tips in this skeleton
    iTips{j} = findTips(skelIn(j).edges,nVerts(j));
    nTips(j) = numel(iTips{j});
        
end

nTipTotal = sum(nTips);
nTipCumSum = vertcat(0,cumsum(nTips));

%Initialize the temporal vertices matrix - these are each point in a track
tVerts = zeros(nTipTotal,2);
%Over-initialize the temporal edges matrix - these are the links between
%points
tEdges = zeros(nTipTotal,2);
nEdgesCumSum = zeros(nFrames,1);
iEdgesCurr = cell(nFrames,1);
iCurrVerts = cell(nFrames,1);
for j = 1:nFrames
    
    iCurrVerts{j} = nTipCumSum(j)+1:nTipCumSum(j+1);
        
    tVerts(iCurrVerts{j},1) = j;
    tVerts(iCurrVerts{j},2) = iTips{j};
    
    
    if j > 1
        
       %Make the distance matrix for t and t-1
       distMat = createDistanceMatrix(skelIn(j-1).vertices(iTips{j-1},:),skelIn(j).vertices(iTips{j},:));       
       %Remove potential links which violate the maximum displacement
       %parameter
       distMat(distMat > p.MaxDisp) = -1;
       
       %Call the damn LAP on it. Use defaults for now.
       [iAB,iBA] = lap(distMat,[],[],1);
       %Get rid of the extra ones from the cost matrix construction so I
       %don't get confused. Fuck you I'm tired.
       iAB = iAB(1:nTips(j-1));       
       isLinked = iAB <= nTips(j);
       nEdgesCurr = nnz(isLinked);
       nEdgesCumSum(j) = nEdgesCumSum(j-1)+nEdgesCurr;
       iEdgesCurr{j-1} = nEdgesCumSum(j-1)+1:nEdgesCumSum(j);       
       tEdges(iEdgesCurr{j-1},1) = iCurrVerts{j-1}(isLinked);
       tEdges(iEdgesCurr{j-1},2) = iCurrVerts{j}(iAB(isLinked));                            
       
       if showPlots
            figure
            hold on
            plotSkel(skelIn(j-1))
            plotSkel(skelIn(j));
            for k = 1:nEdgesCurr
%                plot3([skelIn(j-1).vertices(iTips{j-1}(k),2) skelIn(j).vertices(iTips{j}(iAB(k)),2)],...
%                      [skelIn(j-1).vertices(iTips{j-1}(k),1) skelIn(j).vertices(iTips{j}(iAB(k)),1)],...
%                      [skelIn(j-1).vertices(iTips{j-1}(k),3) skelIn(j).vertices(iTips{j}(iAB(k)),3)],'--k')
            %Get all these indices out to make sure I didn't fuck it up
            %using my incredibly confusing data structure
            currEdge = iEdgesCurr{j-1}(k);%Index of current edge
            currT1 = tVerts(tEdges(currEdge,1),1);%Current t edge start time
            currT2 = tVerts(tEdges(currEdge,2),1);%Current t edge end time
            currV1 = tVerts(tEdges(currEdge,1),2);%Current t edge start skeleton vertex index
            currV2 = tVerts(tEdges(currEdge,2),2);%Current t edge end skeleton vertex index
            
            plot3([skelIn(currT1).vertices(currV1,2) skelIn(currT2).vertices(currV2,2)],...
                 [skelIn(currT1).vertices(currV1,1) skelIn(currT2).vertices(currV2,1)],...
                 [skelIn(currT1).vertices(currV1,3) skelIn(currT2).vertices(currV2,3)],'--k'); 
           end
       end
    end    
end

%Resize the temporal edges array which we over-initialized.
nTEdgesTot = nEdgesCumSum(end);
tEdges = tEdges(1:nTEdgesTot,:);
    
%Make the graph adjacency matrix, find connected components - these are our tracks.
tEdgeMat = sparse(tEdges(:,1),tEdges(:,2),1,nTipTotal,nTipTotal,nTEdgesTot);
[nTracks,trackIndices] = graphconncomp(tEdgeMat,'Weak',true);

%Combine this into output structure
tracks.tEdges = tEdges;
tracks.tVerts = tVerts;
tracks.trackIndices = trackIndices;
tracks.nTracks = nTracks;
tracks.nFrames = nFrames;
tracks = convertTrackStruct(tracks,skelIn,'Imaris');


% if showPlots
%     close all
%     figure
%     hold on
%     %Convert for plotting
%     
%     frameCols = jet(nFrames);
%     for j = 1:nFrames
%         plotSkel(skelIn(j),[],[],[],frameCols(j,:));
%     end
%     
% end




    