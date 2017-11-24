function [ candidateMaskNew,linkMask,EPsPostConnect, pixIdxPostConnect,status,TSFigs] = gcaConnectLinearRidges(EPCandidateSort,labelMat,varargin)
% gcaConnectLinearRidges: This function connects end-points of linear
% ridge candidates by searching for linkage candidates using a KD tree
% and then performing maxWeigthedGraphMatching to determine the optimal
% attachment in the case there are competing possible attachments
% - the cost is currenlty a simple linear summation of the
% local co-linearity of the tangent at the endpoint of the ridge and the vector formed
% by the connection between the two points. Negative values for these
% colinearity terms are not allowed for linkage distances >  and these links will be removed as
% these linkages will result in angle formation > 90 degrees between the
% two ridges. In this function we also disfavor formaing junctions by removing
% connections that pass over the original candidate mask.
%
%% INPUT:
%
%   EPCandidateSort: (REQUIRED) : 1 x c cell-array of 2x4 double arrays
%         where c is the number of candidate ridges.
%         and each cell contains a double array with the
%         endpoint coordinate (x,y) where x is the local direction
%         moving toward the endpoint along the ridge(x,y)
%         sorted by the label
%
%   labelMat: (REQUIRED) : rxc unit8 array
%         where r (row) is the height (ny) and c (col) is the width
%         (nx) of the original input image
%         the candiate ridges labeled by connected
%         component. Labels should correspond to the order of
%         EPCandidateSort-
%         Output of labelmatrix.m function
%   (MB CHECK BEFORE RELEASE: might want to make one input of a candidateMask and calculate
%   these other two inputs here) -
%
% OPTIONAL:
%   img: (OPTIONAL)
%
% PARAMS:
%
%  'MaxRadiusLink' (PARAM): Positive Scalar
%        Maximum radius for linking the endpoints neighboring large-scale ridge candidates
%       Default : 10 (In Pixels)
%       can make 5 for first step of the linking.
%
%  'NoLinkDistanceFromBorder' (PARAM) : Positive Scalar
%        Endpoints within this distance of the image boundary will not be
%        considered for ridge linking.  This parameter was added to ensure that
%        no ridge filter edge effects are linked. Practically
%       Default : 10 (In Pixels)
%       can make 0 for first step of linking
%  'MaxRadiusNoGeoTerm'  (PARAM) : Positive Scalar
%       The max radius for which geometry will be considered when linking
%       ridges - for the large scale ridge linking one often does not want
%       to consider geometry of the ridges as the current scale integration often
%       leads to small and abrupt changes in the NMS when the neurite
%       changes scales and we wish to link these.
%      Default : 3 (In Pixels)
%      can make 0 for first
%  'GeoThresh' (PARAM) : Positive Scalar
%
%
%  Default : 0 for large scale linking
%            0.95 for small scale first step linking)
% OUTPUT:
%   candidateMaskNew: rxc logical array
%         where r (row) is the height (ny) and c (col) is the width
%        (nx) of the original input image of the candidates post linkage
%
%   linkMask: rxc logical array
%         where r (row) is the height (ny) and c (col) is the width
%         (nx) of the original input image marking the optimized interpolated
%         new linkages between candidate ridges
%   
%   pixIdxPostConnect : 1xc cell of rx1 doubles
%         storing the pixel indices for each new ridge candidate
%         c is the number of candidate labels (should be lower than the
%         original number of labels if connections were made) and provides
%         the ridge candidates ID. 
%         r is the number of pixels corresponding to each candidate after
%         ridge candidate building.
%         Note that in this format two pixels can share multiple labels 
%         and hence crossovers are resolved. (This is not the case if we
%         put these in a simple labelMat) 
%   
% status: 1 if viable links were found: 0 if no links were made
%% OLD 
% %%%  EPCandidateSortPostConnect: 1 x c cell-array of 2x2 double arrays
%         where c is still the number of the original connected component
%         ridge labels prior to linking and each cell contains a double array with the
%         endpoint coordinate (x,y) (vectors removed). Labels which have been removed due to
%         linking (merging) are replace by NaNs and the new endpoints of
%         each ridge label are documented
%         (MB CHECK TO MAKE SURE THIS IS OK
%         BEFORE RELEASE) - CURRENTLY DO NOT USE this output - more
%         important when consider the filopodia
%
%%%%%  labelMatPostConnect:  rxc int8 array
%         where r (row) is the height (ny) and c (col) is the width
%         (nx) of the original input image. The label matrix merges the
%         labels of two linked CCs. (MB CHECK TO MAKE SURE THIS IS OK
%         BEFORE RELEASE: CURRENTLY DO NOT USE this output)
%

%% InputParser
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;
%REQUIRED
ip.addRequired('EPCandidateSort');
ip.addRequired('labelMat');
%OPTIONAL
ip.addOptional('img',[]);

% PARAMETERS
ip.addParameter('MaxRadiusLink',10);
ip.addParameter('NoLinkDistanceFromBorder',0);
ip.addParameter('MaxRadiusNoGeoTerm',0) ;
ip.addParameter('GeoThresh',0.9);
ip.addParameter('TSOverlays',true);

ip.parse(EPCandidateSort,labelMat,varargin{:});

%% Initiate
candidateMask = labelMat>0;
[ny,nx] = size(labelMat);
imSize = [ny,nx];
countFigs = 1;
EPsPostConnect = [];
pixIdxPostConnect = []; 
TSFigs = []; 
linkMask = zeros([ny,nx]);
status = 0; 
%%
endPoints = vertcat(EPCandidateSort{:}); % taking these out of a cell array so
endPoints =  endPoints(:,1:2); % take first two columns as added vector 20140913
% create a repmat for indexing
%labels = arrayfun(@(i) repmat(i,2,1),1:numel(EPCandidateSort),'uniformoutput',0);
%labels = vertcat(labels{:});

% typically you want to exclude endpoints from being connected that fall
% within a certain range of the image boundary as these pieces are more
% likely connected to the boundary itself
boundaryMask = zeros(imSize);
boundaryMask(1:imSize(1),1) =1;
boundaryMask(1:imSize(1),imSize(2))=1;
boundaryMask(1,1:imSize(2))= 1;
boundaryMask(imSize(1),1:imSize(2)) =1;
[ boundaryY,boundaryX]= ind2sub(imSize,find(boundaryMask==1));
boundaryCoordXY = [boundaryX,boundaryY];
% first input is the input points second input is the query points
% want to find the endpoints around the boundary
[idxEPsNearBound,~] = KDTreeBallQuery(endPoints,boundaryCoordXY,ip.Results.NoLinkDistanceFromBorder);
toRemove = unique(vertcat(idxEPsNearBound{:}));
%% TSOverlay
endPoints(toRemove,:) = [];
sanityCheck = 0;
if sanityCheck == 1
    figure;
    imshow(candidateMask,[])
    hold on
    scatter(endPoints(:,1),endPoints(:,2),'g','filled');
end
% Added 20150812 
if ~isempty(endPoints); 
    
%% Find all enpoints within x radius of one another- note there will be some
% redundancy that one has to filter and each query point will find itself.
[idx,d] = KDTreeBallQuery(endPoints, endPoints, ip.Results.MaxRadiusLink); % originally 5

%% Sanity Check1: Plot the query and surrounding points found - note it will
% always find 'itself' and this will be marked by a zero distance.
% sanityCheck = 0;
% if sanityCheck == 1
%     nQueries = numel(idx);
%     % show results of search
%     for iQuery =1:nQueries
%         figure;
%         imshow(labelMat>0,[]);
%         hold on
%         % plot all endpoints in yellow
%         scatter(endPoints(:,1),endPoints(:,2),'filled','g');
%         % plot query as red point
%         scatter(endPoints(iQuery,1),endPoints(iQuery,2),'filled','y');
%         idxFoundPerQuery = vertcat(idx{iQuery}); % indices of endPoints
%         dFoundPerQuery = vertcat(d{iQuery});
%         % plot found points within radius of query
%         scatter(endPoints(idxFoundPerQuery,1),endPoints(idxFoundPerQuery,2),'filled','r');
%         % plot the distances calculated for all found points.
%         arrayfun(@(i) text(endPoints(idxFoundPerQuery(i),1),endPoints(idxFoundPerQuery(i),2),num2str(dFoundPerQuery(i)),'color','y'),1:length(idxFoundPerQuery));
%         saveas(gcf,[num2str(iQuery,'%03d') '.tif']);
%         saveas(gcf,[num2str(iQuery,'%03d') '.fig']);
%     end
% end

% format edges: column 1 is indices of endpoint1 (Vertex1) and column 1 is
% the indice of endpoint2
E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:size(endPoints,1), 'UniformOutput', false);
E = vertcat(E{:});
d = vertcat(d{:});
EOrig = E;
% Remove Redundancy
d(E(:,1)<=E(:,2))=[]; % 20140826 Note here caught a small bug don't think I was
% filtering appropirately and this would be reflected ultimately in the
% weights - filter d before truncating
E(E(:,1)<=E(:,2),:)=[];
% remove self associations (ie connection with same label) : this should
% take out those points where it found 'itself' and the distance = 0
% as well as avoid possible edges between endpoints of the same piece.
label1 = arrayfun(@(i) labelMat(sub2ind([ny,nx],endPoints(E(i,1),2), endPoints(E(i,1),1))),1:length(E(:,1)))  ;
label2 = arrayfun(@(i) labelMat(sub2ind([ny,nx],endPoints(E(i,2),2),endPoints(E(i,2),1))),1:length(E(:,1)));

E = E(label1~=label2,:);
d = d(label1~=label2,:);

%% Remove all overlap with previous pixels (this might not be the best option for filo cands might just need to reformat)
% find path between edges
paths=arrayfun(@(i) gcaBresenham([endPoints(E(i,1),1) endPoints(E(i,1),2)], [endPoints(E(i,2),1) endPoints(E(i,2),2)]),...
    1:length(E(:,1)),'uniformoutput',0);

pathsidx = cellfun(@(x) sub2ind([ny,nx],x(2:end-1,2),x(2:end-1,1)),paths,'uniformoutput',0);
% check for overlap with candidateMask
candPix = find(candidateMask==1);
noOverlap = cellfun(@(x) isempty(intersect(x,candPix)),pathsidx); % have no interesction with the candidate mask
E = E(noOverlap,:);
d = d(noOverlap);

%% TSOVerlays : Plot Connections
% (ie after self association,redundancy filter, and overlap filter)

if ip.Results.TSOverlays == true;
    TSFigs(countFigs).h = setFigure(nx,ny,'on');
    TSFigs(countFigs).name =  'Plot Connections';
    TSFigs(countFigs).group = 'Connect_Ridge_Ends'; 
        
    if ~isempty(ip.Results.img);
        imshow(-ip.Results.img,[]);
        hold on
    end
    
    spy(labelMat>0,'b',2);
    hold on
    % should eventually color code by cost... but for now
    %cmap = jet(length(E));
    % easier to keep track of indexing if just use a for loop
    
    idxAll1 = E(:,1);
    idxAll2= E(:,2);
    endPoints1KD = endPoints(idxAll1,:);
    endPoints2KD = endPoints(idxAll2,:);
    
    % scatter the endpoints found
    nEdges = size(E(:,1));
    scatter(endPoints1KD(:,1),endPoints1KD(:,2),5,'b','filled');
    scatter(endPoints2KD(:,1),endPoints2KD(:,2),5,'b','filled');
    
    arrayfun(@(x) plot([endPoints1KD(x,1),endPoints2KD(x,1)],...
        [endPoints1KD(x,2),endPoints2KD(x,2)],'color','g','Linewidth',1),1:nEdges);
    
    countFigs = countFigs +1;
    
    
end % ip.Results

%% get other costs of these edges
% DELTA SCALE - would need to input the scale map
if ~isempty(E)
    % GEOMETRIC -
    % get the local vectors of the tangent at the endpoint
    vect = vertcat(EPCandidateSort{:});
    
    vect = vect(:,3:4); % all endpoints
    vect(toRemove,:) = [];
    
    % vector of connectivity from edge 2 to 1
    deltXCon21 = arrayfun(@(i) endPoints(E(i,1),1)-endPoints(E(i,2),1),1:length(E(:,1)));
    deltYCon21= arrayfun(@(i) endPoints(E(i,1),2)-endPoints(E(i,2),2),1:length(E(:,1)));
    
    % vector of connectivity 1 to 2
    deltXCon12 = arrayfun(@(i) endPoints(E(i,2),1)-endPoints(E(i,1),1),1:length(E(:,1)));
    deltYCon12 = arrayfun(@(i) endPoints(E(i,2),2)-endPoints(E(i,1),2),1:length(E(:,1)));
    
    % for each edge calculated the dot product of the tangent at end point and
    % the vector between the two endpoints.
    
    dotProd12 = arrayfun(@(i) dot([vect(E(i,1),1) vect(E(i,1),2)],[deltXCon12(i) deltYCon12(i)])/d(i),1:length(E(:,1)));
    dotProd21 = arrayfun(@(i) dot(vect(E(i,2),:),[deltXCon21(i) deltYCon21(i)])/d(i),1:length(E(:,1)));
    %%
    if ip.Results.TSOverlays == true;
        TSFigs(countFigs).h = setFigure(nx,ny,'on');
        TSFigs(countFigs).name =  'Plot_Connections_With_Vectors';
        TSFigs(countFigs).group = 'Connect_Ridge_Ends'; 
        
        if ~isempty(ip.Results.img);
            imshow(-ip.Results.img,[]);
            hold on
        end
        
        spy(labelMat>0,'b',1);
        hold on
        % should eventually color code by cost... but for now
        %cmap = jet(length(E));
        % easier to keep track of indexing if just use a for loop
        
        idxAll1 = E(:,1);
        idxAll2= E(:,2);
        endPoints1KD = endPoints(idxAll1,:);
        endPoints2KD = endPoints(idxAll2,:);
        
        % scatter the endpoints found
        nEdges = size(E(:,1),1);
        scatter(endPoints1KD(:,1),endPoints1KD(:,2),5,'b','filled');
        scatter(endPoints2KD(:,1),endPoints2KD(:,2),5,'b','filled');
        
        arrayfun(@(x) plot([endPoints1KD(x,1),endPoints2KD(x,1)],...
            [endPoints1KD(x,2),endPoints2KD(x,2)],'color','g'),1:nEdges);
        
        quiver(endPoints1KD(:,1),endPoints1KD(:,2),vect(idxAll1,1),vect(idxAll1,2),6,'color','b');
        quiver(endPoints2KD(:,1),endPoints2KD(:,2),vect(idxAll2,1),vect(idxAll2,2),6,'color','b');
        deltCon12N = [deltXCon12'./d deltYCon12'./d];
        deltCon21N = [deltXCon21'./d deltYCon21'./d];
        
        
        quiver(endPoints1KD(:,1),endPoints1KD(:,2),deltCon12N(:,1),deltCon12N(:,2),6,'color','g');
        quiver(endPoints2KD(:,1),endPoints2KD(:,2),deltCon21N(:,1),deltCon21N(:,2),6,'color','g');
        
        countFigs = countFigs +1;
        
        
    end % ip.Results
        
    %% Filter for long distance links by geometry
    idxFilt = (d>ip.Results.MaxRadiusNoGeoTerm & (dotProd21' <=ip.Results.GeoThresh | dotProd12' <=ip.Results.GeoThresh)); % could potentially make the linearity threshold
    ERemove = E(idxFilt,:);
    
    %     %% Plot the bad link sanity check if user desires
    %     if sanityCheck2 == 1
    %         hold on
    %         % easier to keep track of indexing if just use a for loop
    %         for iEdge = 1:length(ERemove);
    %             %    % coords of first point of edge
    %             xCoord1 = endPoints(ERemove(iEdge,1),1);
    %             yCoord1 = endPoints(ERemove(iEdge,1),2);
    %             %
    %             %    % coords of second point of edge
    %             xCoord2 = endPoints(ERemove(iEdge,2),1);
    %             yCoord2 = endPoints(ERemove(iEdge,2),2);
    %             %
    %             %
    %             %
    %             hold on
    %             %
    %             %   plot([xCoord1,xCoord2],[yCoord1,yCoord2],'color',cmap(iEdge,:));
    %             plot([xCoord1,xCoord2],[yCoord1,yCoord2],'color','r','Linewidth',2);
    %         end
    %         %plot the remove vectors in red
    %         idx21Bad = find(d>ip.Results.MaxRadiusNoGeoTerm & (dotProd21' <=0));
    %         idx12Bad = find(d>ip.Results.MaxRadiusNoGeoTerm & (dotProd12'<=0));
    %         if (~isempty(idx21Bad) || ~isempty(idx12Bad))
    %             quiver(endPoints(E(idx21Bad,2),1),endPoints(E(idx21Bad,2),2),...
    %                 vect(E(idx21Bad,2),1), vect(E(idx21Bad,2),2),'color','g');
    %
    %             quiver(endPoints(E(idx12Bad,1),1),endPoints(E(idx12Bad,1),2),...
    %                 vect(E(idx12Bad,1),1), vect(E(idx12Bad,1),2),'color','m');
    %
    %
    %             %strDot12 = arrayfun(@(x) num2str(dotProd12,2),1:length(dotProd12),'uniformoutput',0);
    %
    %
    %             arrayfun(@(i) text(endPoints(E(idx12Bad(i),1),1),...
    %                 endPoints(E(idx12Bad(i),1),2),num2str(dotProd12(idx12Bad(i)),2),...
    %                 'FontSize',12,'Color','g'),1:length(idx12Bad));
    %             % plot the distance at the query point
    %             %  %  text(xCoord1,yCoord1,num2str(d(iEdge),3),'color','y');
    %             %arrayfun(@(i) text(endPoints(E(idx21Bad(i),1),1),...
    %             arrayfun(@(i) text(endPoints(E(idx21Bad(i),2),1),...
    %                 endPoints(E(idx21Bad(i),2),2),num2str(dotProd21(idx21Bad(i)),2),...
    %                 'FontSize',12,'Color','m'),1:length(idx21Bad));
    %             %
    %         end
    %     end
    %%
    %
    d = max(d)-d; % max with =0
    d = d./max(d);
    costTotalPreFilt = d+dotProd12'+dotProd21';
    EPreFilt = E;
    E(idxFilt,:) = [];
    d(idxFilt) = [];
    dotProd21(idxFilt) = [];
    dotProd12(idxFilt) = [];
    


%% TSOverlays : Plot Cost by Color
if ip.Results.TSOverlays == true;
    
    TSFigs(countFigs).h = setFigure(nx,ny,'on');
    TSFigs(countFigs).name = ' Color_By_Cost';
    TSFigs(countFigs).group = 'Connect_Ridge_Ends'; 
        
    if ~isempty(ip.Results.img);
        imshow(-ip.Results.img,[]);
        hold on
    end
    
    scatter(endPoints1KD(:,1),endPoints1KD(:,2),5,'b','filled');
    scatter(endPoints2KD(:,1),endPoints2KD(:,2),5,'b','filled');
    
    cMapLength=128; cMap=jet(cMapLength);
    mapper=linspace(min(costTotalPreFilt),max(costTotalPreFilt),cMapLength)';
    
    % get closest colormap index for each feature
    D=createDistanceMatrix(costTotalPreFilt,mapper);
    [sD,idxCMap]=sort(abs(D),2);
    %  %
    %  % for each Edge find the xy coords of the seed point and the candidate
    %  % point and plot
    %          for i = 1:length(E(:,1))
    %              % get the indexes of the edge relative to the original input EPs
    %             idxCand = E(i,2);
    %             idxSeed = E(i,1);
    %             plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)]);
    %
    %         end
    for k=1:cMapLength
        idxCand = EPreFilt(idxCMap(:,1) == k,2);
        idxSeed = EPreFilt(idxCMap(:,1)==k,1);
        for iEdge = 1:length(idxCand) % some can have the same color
            plot([endPoints(idxCand(iEdge),1),endPoints(idxSeed(iEdge),1)],...
                [endPoints(idxCand(iEdge),2),endPoints(idxSeed(iEdge),2)],'color',cMap(k,:),'Linewidth',2);
        end
    end
    spy(labelMat>0,'b');
    text(5,5,'Color Paths By Cost ','FontSize',10,'Color','k');
    text(5,15,'Red : High : Stong Path' ,'FontSize',10,'Color','r');
    text(5,25,'Blue : Low : Weak Path', 'FontSize',10,'Color','b');
    text(5,35,['MaxRadius = ' num2str(ip.Results.MaxRadiusLink) ' Pixels']);
    countFigs = countFigs+1;
    
    %% TSOverlays : Plot Cost by Color
    if ip.Results.TSOverlays == true;
        
        TSFigs(countFigs).h = setFigure(nx,ny,'on');
        TSFigs(countFigs).name = 'KD_Results_AFter_Filter_By_Geometry';
        TSFigs(countFigs).group = 'Connect_Ridge_Ends'; 
        
        
        if ~isempty(ip.Results.img);
            imshow(-ip.Results.img,[]);
            hold on
        end
        
        
        scatter(endPoints1KD(:,1),endPoints1KD(:,2),5,'k','filled');
        scatter(endPoints2KD(:,1),endPoints2KD(:,2),5,'k','filled');
        
        
        idxCMap(idxFilt,:) =[];
        for k=1:cMapLength
            
            idxCand = E(idxCMap(:,1) == k,2);
            idxSeed = E(idxCMap(:,1)==k,1);
            for iEdge = 1:length(idxCand) % some can have the same color
                plot([endPoints(idxCand(iEdge),1),endPoints(idxSeed(iEdge),1)],...
                    [endPoints(idxCand(iEdge),2),endPoints(idxSeed(iEdge),2)],'color',cMap(k,:),'Linewidth',2);
            end
        end
        spy(labelMat>0,'k');
        text(5,5,'After Application of Geometry Threshold','FontSize',10);
        text(5,20,['of ' num2str(ip.Results.GeoThresh)],'FontSize',10);
        
        countFigs = countFigs+1;
        
    end % ip.Results.
    
      
end % ip.Results.
end % isempty E
%%
if ~isempty(E) %check if there are reasonable edges
    % make sure to make d so that minimum d are favored
    
    costTotal = d + dotProd12' + dotProd21'; % again this cost need to be refined a bit
    
    % 20140917 not sure why I need nodelabels here... need to check to
    % see this might have been something specific to the filo...
    [numberNodes,~,nodeLabels] =unique(E(:),'stable');
    numberNodes = length(numberNodes);
    EFinal(:,1) = nodeLabels(1:length(E(:,1)));
    EFinal(:,2) = nodeLabels(length(E(:,1))+1:end);
    
    
    M = maxWeightedMatching(numberNodes, EFinal, costTotal);
    E = E(M,:);
    
    paths=arrayfun(@(i) gcaBresenham([endPoints(E(i,1),1) endPoints(E(i,1),2)], [endPoints(E(i,2),1) endPoints(E(i,2),2)]),...
        1:length(E(:,1)),'uniformoutput',0);
    linkMask = zeros([ny,nx]);
    
    links = vertcat(paths{:});
    
else % no good links
    
    links = []; % since links is empty this will default to mask bad
    linkMask = zeros([ny,nx]);
end % isempty E



if ~isempty(links) % nothing that falls under this criteria
    % Add links to candidate mask
    idxLinks = sub2ind([ny,nx],links(:,2),links(:,1));
    % need to add to mask  iteratively and check for junctions... later
    % see how much comp time this sucks.
    %         linkMaskCheck = zeros(size(maxTh));
    %         for iLink = 1:length(links(:,2))
    %         linkMaskIndCheck(links(iLink,2),links(iLink,1))=1;
    %
    %         end
       
    linkMask(idxLinks) = 1;
    candidateMaskNew = (candidateMask| linkMask);
    % Fix the label mat to unite those pixels that need to be clustered
    % labelMatPostConnect = labelMat; % initiate a labelMatrix that connects the two pieces
    
    % pixIdxLabelsPostConnect = arrayfun(@(x) find(labelMat== x),1:nLabels,'uniformoutput',0);
    labelsAll = labelMat(labelMat~=0);
    labelsAll = unique(labelsAll);
    % for all edges
    
    
    % get the labels of endpoints
    
    % get the label of endpoint 1
    labels1 = labelMat(sub2ind([ny,nx],endPoints(E(:,1),2), endPoints(E(:,1),1)));
    % get the label of endpoint 2 in linkage
    labels2 = labelMat(sub2ind([ny,nx],endPoints(E(:,2),2),endPoints(E(:,2),1)));
    
    %group edges with commonn labesl
    labelsMoreThanOneEdge = intersect(labels1, labels2);
    if ~isempty(labelsMoreThanOneEdge)
        for iGroup = 1:length(labelsMoreThanOneEdge);
            % find idx of edges that that need to be combined
            idxGrpEdges = find(labels1==labelsMoreThanOneEdge(iGroup) | labels2==labelsMoreThanOneEdge(iGroup));
            idxGrpEdgesAll{iGroup} = labels1==labelsMoreThanOneEdge(iGroup) | labels2==labelsMoreThanOneEdge(iGroup); % marks Edges
            % get paths for these edges
            pathsConnect = vertcat( paths{idxGrpEdges});
            pixIdxPaths = sub2ind([ny,nx], pathsConnect(:,2),pathsConnect(:,1));
            pixIdxCands= arrayfun(@(x) find(labelMat==labels1(idxGrpEdges(x))|labelMat==labels2(idxGrpEdges(x))),1:length(idxGrpEdges), ...
                'uniformoutput',0);
            pixIdxNew = [pixIdxPaths; vertcat(pixIdxCands{:})];
            pixIdxPostConnect{iGroup} = pixIdxNew;
            EPsPostConnect{iGroup} =  getEndpoints(pixIdxNew,[ny,nx], 0,1);
            clear pixIdxNew
        end % for iGroup
        
        
        idxNotEdgeGrp = ~sum(horzcat(idxGrpEdgesAll{:}),2);
        nPiecesGroup = numel(pixIdxPostConnect);
    else
        idxNotEdgeGrp = ones(size(E,1),1);
        nPiecesGroup = 0;
    end %isempty % do not have to cluster more than one edge
    %pixIdxPostConnect{iGroup} =
    
    
    %% START FIXING
    % the idx of the edges that need not be grouped 
    notGrpEdges  =  find(idxNotEdgeGrp);
    labels1NG = labels2(notGrpEdges);
    labels2NG = labels1(notGrpEdges);
    pathsNG = paths(notGrpEdges);
     
    % document pieces that form a single connection
    for i = 1:length(labels1NG)
        newPiece{1} = find(labelMat==labels1NG(i));
        newPiece{3} = find(labelMat==labels2NG(i));
        newPiece{2} = sub2ind([ny,nx],pathsNG{i}(:,2),pathsNG{i}(:,1));
        idxC = i + nPiecesGroup;
        pixIdxPostConnect{idxC} = vertcat(newPiece{:});
        clear newPiece;
        EPsPostConnect{idxC} = getEndpoints(pixIdxPostConnect{idxC},[ny,nx],0,1);
    end
    %% get the other labels
    nPiecesConnect = numel(EPsPostConnect);
    % find labels not connected.
    labelsUsed = [labels1;labels2];
    nonConnect = setdiff(labelsAll,labelsUsed);
    % labelsAll
    for iCand = 1:length(nonConnect)
        idxC = nPiecesConnect+iCand;
        pixIdxPostConnect{idxC} = find(labelMat==nonConnect(iCand));
        EPsPostConnect{idxC} = getEndpoints(pixIdxPostConnect{idxC},[ny,nx],0,1); 
    end
    
    % remove the extra pixels per pixIdx 
    pixIdxPostConnect = cellfun(@(x) unique(x),pixIdxPostConnect,'uniformoutput',0); 
    
    %% TSOverlays Figure 
    if ip.Results.TSOverlays == true;
        
        TSFigs(countFigs).h = setFigure(nx,ny,'on');
        TSFigs(countFigs).name = 'Post_Connection_Labels';
        TSFigs(countFigs).group = 'Connect_Ridge_Ends'; 
     
        
        imshow(labelMat>0,[]);
        hold on
        % plot the new EPS
        % plot the new labels
        cmap = lines(numel(EPsPostConnect));
        
        nLabels =  numel(EPsPostConnect);
        
        %cmap = lines(7);
        % sanity check
        [yCand,xCand] = cellfun(@(x) ind2sub([ny,nx],x),pixIdxPostConnect,'uniformoutput',0);
        
        for x = 1:nLabels
            scatter(EPsPostConnect{x}(:,1),EPsPostConnect{x}(:,2),50,cmap(x,:),'filled');
            scatter(xCand{x},yCand{x},20,cmap(x,:),'filled');
            
            % scatter(
        end ;
        % plot shared pixels
        idxAll =  vertcat(pixIdxPostConnect{:});
        [~,idxSingle] = unique(idxAll,'stable');
        shared = idxAll;
        
        shared(idxSingle) = [];
        
        [yShare,xShare] = ind2sub([ny,nx],shared);
        scatter(xShare,yShare,10,'w');
        countFigs = countFigs+1;
    end
     
    status = 1; % there were links
else
    candidateMaskNew = candidateMask;
    %labelMatPostConnect = labelMat;
    
end
else % added 20150812 .. see if can consolidate the if statmements quick fix for now.  
    candidateMaskNew = candidateMask; % for now see we can just keep these candidates what happens - maybe in the future remove hte border cands completely 
    %labelMatPostConnect = labelMat;
    
end 

