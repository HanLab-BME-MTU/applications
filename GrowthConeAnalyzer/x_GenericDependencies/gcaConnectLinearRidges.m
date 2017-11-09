function [ candidateMaskNew,linkMask,EPCandidateSort, labelMatPostConnect,status] = gcaConnectLinearRidges(EPCandidateSort,labelMat,varargin)
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
% PARAMS:         
%
%  'MaxRadiusLargeScaleLink' (PARAM): Positive Scalar
%        Maximum radius for linking the endpoints neighboring large-scale ridge candidates
%       Default : 10 (In Pixels)          
%
%  'NoLinkDistanceFromBorder' (PARAM) : Positive Scalar
%        Endpoints within this distance of the image boundary will not be
%        considered for ridge linking.  This parameter was added to ensure that 
%        no ridge filter edge effects are linked. Practically 
%       Default : 10 (In Pixels)
%      
%  'MaxRadiusNoGeoTerm'  (PARAM) : Positive Scalar
%       The max radius for which geometry will be considered when linking 
%       ridges - for the large scale ridge linking one often does not want
%       to consider geometry of the ridges as the current scale integration often 
%       leads to small and abrupt changes in the NMS when the neurite
%       changes scales and we wish to link these. 
%      Default : 3 (In Pixels)
%
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
%   EPCandidateSortPostConnect: 1 x c cell-array of 2x2 double arrays
%         where c is still the number of the original connected component 
%         ridge labels prior to linking and each cell contains a double array with the 
%         endpoint coordinate (x,y) (vectors removed). Labels which have been removed due to
%         linking (merging) are replace by NaNs and the new endpoints of
%         each ridge label are documented 
%         (MB CHECK TO MAKE SURE THIS IS OK
%         BEFORE RELEASE) - CURRENTLY DO NOT USE this output - more
%         important when consider the filopodia 
% 
%  labelMatPostConnect:  rxc int8 array 
%         where r (row) is the height (ny) and c (col) is the width
%         (nx) of the original input image. The label matrix merges the
%         labels of two linked CCs. (MB CHECK TO MAKE SURE THIS IS OK
%         BEFORE RELEASE: CURRENTLY DO NOT USE this output) 
% 
%  status: 1 if viable links were found: 0 if no links were made 
%% InputParser 
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true; 
%REQUIRED
ip.addRequired('EPCandidateSort');
ip.addRequired('labelMat'); 

% PARAMETERS
ip.addParameter('MaxRadiusLargeScaleLink',10); 
ip.addParameter('NoLinkDistanceFromBorder',10); 
ip.addParameter('MaxRadiusNoGeoTerm',3) ; 
%ip.addParameter('TSOverlays',true,@(x) islogical(x));

ip.parse(EPCandidateSort,labelMat,varargin{:});



%% Initiate 
candidateMask = labelMat>0; 
[ny,nx] = size(labelMat);
imSize = [ny,nx];

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

endPoints(toRemove,:) = [];
sanityCheck = 0;
if sanityCheck == 1
    figure;
    imshow(candidateMask,[])
    hold on
    scatter(endPoints(:,1),endPoints(:,2),'g','filled');
end

%% Find all enpoints within x radius of one another- note there will be some
% redundancy that one has to filter and each query point will find itself.
[idx,d] = KDTreeBallQuery(endPoints, endPoints, ip.Results.MaxRadiusLargeScaleLink); % originally 5



%% Sanity Check1: Plot the query and surrounding points found - note it will
% always find 'itself' and this will be marked by a zero distance.
sanityCheck = 0;
if sanityCheck == 1
    nQueries = numel(idx);
    % show results of search
    for iQuery =1:nQueries
        figure;
        imshow(labelMat>0,[]);
        hold on
        % plot all endpoints in yellow
        scatter(endPoints(:,1),endPoints(:,2),'filled','g');
        % plot query as red point
        scatter(endPoints(iQuery,1),endPoints(iQuery,2),'filled','y');
        idxFoundPerQuery = vertcat(idx{iQuery}); % indices of endPoints
        dFoundPerQuery = vertcat(d{iQuery});
        % plot found points within radius of query
        scatter(endPoints(idxFoundPerQuery,1),endPoints(idxFoundPerQuery,2),'filled','r');
        % plot the distances calculated for all found points.
        arrayfun(@(i) text(endPoints(idxFoundPerQuery(i),1),endPoints(idxFoundPerQuery(i),2),num2str(dFoundPerQuery(i)),'color','y'),1:length(idxFoundPerQuery));
        saveas(gcf,[num2str(iQuery,'%03d') '.tif']);
        saveas(gcf,[num2str(iQuery,'%03d') '.fig']);
        
    end
    
end

% format edges: column 1 is indices of endpoint1 (Vertex1) and column 1 is
% the indice of endpoint2
E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:length(endPoints), 'UniformOutput', false);
E = vertcat(E{:});
d = vertcat(d{:});

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

%% Remove all overlap with previous pixels
% find path between edges
paths=arrayfun(@(i) gcaBresenham([endPoints(E(i,1),1) endPoints(E(i,1),2)], [endPoints(E(i,2),1) endPoints(E(i,2),2)]),...
    1:length(E(:,1)),'uniformoutput',0);

pathsidx = cellfun(@(x) sub2ind([ny,nx],x(2:end-1,2),x(2:end-1,1)),paths,'uniformoutput',0);
% check for overlap with candidateMask
candPix = find(candidateMask==1);
noOverlap = cellfun(@(x) isempty(intersect(x,candPix)),pathsidx); % have no interesction with the candidate mask
E = E(noOverlap,:);
d = d(noOverlap);

%% sanity plot 2 % plot all good possible connections with the distance labels
% (ie after self association,redundancy filter, and overlap filter)
sanityCheck2 = 0;
if sanityCheck2 == 1
    figure;
    imshow(labelMat>0,[]);
    hold on
    % should eventually color code by cost... but for now
    %cmap = jet(length(E));
    % easier to keep track of indexing if just use a for loop
    for iEdge = 1:length(E);
        % coords of first point of edge
        xCoord1 = endPoints(E(iEdge,1),1);
        yCoord1 = endPoints(E(iEdge,1),2);
        
        %    % coords of second point of edge
        xCoord2 = endPoints(E(iEdge,2),1);
        yCoord2 = endPoints(E(iEdge,2),2);
        %
        %
        %
        hold on
        %
        %   plot([xCoord1,xCoord2],[yCoord1,yCoord2],'color',cmap(iEdge,:));
        plot([xCoord1,xCoord2],[yCoord1,yCoord2],'color','y');
        %    % plot the distance at the query point
        %  %  text(xCoord1,yCoord1,num2str(d(iEdge),3),'color','y');
        %
    end
    %scatter(endPoints(E(:,1),1),endPoints(E(:,1),2),'filled','g'); % plot point 1 of edge
    %scatter(endPoints(E(:,2),1),endPoints(E(:,2),2),'filled','y'); % plot all point 2 of
    
end


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
    
    %% Filter for long distance links by geometry
    idxFilt = (d>ip.Results.MaxRadiusNoGeoTerm & (dotProd21' <=0 | dotProd12' <=0));
    ERemove = E(idxFilt,:);
    
    %% Plot the bad link sanity check if user desires
    if sanityCheck2 == 1
        hold on
        % easier to keep track of indexing if just use a for loop
        for iEdge = 1:length(ERemove);
            %    % coords of first point of edge
            xCoord1 = endPoints(ERemove(iEdge,1),1);
            yCoord1 = endPoints(ERemove(iEdge,1),2);
            %
            %    % coords of second point of edge
            xCoord2 = endPoints(ERemove(iEdge,2),1);
            yCoord2 = endPoints(ERemove(iEdge,2),2);
            %
            %
            %
            hold on
            %
            %   plot([xCoord1,xCoord2],[yCoord1,yCoord2],'color',cmap(iEdge,:));
            plot([xCoord1,xCoord2],[yCoord1,yCoord2],'color','r','Linewidth',2);
        end
        %plot the remove vectors in red
        idx21Bad = find(d>ip.Results.MaxRadiusNoGeoTerm & (dotProd21' <=0));
        idx12Bad = find(d>ip.Results.MaxRadiusNoGeoTerm & (dotProd12'<=0));
        if (~isempty(idx21Bad) || ~isempty(idx12Bad))
            quiver(endPoints(E(idx21Bad,2),1),endPoints(E(idx21Bad,2),2),...
                vect(E(idx21Bad,2),1), vect(E(idx21Bad,2),2),'color','g');
            
            quiver(endPoints(E(idx12Bad,1),1),endPoints(E(idx12Bad,1),2),...
                vect(E(idx12Bad,1),1), vect(E(idx12Bad,1),2),'color','m');
            
            
            %strDot12 = arrayfun(@(x) num2str(dotProd12,2),1:length(dotProd12),'uniformoutput',0);
            
            
            arrayfun(@(i) text(endPoints(E(idx12Bad(i),1),1),...
                endPoints(E(idx12Bad(i),1),2),num2str(dotProd12(idx12Bad(i)),2),...
                'FontSize',12,'Color','g'),1:length(idx12Bad));
            % plot the distance at the query point
            %  %  text(xCoord1,yCoord1,num2str(d(iEdge),3),'color','y');
            %arrayfun(@(i) text(endPoints(E(idx21Bad(i),1),1),...
            arrayfun(@(i) text(endPoints(E(idx21Bad(i),2),1),...
                endPoints(E(idx21Bad(i),2),2),num2str(dotProd21(idx21Bad(i)),2),...
                'FontSize',12,'Color','m'),1:length(idx21Bad));
            %
        end
    end
    %
    E(idxFilt,:) = [];
    d(idxFilt) = [];
    dotProd21(idxFilt) = [];
    dotProd12(idxFilt) = [];
end % isempty E
%%
if ~isempty(E) %check if there are reasonable edges
    % make sure to make d so that minimum d are favored
    d = max(d)-d; % max with =0
    d = d./max(d); % normalize - likely not completley correct need to work out..
    
    
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
    labelMatPostConnect = labelMat; % initiate a labelMatrix that connects the two pieces
    
    EPCandidateSortPostConnect = cellfun(@(x) x(:,1:2),EPCandidateSort,'uniformoutput',0); % take out the vectors. 
    % for all edges 
    for i =  1:length(E(:,1))
        % get the label of endpoint 1
        label1 = labelMatPostConnect(sub2ind([ny,nx],endPoints(E(i,1),2), endPoints(E(i,1),1)));
        % get the label of endpoint 2 in linkage
        label2 = labelMatPostConnect(sub2ind([ny,nx],endPoints(E(i,2),2),endPoints(E(i,2),1)));
        test = [label1 ;label2];
        
        % keep the maxvalue label (arbitrary assignment)
        labelKeep  = max(test);
        labelSwitch = min(test);
        % change the label of the portion of the connected piece we chose
        % to switch
        labelMatPostConnect(labelMatPostConnect==labelSwitch) = labelKeep; % 
        
        % add the same label to the path between the two pieces
        labelMatPostConnect(sub2ind([ny,nx],paths{i}(:,2),paths{i}(:,1)))= labelKeep;
        
        % make new mask of the new connected filo and get endpoints
        testMask = double(labelMatPostConnect==labelKeep);
        sumKernel = [1 1 1];
        % find and record the new EPs for the CC
        newEPs = double((testMask.* (conv2(sumKernel, sumKernel', padarrayXT(testMask, [1 1]), 'valid')-1))==1);
        [ye,xe] = find(newEPs~=0);
        %     coords(:,1) = xe;
        %     coords(:,2) = ye;
        EPCandidateSortPostConnect{labelKeep} = [xe ye];
        % keep in same format therefore simply remove all the endpoint
        % information corresponding to the label you eradicated. 
        EPCandidateSortPostConnect{labelSwitch} = [NaN NaN ;NaN NaN]; % take out those endpoints that were previously considered
        
    end
    status = 1; % there were links
else
    candidateMaskNew = candidateMask;
    labelMatPostConnect = labelMat;
    status = 0; % no links
    
    
    
end

