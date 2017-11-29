function [outputMasks, filoInfo,status,pixIdxCandsUnMatched,candFiloEPsUnMatched,TSFigs] = gcaConnectFiloBranch(inputPoints, candFiloEPs, pixIdxCands, labelMatSeedFilo,filoInfo, maxRes,maxTh,img,normalsC,smoothedEdgeC,varargin)
% gcaConnectFiloBranch:

%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUIRED
%
%  %% SEED INFORMATION:
%     (ie High confidence regions to which the candidates will be attached)  %%
%
%  inputPoints: (REQUIRED) Rx2 double array
%     documenting the seed locations
%     where R is the number of points
%     and 2 is the vector of x,y coordinates respectively
%
%  filoInfo: (REQUIRED) 1xC structure array
%     with many flexible fields describing each individual ridge structure
%     of the filopodia/branch network, including each
%     filopodia's correponding coordinates.
%     C is the number of filopodia ridge structures currently documented,
%     hence each c value gives information as to the filo objects label
%     number used in labelMatSeedFilo
%     This structure will be appended in this function as the
%     filopodia/branch network is reconstructed.
%     (Output of gcaRecordFilopodiaSeedInformation.m, and
%      gcaConnectFiloBranch.m (ie current function))
%
%  labelMatSeedFilo: (REQUIRED) RxC double array
%     where R is the size y  of the original image and C
%     where each filopodia is labeled- note this format does not allow for crossings
%     by definition- as each filopodia coordinate can only at  max have one label.
%     It is created by the filoInfo: The veil/stem boundary always
%     has a label of 1, while each filopodia is labeled according to the
%     the number in filoInfo structure.
%
%  %% CANDIDATE FILOPODIA INFORMATION: %%
%
%  candFiloEPs: (REQUIRED) 1xC cell of 2x4 double arrays
%     where C is the number of possible candidate filopdia (ridges) to be attached.
%     and each 2x4 double array holds the
%     [coordX, coordY, vectX, vectY] where coordX/coordY denote the endpoint
%     of the candidate filopodia and vectX,vectY specificies the
%     direction toward the endpoint to be used for the cost
%
% pixIdxCands: (REQUIRED)
%      the corresponding pixIndices of the candidates
%
% maxRes : (maybe able to remove - no longer collecting the response for fitting.
%
% maxTh : (REQUIRED)
%
%% This needs to be replaced by a cellArray for labels 20150604
% labelCandidates: rxc double
%%
%
% Output: outputMasks: structure containing fields
%                      .finalReconstruct (binary mask of the final
%                      reconstruction)
%                      .links (binary mask of the links added)
%                      .candFiloAdded.Body = mask with candidate filo added
%                      to body
%                      .candFiloAdded.Branch = mask with candidate filo
%                      added to another filo
% documentInfo: flag to either collect the new information in a filoInfo
% datastruct or not 1 = collect/ 0 = no collect (useful in the case of
% simple recostructions for internal filo
%         filoInfo: the corresponding filoInfo data structure with the new filopodia information for
%                   each candidate filopodia attached added
%% CHECK INPUT
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('inputPoints');
ip.addRequired('candFiloEPs');
ip.addRequired('pixIdxCands');
ip.addRequired('labelMatSeedFilo');
ip.addRequired('filoInfo');
ip.addRequired('maxRes');
ip.addRequired('maxTh');
ip.addRequired('img');
ip.addRequired('normalsC');
ip.addRequired('smoothedEdgeC');

ip.addParameter('maxRadiusConnectFiloBranch',5);
ip.addParameter('geoThreshFiloBranch',0.5);
ip.addParameter('TSOverlays',false);
ip.addParameter('cMapLength',4096);
ip.addParameter('cMapType','parula');
ip.parse(inputPoints,candFiloEPs,pixIdxCands,labelMatSeedFilo,filoInfo,maxRes,maxTh,img,normalsC,smoothedEdgeC,varargin{:});
p = ip.Results;

%% INITIATE
dims = size(labelMatSeedFilo);
%dims = size(cellBoundary);
out = zeros(dims);
queryPoints = vertcat(candFiloEPs{:});
nq = size(queryPoints,1);
countFig = 1;

candFiloEPsUnMatched = candFiloEPs;
vectCand = cellfun(@(x) x(:,3:4),candFiloEPs,'uniformoutput',0);
% take out the unit vectors
candFiloEPs = cellfun(@(x) x(:,1:2),candFiloEPs,'uniformoutput',0);
pixIdxCandsUnMatched = pixIdxCands;
TSFigs = [];
%% WORKING : To Remove
% 20150604 for now just make a workshift labelCandidates
% note don't want to keep this because we know that the labels need to
% be assigned to more than one pixel which will not be the case here.
labelCandidates = zeros(dims);
for iCand = 1:numel(pixIdxCands)
    labelCandidates(pixIdxCands{iCand})= iCand;
end
%% START
% note the output of KDTreeBall is of the form of a cell array n = the
% number of query points, first parameter is the index of the input points
% around that query within the given search radius. so
% idxInput{iQuery} = [idxInput1, idxInput2, idxIndput3...]);
%candFiloEPs = candFiloEPs(:,1:2);

% Here the InputPoints are all points along the seed.
% We find all input points (seedPoints) within the given maxRadius input by
% the user
% We loop here as each candFilo has more than one end point

for iCanFilo = 1:numel(candFiloEPs)
    [idxInput, dist] = KDTreeBallQuery(inputPoints, candFiloEPs{iCanFilo}, ip.Results.maxRadiusConnectFiloBranch);
    node{iCanFilo} = idxInput; % combine all possible attachments for the candidateFilo
    distNode{iCanFilo}=dist;
end
testMatch =  cellfun(@(x) ~isempty(vertcat(x{:})),node); % if all are empty no matches within Dist

if sum(testMatch) ~= 0 ;
    
    %% Assign Edges: Each candidate filo has N competing attachment sites
    
    % Each candidate is labeled as a node: paths are searched on either end
    % but only a single path is retained
    
    % generate edge map
    % the query point has a possible path or "edge" with each of the input
    % points within the search radius.
    E = arrayfun(@(i) [repmat(i, [length(vertcat(node{i}{:})) 1]) vertcat(node{i}{:})], 1:numel(node), 'UniformOutput', false);
    E = vertcat(E{:});
    % need to label which end of the filo the coordinate corresponds
    idxEnd = arrayfun(@(i) [ repmat(1,numel(node{i}{1}),1) ; repmat(2,numel(node{i}{2}), 1)],1:numel(candFiloEPs),'Uniformoutput',false);
    idxEnd = vertcat(idxEnd{:});
    %
    E = [E idxEnd];
    
    %Collect all distance values
    D = arrayfun(@(i) vertcat(distNode{i}{:}),1:numel(distNode),'uniformoutput',0);
    D = vertcat(D{:});
    
    % calculate the orientation of each candidate filo
    seedPtsx = inputPoints(E(:,2),1);
    seedPtsy = inputPoints(E(:,2),2);
    labelsSeed = labelMatSeedFilo(sub2ind(dims,seedPtsy,seedPtsx));
    
    seedPtsx(labelsSeed==0) =[];
    seedPtsy(labelsSeed==0) = [];
    
    E(labelsSeed==0,:) = [];
    D(labelsSeed==0,:) = [];
    
    labelsSeed(labelsSeed==0)=[];% get rid of inconsitencies between the two inputs
    
    %% get the normals of the cellEdge.
    
    %%%HERE IS QUITE A LARGE   PROBLEM WE HAVE: don't want to have the same
    %%%orienation criteria for ALL attachments endOn is different than body
    % attatchments.
    dotCandAndSeed = zeros(length(labelsSeed),1);
    for iPath = 1 :length(labelsSeed) % over all potential connection paths
        maskDisk = zeros(dims);
        maskDisk(sub2ind(dims,seedPtsy(iPath),seedPtsx(iPath)))= 1;
        maskDisk = imdilate(maskDisk,strel('disk',4));
        maskFilo = zeros(dims);
        maskFilo(labelMatSeedFilo==labelsSeed(iPath)) =1;
        maskTest = maskFilo.*maskDisk;
        [y,x] = ind2sub(dims,find(maskTest==1));
        
        vectInputC = [(x(1)-x(end)) , (y(1)-y(end))];
        vectInput(iPath,:) = vectInputC;
        dInput = sqrt((x(1)-x(end))^2 + (y(1)-y(end))^2);
        vectorCand = vectCand{E(iPath,1)}(E(iPath,3),:);
        dotCandAndSeed(iPath) = abs(dot(vectInputC,vectorCand)./dInput);
        %                  dotCandAndSeed(iPath) = abs(dot(vectInput,vectorCand)./dInput./dCand{E(iPath,1)});
    end
    
    if ip.Results.TSOverlays == true;
        TSFigs(countFig).h= setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Seed_and_Candidates';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        
        imshow(-img,[]);
        hold on
        spy(labelMatSeedFilo,'b');
        spy(labelCandidates,'m');
        text(5,10,'Seed','FontSize',10,'color','b');
        text(5,25,'Candidate Ridges','FontSize',10,'color','m');
        countFig = countFig+1;
        
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Seed_and_Candidates_With_Vectors';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        imshow(-img,[]);
        hold on
        spy(labelMatSeedFilo,'b',5);
        % scatter(seedPtsx(:),seedPtsy(:),5,'y','filled');
        spy(labelCandidates,'r',5);
        xyall = vertcat(candFiloEPs{:});
        vectall = vertcat(vectCand{:});
        quiver(xyall(:,1),xyall(:,2),vectall(:,1),vectall(:,2),0.2,'r');
        scatter(seedPtsx(:),seedPtsy(:),5,'y','filled');
        % old vectors
        quiver(seedPtsx(:),seedPtsy(:),vectInput(:,1),vectInput(:,2),0.2,'y');
        
        countFig = countFig  +1;
    end
    
    % convert indexing to node labels to input into graph matching
    % For each Filo have N Total possible edge matches (K from side 1 and L
    % from side 2) Only possible attachment sites with the search radius are
    % considered. In the end only the maximum cost link will remain (cost of linking will
    % be defined below)
    [candFiloNodes,~,nodeLabels] = unique(E(:,1),'stable'); % reason note: some of the filo will not be candidates as their endpoints are not within the given radius
    NNodeQuery = length(candFiloNodes);
    [inputLinks,~,nodeLabelsInput] = unique(E(:,2),'stable'); % reason note: just in case two filo are competing over the same seed point
    nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
    EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
    numberNodes = length(inputLinks) + length(candFiloNodes);
    
    % DEFINE COST FUNCTION
    % want a link that 1) have a similar orientation to the filo,
    % 2) have a high amount of image intensity, and 3) is a relatively short
    % distance (though this maybe should be slightly less emphasized)
    
    for iPath = 1:size(E,1)
        xCoordQuery = candFiloEPs{E(iPath,1)}(E(iPath,3),1);
        yCoordQuery = candFiloEPs{E(iPath,1)}(E(iPath,3),2) ;
        xCoordInput = inputPoints(E(iPath,2),1);
        yCoordInput = inputPoints(E(iPath,2),2);
        iSeg{iPath} = bresenham([xCoordQuery yCoordQuery], [xCoordInput yCoordInput]);
        int{iPath}  = img(sub2ind(dims, iSeg{iPath}(:,2), iSeg{iPath}(:,1)));
        vectTest = [(xCoordInput-xCoordQuery), (yCoordInput-yCoordQuery)];
        % value of 1 : two segments well aligned (angle) between them = 0
        % as we want alignment of added attatchments
        dTest = D(iPath,1);
        vectorCand = vectCand{E(iPath,1)}(E(iPath,3),:);
        dotProd(iPath,1) = dot(vectorCand,vectTest)/dTest;
    end
    %% COST Definition
    % Distance Term: make so that a lower distance = positive weight
    D = max(D)-D;
    % normalize
    D= D./max(D);
    
    % Intensity Term
    int = cellfun(@(x) mean(x), int); % higher this intensity the more likely to like
    maxIntPath = max(int);
    minIntPath = min(int);
    normInt = (int - minIntPath)./(maxIntPath-minIntPath);
    normInt = normInt';
    %% Final Cost
    % Orientation Term:  above 1 favored 0 unfavored
    costTotal = (0.5*D + normInt +dotProd+dotCandAndSeed);
    %% Perform the sanity checks for the cost terms
    
    if ip.Results.TSOverlays
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Potential_Paths_with_Cost_D';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        
        [idxCMapDist] = gcaPlotLinksByCost(img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,D,ip.Results.cMapLength, ...
            'cMapType',ip.Results.cMapType);
        text(10,5,'Score By Distance');
        idxCMapCell{1} = idxCMapDist;
        figNames{1} = 'Potential_Paths_with_Cost_D';
        countFig = countFig+1;
        
        %%   if ip.Results.TSOverlays == true
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Potential_Paths_with_Cost_Int';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        
        [idxCMapInt] = gcaPlotLinksByCost(img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,normInt,ip.Results.cMapLength, ...
            'cMapType',ip.Results.cMapType);
        idxCMapCell{2} = idxCMapInt;
        figNames{2} = 'Potential_Paths_with_Cost_Int';
        
        text(10,5,'Score By Intensity Mean');
        
        countFig = countFig+1;
        %%
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Potential_Paths_with_Cost_CandLinkerGeo';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        
        [idxCMapCandLink] = gcaPlotLinksByCost(img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,dotProd,ip.Results.cMapLength, ...
            'cMapType',ip.Results.cMapType);
        text(10,5,'Score By Candidate Linker Geometry')
        idxCMapCell{3} = idxCMapCandLink;
        figNames{3} = 'Potential_Paths_with_Cost_CandLinkerGeo';
        
        countFig = countFig+1;
        %% Cand And Seed Value will be from 0 to 1
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Potential_Paths_with_Cost_CandSeedGeo';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        
        [idxCMapCandSeed] = gcaPlotLinksByCost(img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,dotCandAndSeed,ip.Results.cMapLength, ...
            'cMapType',ip.Results.cMapType);
        text(10,5,'Score By Candidate Seed Geometry');
        idxCMapCell{4} = idxCMapCandSeed;
        figNames{4} = 'Potential_Paths_with_Cost_CandSeedGeo';
        
        countFig = countFig+1;
        
        %% Total
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Potential_Paths_with_Cost';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch';
        
        [idxCMap] = gcaPlotLinksByCost(img,labelCandidates,labelMatSeedFilo,candFiloEPs,seedPtsx,seedPtsy,iSeg,costTotal,ip.Results.cMapLength, ...
            'cMapType',ip.Results.cMapType);
        text(10,5,'Final Score');
        idxCMapCell{5} = idxCMap;
        figNames{5} = 'Potential_Paths_with_Cost';
        
        
        countFig = countFig+1;
    end % if TSOverlays
    
    %% Filter links by geometery of candidate and linker
    E = [E costTotal D normInt dotProd dotCandAndSeed];
    
    E = E(dotProd>ip.Results.geoThreshFiloBranch,:);
    EFinal = EFinal(dotProd>ip.Results.geoThreshFiloBranch,:);
    costTotal = costTotal(dotProd>ip.Results.geoThreshFiloBranch);
    if ip.Results.TSOverlays
        idxCMap =  idxCMap(dotProd>ip.Results.geoThreshFiloBranch,:);
        idxCMapCell = cellfun(@(x) x(dotProd>ip.Results.geoThreshFiloBranch,:),idxCMapCell,'uniformoutput',0);
        
    end
    iSeg= iSeg(:,dotProd>ip.Results.geoThreshFiloBranch);
    
    %% Always take out links that cross the veil stem
    
    iSegLinIdx = cellfun(@(x) sub2ind(dims,x(:,2),x(:,1)),iSeg,'uniformoutput',0);
    
    veilMask = zeros(dims);
    linIdx = sub2ind(dims,inputPoints(:,2),inputPoints(:,1));
    % make sure not NaN
    linIdx = linIdx(~isnan(linIdx));
    veilMask(linIdx) = 1;
    veilMask = imfill(veilMask,'holes');
    veilMask(linIdx) = 0;
    linIdxVeilStem= find(veilMask);
    noOverlapLinks1 = cellfun(@(x) isempty(intersect(x,linIdxVeilStem)),iSegLinIdx);
    
    E = E(noOverlapLinks1,:);
    EFinal = EFinal(noOverlapLinks1,:);
    costTotal = costTotal(noOverlapLinks1);
    if ip.Results.TSOverlays == 1
        idxCMap =  idxCMap(noOverlapLinks1,:);
        idxCMapCell = cellfun(@(x) x(noOverlapLinks1,:),idxCMapCell,'uniformoutput',0);
    end
    iSeg= iSeg(noOverlapLinks1);
    
    %% OPTIONAL : Do NOT allow linkers to cross other candidate filopodia (can only introduce a cross by crossing a seed)
    noOverlap = 1;
    
    if noOverlap == 1
        
        iSegLinIdx = cellfun(@(x) sub2ind(dims,x(:,2),x(:,1)),iSeg,'uniformoutput',0);
        
        candPix = vertcat(pixIdxCands{:});
        % Fill along the cadidate pixels to make sure overlapping
        % links are removed
        diagFill = zeros(dims);
        diagFill(candPix) = 1;
        diagFill = bwmorph(diagFill,'diag');
        candPix = find(diagFill);
        % Remove endpoints (these are connected to the linker in all cases)
        cEPsAll = vertcat(candFiloEPs{:});
        [linIdxCandEPs]  = sub2ind(dims,cEPsAll(:,2),cEPsAll(:,1));
        candPix   = setdiff(candPix,linIdxCandEPs);
        
        noOverlapLinks = cellfun(@(x) isempty(intersect(x,candPix)),iSegLinIdx); % have no interesction with the candidate mask
        E = E(noOverlapLinks,:);
        EFinal = EFinal(noOverlapLinks,:);
        costTotal = costTotal(noOverlapLinks);
        if ip.Results.TSOverlays == 1
            idxCMap =  idxCMap(noOverlapLinks,:);
            idxCMapCell = cellfun(@(x) x(noOverlapLinks,:),idxCMapCell,'uniformoutput',0);
        end
        iSeg= iSeg(noOverlapLinks);
        
    end % if no overlap
    %%
    [candFiloNodes,~,nodeLabels] = unique(EFinal(:,1),'stable'); % reason note: some of the filo will not be candidates as their endpoints are not within the given radius
    NNodeQuery = length(candFiloNodes);
    [inputLinks,~,nodeLabelsInput] = unique(EFinal(:,2),'stable'); % reason note: just in case two filo are competing over the same seed point
    nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
    EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
    numberNodes = length(inputLinks) + length(candFiloNodes);
    
    %% Perform Weighted Graph Matching if Edges Exist
    
    if ~isempty(EFinal)
        %% TS Overlays After Geometry Thresholds : Histograms of Cost Components
        if ip.Results.TSOverlays
            TSFigs(countFig).h = setAxis('on');
            TSFigs(countFig).name = 'Histograms_of_Cost_Parameters_AfterGeoThresh';
            TSFigs(countFig).group = 'Reconstruct_FiloBranch';
            nPaths = size(E,1);
            for i = 1:5
                [n{i},center{i}] = hist(E(:,3+i),50);
                minVal(i) = min(E(:,3+i));
                maxVal(i) = max(E(:,3+i));
            end
            
            allN = horzcat(n{:});
            maxYVal = max(allN./nPaths);
            %maxYVal = allN(:);
            
            subplot(5,1,1);
            
            bar(center{1},n{1}/nPaths);
            xlabel('Cost Total');
            axis([ip.Results.geoThreshFiloBranch,3.5,0,maxYVal]);
            title(['Min ' num2str(minVal(1),3) ' Max ' num2str(maxVal(1),3)]);
            
            subplot(5,1,2);
            
            bar(center{2},n{2}/nPaths);
            xlabel('Distance');
            axis([0,1,0,maxYVal]);
            title(['Min ' num2str(minVal(2),3) ' Max ' num2str(maxVal(2),3)]);
            
            subplot(5,1,3);
            
            bar(center{3},n{3}/nPaths);
            xlabel('Mean Intensity');
            axis([0,1,0,maxYVal]);
            title(['Min ' num2str(minVal(3),3) ' Max ' num2str(maxVal(3),3)]);
            
            subplot(5,1,4);
            
            bar(center{4},n{4}/nPaths);
            xlabel('Geometry With Linker');
            axis([ip.Results.geoThreshFiloBranch,1,0,maxYVal]);
            title(['Min ' num2str(minVal(4),3) ' Max ' num2str(maxVal(4),3)]);
            
            subplot(5,1,5);
            
            bar(center{5},n{5}/nPaths);
            xlabel('Geometry Candidate and Seed');
            axis([0,1,0,maxYVal]);
            title(['Min ' num2str(minVal(5),3) ' Max ' num2str(maxVal(5),3)]);
            countFig = countFig+1;
        end
        %% TS Overlays After Geometry Thresholds : Color Code By Cost
        
        if ip.Results.TSOverlays
            
            for iCost = 1:5
                
                TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
                TSFigs(countFig).name = [figNames{iCost} '_filterGeo'];
                TSFigs(countFig).group = 'Reconstruct_FiloBranch';
                
                imshow(-img,[]);
                
                hold on
                text(5,10,['Geometry Threshold' num2str(ip.Results.geoThreshFiloBranch)],'FontSize',10,'Color','k');
                
                spy(labelMatSeedFilo>0,'k');
                spy(labelCandidates>0,'k',5);
                scatter(seedPtsx(:),seedPtsy(:),'k','filled');
                allCandEPs = vertcat(candFiloEPs{:});
                scatter(allCandEPs(:,1),allCandEPs(:,2),'k','filled');
                cMapType = str2func(ip.Results.cMapType);
                %
                cMapLength=ip.Results.cMapLength; cMap= cMapType(cMapLength);
                
                for k = 1:length(cMap);
                    if sum(idxCMapCell{iCost}(:,1)==k)~=0
                        toPlot = iSeg(idxCMapCell{iCost}(:,1) == k);
                        
                        cellfun(@(x) plot([x(1,1),x(end,1)],[x(1,2),x(end,2)],'color',cMap(k,:)),toPlot);
                        clear toPlot
                    else
                    end
                    % make colormap of costs.
                    
                    % show each segment in iSeg cell plotted by the costTotal colo
                    %
                end
                countFig = countFig +1;
            end % for iCost
        end % ip.Results.TSOverlays
        %% Matching
        M = maxWeightedMatching(numberNodes, EFinal, costTotal);
        % check for double labels
        % convertBack
        % E = [candFiloNodes(nodeLabels) inputLinks(nodeLabelsInput)]; % convert back to original indices of input and query points
        E = E(M,:);% get those edges that matched (from original indexing)
    end   % isempty
    
    %% Update the unmatched candidate list and prepare output
    
    if ~isempty(E) % might be empty now if all were repeats
        pixIdxCandsUnMatched(E(:,1)) = []; % take out the matched candidates
        candFiloEPsUnMatched(E(:,1)) =[];
        % add linear segments corresponding to linked endpoints
        goodConnect = iSeg(M);
        % convert to logical indexing
        pixGoodConnect = cellfun(@(i) sub2ind(dims,i(:,2),i(:,1)), goodConnect,'uniformoutput',0);
        out(vertcat(pixGoodConnect{:}))= 1;
    end
    %% Update the final results
    links = out;
    outputMasks.links = links;
    allInputMask= (labelCandidates>0|labelMatSeedFilo>0);
    out = double(out|allInputMask);
    outputMasks.finalReconstruct = out;
    
    %% start documenting data
    labelInputCon = zeros(length(E(:,1)),1);
    % get the filoIdx of those filo to which attachments have been made
    % NOTE : remember that if E is empty length will
    % be zero and by small favors this doesn't error in matlab when one has
    % 1:0 (though this is likely not good practice) it will just skip as we
    % would like it to 
    for iMatch = 1:length(E(:,1))
        labelInputCon(iMatch) = labelMatSeedFilo(sub2ind(dims,inputPoints(E(iMatch,2),2),inputPoints(E(iMatch,2),1)));
       
        labelCandCon(iMatch) = E(iMatch,1);
        
        if E(iMatch,3) == 1;
            EPsCand(iMatch) = 2 ;
        else
            EPsCand(iMatch) = 1;
        end
    end % iMatch
    
    % veil/stem is labeled with 1 so
    % indexing is off by 1
    labelInputCon = labelInputCon -1;% change back to 1
    
    %% Document Veil Stem Attachments
    idxBodyAttach = find(labelInputCon==0);
    
    if ~isempty(idxBodyAttach)
        num = length(idxBodyAttach);
        % record the new attachments to the body
        for iFilo = 1:num
            % quickfix for bug
            if labelCandCon(idxBodyAttach(iFilo))~=0 % proceed otherwise don't count
                EPNum = EPsCand(idxBodyAttach(iFilo));
                xEP = candFiloEPs{E(idxBodyAttach(iFilo),1)}(EPNum,1);
                yEP = candFiloEPs{E(idxBodyAttach(iFilo),1)}(EPNum,2);
                
                testMask = zeros(dims);
                %  testMask(labelMatSeedFilo==1|labelMatSeedFilo==2) = 1; % get neuriteBodyAll
                % put together a mask of
                % The seed
                testMask(labelMatSeedFilo==1)=1;
                
                testMask(vertcat(pixIdxCands{labelCandCon(idxBodyAttach(iFilo))})) = 1;
                
                % The Connection
                testMask(pixGoodConnect{idxBodyAttach(iFilo)}) = 1;
                
                testMask = bwmorph(testMask,'thin','inf'); %
                testMask = logical(testMask);
                transform = bwdistgeodesic(testMask,xEP,yEP);
                % now need to get pixels
                pixIdxBack = nan(50,1); % overinitialize to make happy
                
                iPix = 1;
                while length(find(transform==iPix)) == 1
                    pixIdxBack(iPix) = find(transform==iPix); %
                    iPix = iPix +1;
                end
                
                pixIdxBack = pixIdxBack(~isnan(pixIdxBack));
                if ~isempty(pixIdxBack)
                    
                    % get coords of the these pixels
                    [yBack,xBack]= ind2sub(size(out),pixIdxBack);
                    
                    verticesEP = [yEP xEP];
                    edgePathCoord{1} = [yBack, xBack];
                    
                    x = walkFiloForAndBack([], verticesEP,edgePathCoord,maxTh,maxRes,img,0,10);
                    x.type = 0; % label a primary filo
                    x.cross = 0;
                    x.groupCount = max(vertcat(filoInfo(:).groupCount)) +1; % start a group
                    %x.bodyType = iBody;
                    if isfield(filoInfo,'conIdx')
                        x.conIdx = []; % just to make sure fields are consistent
                        x.conXYCoords= [];
                    end
                    
                    if ~isempty(smoothedEdgeC)
                        
                        %%%% CALCULATE THE LOCAL ORIENTATION %%%%
                        baseFilo = [xBack(end),yBack(end)];
                        % region identification)
                        [idx, dist] = KDTreeBallQuery(smoothedEdgeC,baseFilo,3);
                        idx = idx{:};
                        
                        % Note can either get orientation of filo from the maxTh output of the steerable filter or from just calculating
                        % a small local vector. we will see which one is cleaner ... so far I tend
                        % to favor the small vector...
                        
                        % would potentially add a normalsC rotated.
                        %% TESTING 20150616
                        %20150616 POTENTIALLY CHANGE? Note maybe want to
                        %the rotation on the average vector?
                        avgNormLocal = mean(normalsC(idx,:),1);% might want to change to a majority?
                        pathCoords = [yBack xBack];
                        
                        %%
                        sanityCheck =0;
                        if sanityCheck == 1
                            
                            imshow(-img,[]);
                            hold on
                            
                            roiYX = bwboundaries(labelMatSeedFilo==1);
                            cellfun(@(x) plot(x(:,2),x(:,1),'b'),roiYX);
                            scatter(smoothedEdgeC(idx,1),smoothedEdgeC(idx,2),10,'b','filled');
                            quiver(smoothedEdgeC(idx(3),1),smoothedEdgeC(idx(3),2), avgNormLocal(1),avgNormLocal(2),10,'filled','color','c','Linewidth',2);
                            quiver(smoothedEdgeC(idx,1),smoothedEdgeC(idx,2),normalsC(idx,1),normalsC(idx,2),'color','g');
                        end
                        %%
                        % test length of pathCoords
                        pixFilo = size(pathCoords,1);
                        if pixFilo <= 4
                            back = pixFilo-1;
                        else
                            back = 4;
                        end
                        
                        if back ~=0
                            localVectFilo = [pathCoords(end-back,2)-pathCoords(end-1,2),pathCoords(end-back,1)-pathCoords(end-1,1)];
                            vectLength = sqrt((pathCoords(end-back,2)-pathCoords(end-1,2)) ^2 + (pathCoords(end-back,1) - pathCoords(end-1,1))^2);
                            normLength = sqrt(avgNormLocal(1)^2 + avgNormLocal(2)^2);
                            cosAngle = dot(avgNormLocal(1:2),localVectFilo)/vectLength/normLength;
                            angleToBody = acosd(cosAngle);
                            %angleToBody = 180- angle -90;
                        else
                            angleToBody = NaN;
                            localVectFilo = NaN;
                        end
                        x.orientation = angleToBody;
                        x.localVectAttach = avgNormLocal; % for now just save the normal vector
                        x.localVectFilo= localVectFilo;
                    else
                        x.orientation = [];
                        x.localVectAttach = [];
                        x.localVectFilo = [];
                    end
                    
                    idxBodyAttachFiloInfo = numel(filoInfo) +1;
                    
                    %Fix fieldnames to match (might not be a problem in later versions
                    %if initialize data structure more appropriately: we'll see if it's
                    %worth fixing
                    fieldsx = fieldnames(x);
                    fieldsFiloInfo = fieldnames(filoInfo);
                    fieldsAddtoX = setdiff(fieldsFiloInfo,fieldsx);
                    
                    for i = 1:length(fieldsAddtoX)
                        x.(char(fieldsAddtoX(i))) = NaN;
                    end
                    
                    x = orderfields(x); % put in alphabetical order
                    filoInfo = orderfields(filoInfo);
                    
                    filoInfo(idxBodyAttachFiloInfo) = x;
                    
                    candFiloAdded.Body{iFilo} = pixIdxBack;
                    clear x edgePathCoord xBack yBack
                    
                end % don't record if isempty (this sometimes happens if there are weird cycles in the response
            end % if labelCandCon
        end % iFilo
    end % ~isempty
    %% Document Filopodia Attachments
    idxFiloAttach = find(labelInputCon >0); % idx of those matches that are attached to a seed filo
    if ~isempty(idxFiloAttach)
        for i = 1:length(idxFiloAttach)
            test{i}   =  filoInfo(labelInputCon(idxFiloAttach(i))).Ext_pixIndicesBack(end);
        end
        pixEnd = vertcat(test{:});
        [y,x] = ind2sub(dims,pixEnd);
        coordsEnd = [x,y];
        % coordsEnd = vertcat(filoInfo(labelInputCon(idxFiloAttach)).endpointCoord);
        %  coordsEnd = [coordsEnd(:,2) coordsEnd(:,1)];
        %
        attachSite = inputPoints(E(idxFiloAttach,2),:);
        for i = 1:length(coordsEnd(:,1))
            idxEndAttach(i) = isequal(coordsEnd(i,:),attachSite(i,:));
        end
        
        %%%% REATTACH AND DOCUMENT END-ON ATTACHMENTS %%%%
        idxEndOnAttach = idxFiloAttach(idxEndAttach);
        
        attachSitesEndOn = attachSite(idxEndAttach,:);
        labelsEndon = labelMatSeedFilo(sub2ind(dims,attachSitesEndOn(:,2),attachSitesEndOn(:,1)));
        labelsEndon = labelsEndon-1;
        subtractEndOn = 0;
        
        for iEndon = 1:length(attachSitesEndOn(:,1))
            
            EPNum = EPsCand(idxEndOnAttach(iEndon));
            xEP = candFiloEPs{E(idxEndOnAttach(iEndon),1)}(EPNum,1);
            yEP = candFiloEPs{E(idxEndOnAttach(iEndon),1)}(EPNum,2);
            % first test if you are getting crap
            testMask = zeros(dims);
            if labelCandCon(idxEndOnAttach(iEndon)) ~=0
                testMask(vertcat(pixIdxCands{labelCandCon(idxEndOnAttach(iEndon))})) = 1;
                
                testMask(pixGoodConnect{idxEndOnAttach(iEndon)}) = 1; % get the links
                testMask(filoInfo(labelsEndon(iEndon)).Ext_pixIndicesBack) =1;
                
                nn = padarrayXT(double(testMask~=0), [1 1]);
                sumKernel = [1 1 1];
                nn = conv2(sumKernel, sumKernel', nn, 'valid');
                nn1 = (nn-1) .* (testMask~=0);
                junctMask = nn1>2;
                if sum(junctMask(:)) == 0;
                    
                    testMask(labelMatSeedFilo==1) = 1;
                    testMask = logical(testMask);
                    testMask = bwmorph(testMask,'thin','inf');
                    % test to make sure it doesn't still introduce a
                    % connectivity problem
                    transform = bwdistgeodesic(testMask,xEP,yEP);
                    % now need to get pixels
                    pixIdxBack = nan(50,1); % overinitialize to make happy
                    
                    iPix = 1;
                    while length(find(transform==iPix)) == 1
                        pixIdxBack(iPix) = find(transform==iPix);
                        iPix = iPix +1;
                    end
                    pixIdxBack = pixIdxBack(~isnan(pixIdxBack));
                    % get coords of the these pixels
                    [yBack,xBack]= ind2sub(size(out),pixIdxBack);
                    
                    verticesEP = [yEP xEP];
                    edgePathCoord{1} = [yBack, xBack];
                    x = walkFiloForAndBack([], verticesEP,edgePathCoord,maxTh,maxRes,img,0,10);
                    x.type = filoInfo(labelsEndon(iEndon)).type; % keep it whatever type it was previously
                    x.cross = 1; % mark it as a likely crossover.
                    
                    % keep any old field the same (ie internal filo, conIdx etc)
                    fieldsx = fieldnames(x);
                    fieldsFiloInfo = fieldnames(filoInfo);
                    fieldsAddtoX = setdiff(fieldsFiloInfo,fieldsx);
                    for iField = 1:length(fieldsAddtoX)
                        x.(char(fieldsAddtoX(iField))) = filoInfo(labelsEndon(iEndon)).(char(fieldsAddtoX(iField)));
                    end
                    
                    x = orderfields(x); % put in alphabetical order
                    filoInfo = orderfields(filoInfo);
                    
                    filoInfo(labelsEndon(iEndon)) = x; % rewritethe data
                    pixAdd =  pixIdxCands{labelCandCon((idxEndOnAttach(iEndon)))};
                    
                    linkAdd = pixGoodConnect{idxEndOnAttach(iEndon)};
                    candFiloAdded.EndOn{iEndon} = [pixAdd;linkAdd];
                    clear x edgePathCoord xBack yBack
                else % if % likely garbage don't include
                    % remove from links
                    
                    pixIndicesToRemove = pixGoodConnect{idxEndOnAttach(iEndon)};
                    pixGoodConnectTest = pixGoodConnect;
                    pixGoodConnectTest(idxEndOnAttach(iEndon)) = [];
                    
                    % check for attachment
                    savePix = cellfun(@(x) intersect(x,pixIndicesToRemove),pixGoodConnectTest,'uniformoutput',0);
                    
                    links(pixGoodConnect{idxEndOnAttach(iEndon)}) =0;
                    
                    if ~isempty(savePix)
                        % find those pixel that are in more than one link and save.
                        savePix = vertcat(savePix{:});
                        links(savePix) = 1;
                        %                         display('saving pixels')
                    end
                    
                    subtractEndOn = subtractEndOn+1;
                    % updata output masks
                    outputMasks.links = links;
                end
            else %
            end % if labelCon (label Problem fix)
        end % for iEndon
        
        %%%% ADD FILOPODIA BRANCHES TO THE STRUCTURE %%%%
        idxFiloAttach = idxFiloAttach(idxEndAttach==0);
        
        % record the reattachment on these filo
        for iFilo = 1:length(idxFiloAttach)
            idxSeedFilo = labelInputCon(idxFiloAttach(iFilo));
            if filoInfo(idxSeedFilo).type == 0 % if single make a primary
                % change the type label to mark as a primary branch stem
                filoInfo(idxSeedFilo).type  = 1; % mark it as the main part of a branch
            end
            
            idxBranch = numel(filoInfo) +1;
            % mark the branches connectivity information
            % matrix of 1 idx of partner, 2 branchpoint XY coords
            if ~isfield(filoInfo,'conIdx')
                filoInfo(idxSeedFilo).conIdx = idxBranch;
            else
                if isempty(filoInfo(idxSeedFilo).conIdx)
                    filoInfo(idxSeedFilo).conIdx = idxBranch;
                else
                    filoInfo(idxSeedFilo).conIdx(end+1,1) = idxBranch; % can have more than one attachment
                end
            end
            
            % get the link of the input
            branchPointX = inputPoints(E(idxFiloAttach(iFilo),2),1);
            branchPointY = inputPoints(E(idxFiloAttach(iFilo),2),2);
            %filoInfo(idxSeedFilo).conXYCoords = [branchPointX, branchPointY];
            % here it  might be advantageous to record the distance from the
            % base point of the filo (in pixels) NOTE might want to fit this to
            % get a better distance
            
            filoSeedPix = filoInfo(idxSeedFilo).Ext_pixIndicesBack;
            idxBranchPt = sub2ind(dims,branchPointY,branchPointX);
            distPix =  find(idxBranchPt==filoSeedPix);
            if isempty(distPix)
                distPix =NaN;
            end
            
            if ~isfield(filoInfo,'conXYCoords');
                filoInfo(idxSeedFilo).conXYCoords = [branchPointX,branchPointY,distPix];
            else
                if isempty(filoInfo(idxSeedFilo).conXYCoords)
                    filoInfo(idxSeedFilo).conXYCoords = [branchPointX,branchPointY,distPix];
                else
                    % first test length to make sure ok
                    
                    filoInfo(idxSeedFilo).conXYCoords(end+1,:) = [branchPointX,branchPointY,distPix];
                end
            end
            
            
            %structure the distance from base of attachment
            % get the local orientation of the seed filo (note could also get
            % this from the maxTh info..though there was some instability
            % here so decided just to do directly.
            testPoints = filoInfo(idxSeedFilo).Ext_coordsXY;
            testPoints = testPoints(~isnan(testPoints(:,1)),:); % make sure no NaNs
            % NaNs make the KDTree not work correctly all the time... better to just avoid.
            
            [idxBranchRegion,dist]= KDTreeBallQuery(testPoints,[branchPointX,branchPointY],3);
            if ~isempty(vertcat(idxBranchRegion{:})) &&  length(vertcat(idxBranchRegion{:}))>2
                idxBranchRegion = idxBranchRegion{:}(end-1:end,:);
                idxBranchRegion = sort(idxBranchRegion);  % sort so you know which one is nearer to the base of the filo
                
                seedFiloLocBranchRegXY = testPoints(idxBranchRegion,:); % last two points should be the furthest two from centerpoint
                vectSeedFiloLocBranchReg = [seedFiloLocBranchRegXY(2,1)-seedFiloLocBranchRegXY(1,1),seedFiloLocBranchRegXY(2,2)-seedFiloLocBranchRegXY(1,2)];
                magSeedVect = sqrt(vectSeedFiloLocBranchReg(1)^2+vectSeedFiloLocBranchReg(2)^2);
            else
                magSeedVect = NaN;
                vectSeedFiloLocBranchReg = [NaN NaN];
            end
            
            % walk the structure back
            EPNum = EPsCand(idxFiloAttach(iFilo));
            xEP = candFiloEPs{E(idxFiloAttach(iFilo),1)}(EPNum,1);
            yEP = candFiloEPs{E(idxFiloAttach(iFilo),1)}(EPNum,2);
            
            verticesEP = [yEP xEP];
            
            if labelCandCon(idxFiloAttach(iFilo))~=0 ; % something  is wrong ! skip
                testMask = zeros(size(out));
                % put the pixindices of the filo in the mask
                testMask(filoInfo(idxSeedFilo).Ext_pixIndicesBack)=1;
                testMask(vertcat(pixIdxCands{labelCandCon(idxFiloAttach(iFilo))})) = 1;
                testMask(pixGoodConnect{idxFiloAttach(iFilo)})=1;
                testMask = logical(testMask);
                
                testMask2 = bwmorph(testMask,'thin','inf'); % thin to avoid junctions that might be made upon the linkage
                
                %test for more than one junction if more
                % than one spur
                nn = padarrayXT(double(testMask2~=0), [1 1]);
                sumKernel = [1 1 1];
                nn = conv2(sumKernel, sumKernel', nn, 'valid');
                nn1 = (nn-1) .* (testMask2~=0);
                junctTest = nn1>2;
                CCJunct = bwconncomp(junctTest);
                
                if CCJunct.NumObjects >1
                    testMask = bwmorph(testMask2,'spur');
                    % add back the end point
                    testMask(yEP,xEP) = 1;
                elseif CCJunct.NumObjects == 1
                    testMask = testMask2; % else don't thin the test mask if it gets rid of all the junctions...
                end
                
                transform = bwdistgeodesic(testMask,xEP,yEP);
                % now need to get pixels and do fit on the branch
                pixIdxBack = nan(50,1); % overinitialize to make happy
                
                iPix = 1;
                while length(find(transform==iPix)) == 1
                    pixIdxBack(iPix) = find(transform==iPix); %
                    iPix = iPix +1;
                end
                pixIdxBack = pixIdxBack(~isnan(pixIdxBack));
                % sanity check here
                %         x = labelCandidates==labelCandCon(idxFiloAttach(iFilo));
                %         candNum = sum(x(:));
                %         linkNum = length(pixGoodConnect{idxFiloAttach(iFilo)});
                %         orig = filoInfo(idxFiloAttach).Ext_pixIndicesBack;
                %         test = candNum+linkNum+length(orig);
                
                % make sure to stop at branch point as sometime go through !
                idxStop =  find(pixIdxBack == idxBranchPt);
                if ~isempty(idxStop)
                    pixIdxBack = pixIdxBack(1:idxStop);
                end
                
                % get coords of the these pixels
                [yBack,xBack]= ind2sub(size(out),pixIdxBack);
                if length(xBack) >4
                    add = 4;
                else add = length(xBack)-1;
                end
                if add >0 % means only 1 pix should probably just prune these now as these are branche stubs and are likely just noise
                    % NOTE should problaby clean this up at some point
                    vectBranch = [xBack(end-add)-xBack(end),yBack(end-add) - yBack(end)];
                    magBranchVect = sqrt(vectBranch(1)^2+vectBranch(2)^2);
                else
                    vectBranch = [NaN,NaN];
                    magBranchVect = NaN;
                end
                edgePathCoord{1} = [yBack, xBack];
                % eventually replace
                % x = gcaProjectAndRecordFiloCoords(verticesEP,edgePathCoord,maxTh,maxRes,img,'numPixSearchForward',10);
                x = walkFiloForAndBack([], verticesEP,edgePathCoord,maxTh,maxRes,img,0,10);
                x.type = filoInfo(idxSeedFilo).type +1; % label type as subsidiary
                
                % make it such that this is kept empty the connection is
                % only stored once on the main branch.
                x.conIdx = [];
                x.conXYCoords = [];
                
                x.cross = 0;
                orientBranch = acosd(dot(vectSeedFiloLocBranchReg ,vectBranch)/magBranchVect/magSeedVect);
                x.orientation = orientBranch; % in degrees.
                x.localVectFilo = vectBranch;
                x.localVectAttach = vectSeedFiloLocBranchReg;
                
                fieldsx = fieldnames(x);
                fieldsFiloInfo = fieldnames(filoInfo);
                fieldsAddtoX = setdiff(fieldsFiloInfo,fieldsx);
                for i = 1:length(fieldsAddtoX)
                    x.(char(fieldsAddtoX(i))) = NaN; % no internal filopodia if branch
                end
                x.groupCount = filoInfo(idxSeedFilo).groupCount; % propogate this value forward
                % x.bodyType = NaN;
                x = orderfields(x); % put in alphabetical order
                filoInfo = orderfields(filoInfo);
                filoInfo(idxBranch) = x;
                % add coordinates to labelMatSeed so can find attachments in next
                % iteration
                
                candFiloAdded.Filo{iFilo} = pixIdxBack;
                clear x edgePathCoord pixIdxBack distPix;
            end % if  labelCanCon
        end % for iFilo
    else
        idxEndOnAttach = [];
    end % isempty
    
    attachMask1 = zeros(size(img));
    attachMask2 =zeros(size(img));
    attachMask3 = zeros(size(img));
    if ~isempty(idxFiloAttach)
        attachMask1(vertcat(candFiloAdded.Filo{:}))=1 ;
    end
    if ~isempty(idxBodyAttach)
        attachMask2(vertcat(candFiloAdded.Body{:})) = 1;
    end
    if ~ isempty(idxEndOnAttach) && length(idxEndOnAttach)>subtractEndOn;
        attachMask3(vertcat(candFiloAdded.EndOn{:})) = 1;
    end
    outputMasks.candFiloAdded.Body = attachMask2;
    outputMasks.candFiloAdded.Branch = attachMask1;
    outputMasks.candFiloAdded.EndOn = attachMask3;
    
    clear out links
    
    % check status
    if  sum(sum([outputMasks.candFiloAdded.Body(:)  ...
            outputMasks.candFiloAdded.Branch(:)   ...
            outputMasks.candFiloAdded.EndOn(:)])) == 0
        status = 0;
    else
        status = 1;
    end
else % no candidates were within the distance so exit and record zeros for output
    outputMasks.finalReconstruct = zeros(size(img));
    outputMasks.candFiloAdded.Body = zeros(size(img));
    outputMasks.candFiloAdded.Branch = zeros(size(img));
    outputMasks.links = zeros(size(img));
    status =0;
end % if sum

end % function 