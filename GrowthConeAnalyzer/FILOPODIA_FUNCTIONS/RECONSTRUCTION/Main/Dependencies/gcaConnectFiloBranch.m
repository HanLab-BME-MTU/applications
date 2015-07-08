function [outputMasks, filoInfo,status,pixIdxCandsUnMatched,candFiloEPsUnMatched,TSFigs] = gcaConnectFiloBranch(inputPoints, candFiloEPs, pixIdxCands, labelMatSeedFilo,filoInfo, maxRes,maxTh,img,normalsC,smoothedEdgeC,varargin)
% gcaConnectFilopodia:
% NOTES: used to be connectFilo until 20141022
% someday should consolidate some of these input parameters into structures
% I think right now it's annoying

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
%  
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
%
% maxRes : (maybe able to remove - no longer collecting the response for fitting. 
% 
% maxTh : (REQUIRED) 
%
%
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
%ip.addRequired('labelCandidates');
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
ip.parse(inputPoints,candFiloEPs,pixIdxCands,labelMatSeedFilo,filoInfo,maxRes,maxTh,img,normalsC,smoothedEdgeC,varargin{:});
p = ip.Results;

%% INITIATE
dims = size(labelMatSeedFilo); 
%dims = size(cellBoundary);
out = zeros(dims);
queryPoints = vertcat(candFiloEPs{:});
nq = size(queryPoints,1);
countFig = 1;
% replaces line 122
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

                            
%% 

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
    %
    E = [E idxEnd];
    
    %Collect all distance values
    D = arrayfun(@(i) vertcat(distNode{i}{:}),1:numel(distNode),'uniformoutput',0);
    D = vertcat(D{:});
    
    % calc the vector of each candidate and it's length for the orientation
    % part of the costMat
    % vectCand =  cellfun(@(x) [x(1,1)-x(2,1), x(1,2) - x(2,2)], candFiloEPs ,'uniformoutput',0);
    %dCand = cellfun(@(x) sqrt((x(1,1)-x(2,1))^2 + (x(1,2)-x(2,2))^2),candFiloEPs,'uniformoutput',0);
    
    % calculate the orientation of each candidate filo
    seedPtsx = inputPoints(E(:,2),1);
    seedPtsy = inputPoints(E(:,2),2);
    labelsSeed = labelMatSeedFilo(sub2ind(dims,seedPtsy,seedPtsx));
    % make sure input seed and labels are the same as now you did some crazy stuff
    % with teh
    
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
        spy(labelCandidates,'r');
        text(5,10,'Seed','FontSize',10,'color','r');
        text(5,25,'Candidate Ridges','FontSize',10,'color','b');
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
    
    
    %%
    
    
    
    
    
    % covert indexing to node labels to input into graph matching
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
    
    %% TS Overlay : Plot
   % figure('Visible','on');
    %img = img./((2^16)-1); % ask ludo the bitdepth of the camera should do in beginning
    plotPaths =0;
%     c = colormap(lines(size(E,1)));
%     if plotPaths ==1
%         imshow(img,[])
%         hold on
%         spy(labelMatSeedFilo,'b');
%         spy(labelCandidates,'y');
%         
%     end
    % remember 1st Col E = idx filo Cand, 2nd col = idx input point, and 3rd
    % col = sideOfCandFilo Idx
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
        if plotPaths ==1
            plot([xCoordQuery;xCoordInput],[yCoordQuery;yCoordInput],'color',c(iPath,:));
        end
        dotProd(iPath,1) = dot(vectorCand,vectTest)/dTest;
    end
    
    % plot for paper
    % if examplePlot ==1 ;
    %     candFiloNum = 1;% the number of the candidate filo you would like to plot
    %     EPlot = E(E(:,1)==candFiloNum,:);
    %
    %
    %
    %     pathidxPlot = find(E(:,1)==candFiloNum);
    %      c = colormap(lines(length(size(E,1))));
    %     scatter(candFiloEPs{candFiloNum}(:,1),candFiloEPs{candFiloNum}(:,2),'y','filled');
    %     scatter(inputPoints(EPlot(:,2),1),inputPoints(EPlot(:,2),2),'r','filled');
    %
    %
    %
    %     for iPath =  1:length(pathidxPlot)
    %        xCoordQuery = candFiloEPs{E(iPath,1)}(E(iPath,3),1);
    %               yCoordQuery = candFiloEPs{E(iPath,1)}(E(iPath,3),2) ;
    %               xCoordInput = inputPoints(E(iPath,2),1);
    %               yCoordInput = inputPoints(E(iPath,2),2);
    %               x = [xCoordQuery;xCoordInput]';
    %               y = [yCoordQuery;yCoordInput]';
    %               plot(x,y,'color',c(iPath,:))
    %     end
    % end
    %% COST Definition
    % Distance Term: make so that a lower distance = positive weight
    D = max(D)-D;
    % normalize
    D= D./max(D);
    
    %% Color Code by distance term 
    

    
    %% Color Code by Intensity Term 
    
    
    
    % Intensity Term: defined as the average intensity from end point to body attachment: higher intensity = higher weights
    % OLD Normaliztion Method : based on normalization by intensity of per
    % pixel max 
    %intMax = max(vertcat(int{:}));
    %int = cellfun(@(x) x./intMax,int,'uniformoutput',0); % normalize by max intensity
    
    % should probably do per path max (20150708) 
    % get the average intensity. 
    int = cellfun(@(x) mean(x), int); % higher this intensity the more likely to like
    maxIntPath = max(int);  
    minIntPath = min(int); 
    normInt = (int - minIntPath)./(maxIntPath-minIntPath); 
    normInt = normInt'; 
    
   
    
  
    
    % Orientation Term:  above 1 favored 0 unfavored
    costTotal = (0.5*D + normInt +dotProd+dotCandAndSeed); % note to self: let's see what we get by just combining these values
    % linearly at first: in the end might want to disfavor the distance term
    % and favor more the int and dotProd term.
    
   
    
    
    
%     setAxis
%     hist(costTotal)
    %% Plot all the paths by the costTotal
    % %% FIRST DO INDIVIDUAL COMPONENTS 
%     if ip.Results.TSOverlays == true 
%         TSFigs(countFig).h = setFigure(dims(1),dims(2),'on'); 
%         TSFigs(countFig).name = 'Histograms_of_Cost_Parameters_BeforeGeoThresh'; 
%         TSFigs(countFig).group = 'Resconstruct_FiloBranch'; 
%         
%         setAxis('on'); 
%         subplot(5,1,1); 
%         hist(D,100); 
%         
%         subplot(5,1,2); 
%         hist(normInt,100); 
%         
%         subplot(5,1,3); 
%         hist(dotProd,100); 
%         
%         
%         subplot(5,1,4); 
%         hist(dotCandAndSeed,100); 
%         
%         subplot(5,1,5); 
%         hist(costTotal);  
%     end 
       
        
        
      
    
    
    
    
    
     %% Color Code by distance term 
     
     if ip.Results.TSOverlays == true 
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on'); 
        TSFigs(countFig).name = 'Potential_Paths_with_Cost_D'; 
        TSFigs(countFig).group = 'Reconstruct_FiloBranch'; 
    
      
        
        
        
    % % Start Plot
        imshow(-img,[]); 
        hold on 
        spy(labelMatSeedFilo>0,'k'); 
        spy(labelCandidates>0,'k',5); 
        scatter(seedPtsx(:),seedPtsy(:),'k','filled'); 
        allCandEPs = vertcat(candFiloEPs{:}); 
        scatter(allCandEPs(:,1),allCandEPs(:,2),'k','filled'); 
        
        % create distance mapper
        cMapLength=128; cMap=jet(cMapLength);
        mapper=linspace(min(D),max(D),cMapLength)';
        
        DMap=createDistanceMatrix(costTotal,mapper);
        [sD,idxCMap]=sort(abs(DMap),2);
        
        for k = 1:length(cMap);
            if sum(idxCMap(:,1)==k)~=0
                toPlot = iSeg(idxCMap(:,1) == k);
                
                cellfun(@(x) plot([x(1,1),x(end,1)],[x(1,2),x(end,2)],'color',cMap(k,:)),toPlot);
                clear toPlot
            else
            end
            % make colormap of costs.
            
            % show each segment in iSeg cell plotted by the costTotal color
            %
        end 
     
    countFig = countFig+1; 
     end
 %%   
    
    
    
    
 %% TOTAL    
    if ip.Results.TSOverlays == true 
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on'); 
        TSFigs(countFig).name = 'Potential_Paths_with_Cost'; 
        TSFigs(countFig).group = 'Reconstruct_FiloBranch'; 
        
        % Start Plot
        imshow(-img,[]); 
        hold on 
        spy(labelMatSeedFilo>0,'k'); 
        spy(labelCandidates>0,'k',5); 
        scatter(seedPtsx(:),seedPtsy(:),'k','filled'); 
        allCandEPs = vertcat(candFiloEPs{:}); 
        scatter(allCandEPs(:,1),allCandEPs(:,2),'k','filled'); 
        
        % create distance mapper
        cMapLength=128; cMap=jet(cMapLength);
        mapper=linspace(min(costTotal),max(costTotal),cMapLength)';
        
        DMap=createDistanceMatrix(costTotal,mapper);
        [sD,idxCMap]=sort(abs(DMap),2);
        
        for k = 1:length(cMap);
            if sum(idxCMap(:,1)==k)~=0
                toPlot = iSeg(idxCMap(:,1) == k);
                
                cellfun(@(x) plot([x(1,1),x(end,1)],[x(1,2),x(end,2)],'color',cMap(k,:)),toPlot);
                clear toPlot
            else
            end
            % make colormap of costs.
            
            % show each segment in iSeg cell plotted by the costTotal color
            %
        end 
    
    countFig = countFig+1; 
    end % if TSOverlays
    
    
 
    %% Perform Matching to Resolve Graph
    E = [E costTotal D normInt dotProd dotCandAndSeed];
    
    E = E(dotProd>ip.Results.geoThreshFiloBranch,:); 
    EFinal = EFinal(dotProd>ip.Results.geoThreshFiloBranch,:);
    costTotal = costTotal(dotProd>ip.Results.geoThreshFiloBranch);
    idxCMap =  idxCMap(dotProd>ip.Results.geoThreshFiloBranch,:);
    iSeg= iSeg(:,dotProd>ip.Results.geoThreshFiloBranch);
     %numberNodes =  sum(dotProd); 
    
    %numberNodes = size(EFinal,1); 
    [candFiloNodes,~,nodeLabels] = unique(EFinal(:,1),'stable'); % reason note: some of the filo will not be candidates as their endpoints are not within the given radius
    NNodeQuery = length(candFiloNodes);
    [inputLinks,~,nodeLabelsInput] = unique(EFinal(:,2),'stable'); % reason note: just in case two filo are competing over the same seed point
    nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
    EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
    numberNodes = length(inputLinks) + length(candFiloNodes);
    %numberOfNodes = 98; 
    
    
    %% TS Overlay : Geometry Thresholds
    if ip.Results.TSOverlays == true 
       TSFigs(countFig).h = setFigure(dims(1),dims(2),'on'); 
       TSFigs(countFig).name = 'Histograms_of_Cost_Parameters_AfterGeoThresh'; 
       TSFigs(countFig).group = 'Reconstruct_FiloBranch'; 
        setAxis('on'); 
        subplot(5,1,1); 
        [n,center] = hist(E(:,4),100); 
        bar(center,n/max(n)); 
        xlabel('Cost Total'); 
        axis([-1,3.5,0,1]); 
        
        
        subplot(5,1,2); 
        [n,center] = hist(E(:,5),100); 
        bar(center,n/max(n)); 
        xlabel('Distance'); 
        axis([0,1,0,1]); 
        
        subplot(5,1,3); 
        [n,center] = hist(E(:,6),100); 
        bar(center,n/max(n)); 
        xlabel('Mean Intensity'); 
        axis([0,1,0,1]); 
        
        subplot(5,1,4); 
        [n,center] = hist(E(:,7),100);
        bar(center,n/max(n)); 
        xlabel('Geometry With Linker'); 
        axis([-1,1,0,1]); 
        
        subplot(5,1,5); 
        [n,center] =  hist(E(:,8),100);  
        bar(center,n/max(n)); 
        xlabel('Geometry Candidate and Seed'); 
        axis([0,1,0,1]); 
    end 
    
    
    if ip.Results.TSOverlays == true
        TSFigs(countFig).h = setFigure(dims(1),dims(2),'on');
        TSFigs(countFig).name = 'Potential_Paths_After_Geometry_Threshold';
        TSFigs(countFig).group = 'Reconstruct_FiloBranch'; 
     
        imshow(-img,[]);
        
        hold on
           text(5,10,['Geometry Threshold' num2str(ip.Results.geoThreshFiloBranch)],'FontSize',10,'Color','k'); 
        
        spy(labelMatSeedFilo>0,'k');
        spy(labelCandidates>0,'k',5);
        scatter(seedPtsx(:),seedPtsy(:),'k','filled');
        allCandEPs = vertcat(candFiloEPs{:});
        scatter(allCandEPs(:,1),allCandEPs(:,2),'k','filled');
        % create distance mapper
        cMapLength=128; cMap=jet(cMapLength);
%         mapper=linspace(min(costTotal),max(costTotal),cMapLength)';
%         
%         D=createDistanceMatrix(costTotal,mapper);
%         [sD,idxCMap]=sort(abs(D),2);
        
        for k = 1:length(cMap);
            if sum(idxCMap(:,1)==k)~=0
                toPlot = iSeg(idxCMap(:,1) == k);
                
                cellfun(@(x) plot([x(1,1),x(end,1)],[x(1,2),x(end,2)],'color',cMap(k,:)),toPlot);
                clear toPlot
            else
            end
            % make colormap of costs.
            
            % show each segment in iSeg cell plotted by the costTotal color
            %
        end
        
    end % ip.Results.TSOverlays
    %% Perform Weighted Graph Matching
    
    
    
    if ~isempty(EFinal)
        %     idx = E(:,1) < E(:,2);
        
        %     E = E(idx,:); % remove redundancy
        %% NOTE Currently not counting the number of nodes correctly! 
        % might be because I never considered the seed coords might have
        % the same labels. 
        %% TESTING 
        %numberNodes = size(EFinal,1)-1;  
        %% 
        M = maxWeightedMatching(numberNodes, EFinal, costTotal);
        % check for double labels
        % convertBack
        % E = [candFiloNodes(nodeLabels) inputLinks(nodeLabelsInput)]; % convert back to original indices of input and query points
        E = E(M,:);% get those edges that matched (from original indexing)
    end   % isempty
    
    
    
    %
    if ~isempty(E) % might be empty now if all were repeats
        pixIdxCandsUnMatched(E(:,1)) = []; % take out the matched candidates 
        candFiloEPsUnMatched(E(:,1)) =[]; 
        % add linear segments corresponding to linked endpoints
        % actually this gets me into trouble if really want to link
        % these effectively need to consider the fluorescence intensity
        % into the cost of linking each point this can throw off my fits
        % it's a little stupid because we have some of this junction info
        % before the NMS but I throw it away
        goodConnect = iSeg(M);
        
%% Troubleshoot overlay 

        
        
        
        
        
        % convert to pixIdx
        pixGoodConnect = cellfun(@(i) sub2ind(dims,i(:,2),i(:,1)), goodConnect,'uniformoutput',0);
        out(vertcat(pixGoodConnect{:}))= 1;
        
    end
   %% If only keep smooth continuations.  
    
    
    
    %% Plot final results.
    
    
    links = out;
    outputMasks.links = links;
    allInputMask= (labelCandidates>0|labelMatSeedFilo>0); 
    out = double(out|allInputMask); 
    %out = double(out | cellBoundary);
    outputMasks.finalReconstruct = out;
    
    %% start documenting data
    labelInputCon = zeros(length(E(:,1)),1);
    % get the filoIdx of those filo to which attachments have been made
    for iMatch = 1:length(E(:,1))
        labelInputCon(iMatch) = labelMatSeedFilo(sub2ind(dims,inputPoints(E(iMatch,2),2),inputPoints(E(iMatch,2),1)));
        %% Need to Fix 20150604 Change to a cell array
        % NOTE 201500704 this was previously a bit redundant - your
        % info as to the label of the canidate was already in E(iMatch,1)
        %labelCandCon(iMatch) = labelCandidates(sub2ind(dims,candFiloEPs{E(iMatch,1)}(1,2),candFiloEPs{E(iMatch,1)}(1,1)));
        % new added 20150704
        
        labelCandCon(iMatch) = E(iMatch,1); 
%         
%         
%         %%
%         % NOte 2015070704 what? why would I need to switch this?
         if E(iMatch,3) == 1;
             EPsCand(iMatch) = 2 ;
         else
             EPsCand(iMatch) = 1;
         end
    end
    
    
    
    
    % other end of candidate
    
    
    
    % neurite edge mask is labeled with 1 (thick parts) and  2 (thin parts), so
    % indexing is off by 2 : FIX
    labelInputCon = labelInputCon -1;% change back to 1
    
 %% Document Veil Stem Attachments    
    % some might be attached to mask some will make new branch structures
    %     idxBodyAttachThick = find(labelInputCon ==-1); % idx of those matches that are attached to the body
    %     idxBodyAttachThin = find(labelInputCon ==0);
    idxBodyAttach = find(labelInputCon==0);
    
    %for iBody = 1:2
    %         if iBody == 1;
    %             idxBodyAttach = idxBodyAttachThick;
    %
    %         else
    %             idxBodyAttach = idxBodyAttachThin;
    %         end
    
    
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
                
               % The Candidate Filo 
                % To replace 
               %testMask(labelCandidates == labelCandCon(idxBodyAttach(iFilo)))=1;
               
               % Fix 20150704
               testMask(vertcat(pixIdxCands{labelCandCon(idxBodyAttach(iFilo))})) = 1; 
                
               % The Connection 
                testMask(pixGoodConnect{idxBodyAttach(iFilo)}) = 1;
                
                testMask = bwmorph(testMask,'thin','inf'); % added 20141026 NEED TO THIN
                testMask = logical(testMask);
                transform = bwdistgeodesic(testMask,xEP,yEP);
                % now need to get pixels and do fit on the branch
                pixIdxBack = nan(50,1); % overinitialize to make happy
                %endpoint of candidate coord find coord NOT part of connection
                
                
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
                        avgNormLocal = mean(normalsC(idx,:),1);% might want to change to a majority? (i don't think these vectors are really normalized)
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
                             angleToBody = rad2deg(acos(cosAngle));
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
                    % x = fitLinescansNew(x,dims,0,11,0,0); % add the new filo fit
                    
                    idxBodyAttachFiloInfo = numel(filoInfo) +1;
                    
                    %Fix fieldnames to match (might not be a problem in later versions
                    %if initialize data structure more appropriately: we'll see if it's
                    %worth fixing
                    fieldsx = fieldnames(x);
                    fieldsFiloInfo = fieldnames(filoInfo);
                    fieldsAddtoX = setdiff(fieldsFiloInfo,fieldsx);
                    
                    %         for i = 1:length(fieldsAddtoFiloInfo)
                    %             filoInfo(idxBodyAttachFiloInfo).(char(fieldsAddtoFiloInfo(i))) = [];
                    %         end
                    for i = 1:length(fieldsAddtoX)
                        x.(char(fieldsAddtoX(i))) = NaN; % In future actually should look for internal filo here.
                    end
                    
                    x = orderfields(x); % put in alphabetical order
                    filoInfo = orderfields(filoInfo);
                    
                    filoInfo(idxBodyAttachFiloInfo) = x;
                    
                    candFiloAdded.Body{iFilo} = pixIdxBack;
                    clear x edgePathCoord xBack yBack
                    %
                    
                end % don't record if isempty (this sometimes happens if there are weird cycles in the response % 20140426
            end % if labelCandCon
        end % iFilo
    end % ~isempty
    % end % for iBody
%% Document Filopodia Attachments     
    idxFiloAttach = find(labelInputCon >0); % idx of those matches that are attached to a seed filo
    %
    
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
            if labelCandCon(idxEndOnAttach(iEndon)) ~=0 % Again quick fix for the problem: NOTE 20141017 this 'quick fix is causing some problems
                % of its own.
                %testMask(labelCandidates == labelCandCon(idxEndOnAttach(iEndon)))=1; % get the candidate pixels
                % Fixed 20150704
                testMask(vertcat(pixIdxCands{labelCandCon(idxEndOnAttach(iEndon))})) = 1; 
                
                testMask(pixGoodConnect{idxEndOnAttach(iEndon)}) = 1; % get the links
                testMask(filoInfo(labelsEndon(iEndon)).Ext_pixIndicesBack) =1;
                % testMask(filoInfo(labelsEndon(iEndon)).Ext_pixIndicesFor(1,1)) =1; % 20141017 : MAINLY DUE TO THIS PIECE HERE I THINK
                % THIS WAS NOT THE ACTUAL PIX USED FOR THE ENDPOINT..just the
                % pixIndicesBack should techically be sufficient.
                % if this connection introduces points of intersection it is
                % not cool
                
                close gcf
                
                nn = padarrayXT(double(testMask~=0), [1 1]);
                sumKernel = [1 1 1];
                nn = conv2(sumKernel, sumKernel', nn, 'valid');
                nn1 = (nn-1) .* (testMask~=0);
                junctMask = nn1>2;
                if sum(junctMask(:)) == 0;
                    
                    % proceed
                    
                    
                   % testMask(labelMatSeedFilo==1|labelMatSeedFilo==2) = 1; % get neuriteBody
                   testMask(labelMatSeedFilo==1) = 1; % note changed back 20150704 because no longer document 'thin body'
                    %         testMask(labelCandidates == labelCandCon(idxEndOnAttach(iEndon)))=1; % get the candidate pixels
                    %         testMask(pixGoodConnect{idxEndOnAttach(iEndon)}) = 1; % get the links
                    %         testMask(filoInfo(labelsEndon(iEndon)).Ext_pixIndicesBack) =1;
                    %         testMask(filoInfo(labelsEndon(iEndon)).Ext_pixIndicesFor(1,1)) =1;
                    % make sure to thin!! % might not be necessary if don't add noisy stuff
                    %
                    
                    testMask = logical(testMask);
                    testMask = bwmorph(testMask,'thin','inf');
                    % test to make sure it doesn't still introduce a connectivity probl
                    transform = bwdistgeodesic(testMask,xEP,yEP);
                    % now need to get pixels and do fit on the branch
                    pixIdxBack = nan(50,1); % overinitialize to make happy
                    %endpoint of candidate coord find coord NOT part of connection
                    
                    
                    iPix = 1;
                    while length(find(transform==iPix)) == 1
                        pixIdxBack(iPix) = find(transform==iPix); %
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
                    if isfield(filoInfo,'conIdx')
                        x.conIdx = []; % just to make sure fields are consistent
                        x.conXYCoords= [];
                    end
                    % refit
                    %x = fitLinescansNew(x,dims,0,11,0,0); % add the new filo fit
                    % make sure to take off the flag that doesn't allow you to re-do
                    
                    % keep any old field the same (ie internal filo etc)
                    fieldsx = fieldnames(x);
                    fieldsFiloInfo = fieldnames(filoInfo);
                    fieldsAddtoX = setdiff(fieldsFiloInfo,fieldsx);
                    for iField = 1:length(fieldsAddtoX)
                        x.(char(fieldsAddtoX(iField))) = filoInfo(labelsEndon(iEndon)).(char(fieldsAddtoX(iField)));
                    end
                    
                    
                    
                    
                    x = orderfields(x); % put in alphabetical order
                    filoInfo = orderfields(filoInfo);
                    
                    filoInfo(labelsEndon(iEndon)) = x; % rewritethe data
                    pixAdd =  pixIdxCands{labelCandCon((idxEndOnAttach(iEndon)))}; % fixed 20150704
                   % pixAdd  = find(labelCandidates == labelCandCon(idxEndOnAttach(iEndon)));
                   % ok = isequal(pixAddTest,pixAdd); 
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
                        display('saving pixels')
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
                
            else
            end
            
            
            idxBranch = numel(filoInfo) +1;
            % filoInfo(idxBranch).type =5; % record the new filo marking it as a branch subsidiary
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
            filoInfo(idxSeedFilo).conXYCoords = [branchPointX, branchPointY, distPix]; % for now save in same part of
            %structure the distance from base of attachment
            % get the local orientation of the seed filo (note could also get
            % this from the maxTh info..
            testPoints = filoInfo(idxSeedFilo).Ext_coordsXY;
            testPoints = testPoints(~isnan(testPoints(:,1)),:); % make sure no NaNs added 20141201 for instance at the border...
            % NaNs make the KDTree not work correctly all the time (check out
            % why - but I noticed at smaller radii it will not find the points
            % you are searching for... better to just avoid.
            
            [idxBranchRegion,dist]= KDTreeBallQuery(testPoints,[branchPointX,branchPointY],3);
            if ~isempty(vertcat(idxBranchRegion{:})) &&  length(vertcat(idxBranchRegion{:}))>2 % added if statement 08-25-2013 and 09-01
                idxBranchRegion = idxBranchRegion{:}(end-1:end,:); % 09-01-2013 why not using more points here for vector?
                idxBranchRegion = sort(idxBranchRegion);  % sort so you know which one is nearer to the base of the filo
                
                seedFiloLocBranchRegXY = testPoints(idxBranchRegion,:); % last two points should be the furthest two from centerpoint
                vectSeedFiloLocBranchReg = [seedFiloLocBranchRegXY(2,1)-seedFiloLocBranchRegXY(1,1),seedFiloLocBranchRegXY(2,2)-seedFiloLocBranchRegXY(1,2)];
                magSeedVect = sqrt(vectSeedFiloLocBranchReg(1)^2+vectSeedFiloLocBranchReg(2)^2);
            else
                magSeedVect = NaN;
                vectSeedFiloLocBranchReg = [NaN NaN]; % added 201411201
                
            end
            
            
            
            
            % walk the structure back
            EPNum = EPsCand(idxFiloAttach(iFilo));
            xEP = candFiloEPs{E(idxFiloAttach(iFilo),1)}(EPNum,1);
            yEP = candFiloEPs{E(idxFiloAttach(iFilo),1)}(EPNum,2);
            
            verticesEP = [yEP xEP];
            %out = logical(out); NOTE don't want to feed in whole mask here as
            %you can have overlaps with other filo. might be good for documenting
            % cross overs eventually but for now we want ALL the coordinates
            % going back to the site of attachment
            % having a problem here (rhoA_02 KD new day frame 23...for some
            % reason one of the labelCandCon(idxFiloAttach(iFilo)) is equal to
            % 0 !! so it is not a label.. TROUBLESHOOT LATER just need to get
            % orientation info now for meeting 06-15-2013
            if labelCandCon(idxFiloAttach(iFilo))~=0 ; % something  is wrong ! skip
                testMask = zeros(size(out));
                % put the pixindices of the filo in the mask
                testMask(filoInfo(idxSeedFilo).Ext_pixIndicesBack)=1;
                %testMask(labelCandidates == labelCandCon(idxFiloAttach(iFilo))) = 1;
               % fixed 20150704
                testMask(vertcat(pixIdxCands{labelCandCon(idxFiloAttach(iFilo))})) = 1; 
                
                
                testMask(pixGoodConnect{idxFiloAttach(iFilo)})=1;
                testMask = logical(testMask);
                testMask = bwmorph(testMask,'thin','inf'); % thin to avoid junctions that might be made upon the linkage
                
                % attempt to fix  20141105 -test for more than one junction if more
                % than one spur
                nn = padarrayXT(double(testMask~=0), [1 1]);
                sumKernel = [1 1 1];
                nn = conv2(sumKernel, sumKernel', nn, 'valid');
                nn1 = (nn-1) .* (testMask~=0);
                junctTest = nn1>2;
                CCJunct = bwconncomp(junctTest);
                if CCJunct.NumObjects >1
                    testMask = bwmorph(testMask,'spur');
                    % add back the end point
                    testMask(yEP,xEP) = 1;
                    
                end
                
                % MARIA COME BACKE TO THIS! sometimes here still have small junctions that are
                % introduced by the linker ... One way to avoid this is to perform
                % a spur operation if have more than one junction introduced.
                % though that is not very elegant :( - the problem is simply I
                % record pixels back until the junction point. If there are two
                % junctions it will truncate the recording prematurely. -
                % other thing you can do is test for junctions in the linker - filo
                % combo.. again not elegant. -- above is one solution 20141105- it
                % worked for RhoA 01 20130827 cell 01 filopodia 40 where this was a
                % problem
                transform = bwdistgeodesic(testMask,xEP,yEP);
                % now need to get pixels and do fit on the branch
                pixIdxBack = nan(50,1); % overinitialize to make happy
                %endpoint of candidate coord find coord NOT part of connection
                
                
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
                x.type = filoInfo(idxSeedFilo).type +1; % label type as subsidiary (might want to change this)
                x.conIdx = idxFiloAttach(iFilo);
                x.conXYCoords =  [branchPointX branchPointY distPix];
                
                
                x.cross = 0;
                orientBranch = rad2deg(acos(dot( vectSeedFiloLocBranchReg ,vectBranch)/magBranchVect/magSeedVect));
                x.orientation = orientBranch; % in degrees.
                x.localVectFilo = vectBranch;
                x.localVectAttach = vectSeedFiloLocBranchReg;
                %         x = fitLinescansNew(x,dims,0,11,0,0); % add the new filo fit
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
            end % if  labelCanCon (check will likely remove later)
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


end


%
%
%





