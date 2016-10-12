
function [maskPostConnect,linkMask,status,TSFigs] = gcaConnectEmbeddedRidgeCandidates(internalCandEPs,internalSeedEPs,seedMask,labelMatRidgeCandEmbed,varargin)
%% gcaConnectEmbeddedRidgeCandidates

% SUMMARY: Connects Co-linear Candidate Internal Ridges (a filoTracker Filo Reconstruction Function)
% This function will take the internal filopodia endpoints and attempt to
% reconnect them to the user defined seed points (typically high-confidence external/internal filo)
% Candidates are filtered
% linking the two candidate ridges


%% INPUT:

% internalCanEPs: an rx2 array of candidate endpoints: these will be input
%                 points (ie you are searching for points in this population
%                 that meet the specified dist criteria ) (Post Cleaning)
%
%
% internalSeedEPs: (REQUIRED) : an rx2 array
%    marking the positions of the endpoints of high confidence
%    ridge signal (filopodia) outside of veilStemMask which will serve
%    as query pts for the KDTree Ball Query search
%    where r is the number of endpoints and 2 is the xy coords
%    of the endpoint. These ridge endpoints are typically pre-filterd such
%    that the endpoint of the ridge nearest the veilStem estimate is chosen,
%    as this is typically the only reasonable value of linking- though
%
%
% seedMask: (REQUIRED) : an rxc logical array (binary mask)
%    marking the full mask of the ridge candidates (Filopodia) post-cleaning
%    outside the veil that will serve as the seed for the internal filopodia linking
%
%
% labelMatRidgeCandEmbed: (REQUIRED) : an rxc double array
%    of the embedded ridge candidates (actin bundles) where each ridge
%    candidate connected component is given an independent numeric label
%    1:number of CCs
%
%% OPTIONAL : For plots
% img : (OPTIONAL) : an rxc double array of the original image
%    for plotting of troubleshoot overlays
%
% veilStem : (OPTIONAL) : an
%
%% PARAM
% 'maxRadiusLinkEmbedded' (PARAM) : Scalar
%    Only embedded ridge candidate end points that are within this max
%    search radius around each seed ridge endpoint are considered for matching.
%    Default: 10 Pixels
%
% 'geoThreshEmbedded' : (PARAM)  : Scalar
%    Hard threshold for the geometry of links: geometric linkages with
%    dot products of the vector of the ridge at the linking point and the
%    linear connection between the two ridges.
%    1 indicates perfect colinearity (0 degree angle), 0 indicates complete
%    orthoganality (90 angle angle), -1 indicates complete colinearity (-180
%    degree angle)
%    Default = 0.9 (~ 25 Degree Angle)
%
%% OUTPUT:
% maskPostConnect:
%   2D double matrix the size of the image, binary mask of all ridges
%   (ie filopdia) post connection
%
% linkMask:
%   2D double matrix the size of the image , binary mask of all viable links between ridge candidates made
%
% status: logical
%    1 if links were found, 0 if no viable links were found
%% INPUT Parser
ip = inputParser;

ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('internalCandEPs');
ip.addRequired('internalSeedEPs');
ip.addRequired('seedMask');
ip.addRequired('labelMatRidgeCandEmbed');


% Optional for plotting
ip.addOptional('img',[]);
ip.addOptional('edgeMask',[]);


ip.addParameter('maxRadiusLinkEmbedded',10,@(x) isscalar(x));
ip.addParameter('geoThreshEmbedded',0.9,@(x) isscalar(x));
ip.addParameter('TSOverlays',true);

% Plotting Parameters
ip.parse(internalCandEPs,internalSeedEPs,seedMask,labelMatRidgeCandEmbed,varargin{:});

p = ip.Results;

%% Initiate
imSize = size(seedMask);
countFigs = 1;
% extract the corresponding vectors
vectSeed = internalSeedEPs(:,3:4);
vectCand = internalCandEPs(:,3:4);
TSFigs = []; 
%%
% find all query points (internal filo EP coordinates) within a search radius surrounding the seed points (external filo EP coords)
[idx,d] = KDTreeBallQuery(internalCandEPs(:,1:2), internalSeedEPs(:,1:2), ip.Results.maxRadiusLinkEmbedded);

% Get Edges: tranfer data from cell to double:  column 1 = index of the seed while
% column 2 = the index of the candidate, rows =  n number of
% potential paths
% first col E = idx filo Cand, 2nd col = idx input point,
E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:length(internalSeedEPs(:,1)), 'UniformOutput', false);
E = vertcat(E{:});

%% sanity check
% plot all paths
%examplePlot =1 ;
%  if examplePlot ==1 ;
%candFiloNum = 1;% the number of the candidate filo you would like to plot
%EPlot = E(E(:,1)==candFiloNum,:);
% %
% %
% %
% %      pathidxPlot = find(E(:,1)==candFiloNum);
% %       c = colormap(lines(length(size(E,1))));
% %      scatter(candFiloEPs{candFiloNum}(:,1),candFiloEPs{candFiloNum}(:,2),'y','filled');
% %    scatter(inputPoints(EPlot(:,2),1),inputPoints(EPlot(:,2),2),'r','filled');
% %
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
%% Give each edge a cost based on geometry of linking the two points

if ~isempty(E) % continue
    
    idxCandAll = E(:,2);
    idxSeedAll = E(:,1);
    
    nEdges = size(idxCandAll,1);
    vectSeedPostKD =  vectSeed(idxSeedAll,:); % xy vect of the local direction at seed position
    vectCandPostKD = vectCand(idxCandAll,:);
    
    EPsCandPostKD = internalCandEPs(idxCandAll,1:2); % xy vect of all possible cand points
    EPsSeedPostKD = internalSeedEPs(idxSeedAll,1:2); % xy vect of all possible seed points
    
    
    %% small sanity
    % Seed To Cand
    
    %     if ip.Results.TSOverlays == true;
    %
    %         TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'on'); % reget the handle
    %         TSFigs(countFigs).name = 'Vectors';
    %
    %         if ~isempty(ip.Results.img);
    %             imshow(-ip.Results.img,[]);
    %         end
    %         hold on
    %
    %
    %         scatter(internalSeedEPs(E(:,1),1),internalSeedEPs(E(:,1),2),10,'r');
    %
    %         scatter(internalCandEPs(E(:,2),1),internalCandEPs(E(:,2),2),10,'b');
    %         spy(seedMask,'r');
    %
    %
    %         xSeed = internalSeedEPs(E(:,1),1);
    %         ySeed = internalSeedEPs(E(:,1),2);
    %
    %         % quiver(xSeed,ySeed,internalSeedsEPs(E,(:,1),, deltYConSeedToCand(:),'color','g');
    %
    %
    %         text(5,10,'Ridge Seeds','FontSize',10,'Color','r');
    %         text(5,20,'Embedded Candidates','FontSize',10,'Color','b');
    %
    %     end
    %%
    % vector from Cand to Seed assuming a linear connection for all edges
    % found by the KD tree
    deltXSeedToCand=  EPsCandPostKD(:,1) - EPsSeedPostKD(:,1);
    deltYSeedToCand = EPsCandPostKD(:,2) - EPsSeedPostKD(:,2);
    
    % norm factor
    dSeedToCand = sqrt((deltXSeedToCand').^2+ (deltYSeedToCand').^2)';
    
    % normalize to get connection vector for all edges (going from cand to
    % seed point
    vectSeedToCandNorm = ([deltXSeedToCand./dSeedToCand,deltYSeedToCand./dSeedToCand]);
    
    % get the dot product to test for colinearity of the local direction of the edge at the seed point
    % and direction of linear connection from seed to candidate
    dotProdSeedToCand = arrayfun(@(i) dot(vectSeedToCandNorm(i,:),vectSeedPostKD(i,:)),1:nEdges);
    costGeo1= dotProdSeedToCand';
    %%
    % vector from Cand to Seed assuming a linear connection for all edges
    % found by the KD tree
    deltXCandToSeed=  EPsSeedPostKD(:,1) - EPsCandPostKD(:,1);
    deltYCandToSeed = EPsSeedPostKD(:,2) - EPsCandPostKD(:,2);
    
    % norm factor
    dCandToSeed = sqrt((deltXCandToSeed').^2+ (deltYCandToSeed').^2)';
    
    % normalize to get connection vector for all edges (going from cand to
    % seed point
    vectCand2SeedNorm = ([deltXCandToSeed./dCandToSeed,deltYCandToSeed./dCandToSeed]);
    
    % get the dot product to test for colinearity between the local vector
    % of candidate and the linear link between candidate and seed.
    dotProdCandToSeed = arrayfun(@(i) dot(vectCand2SeedNorm(i,:),vectCandPostKD(i,:)),1:nEdges);
    
    costGeo2 = dotProdCandToSeed';
    
    costTotal = costGeo1 + costGeo2; % again this cost likely needs modification.
    
    %%
    if ip.Results.TSOverlays == true;
        
        TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'off');
        TSFigs(countFigs).name = 'Before_Matching_with_Vectors';
        TSFigs(countFigs).group = 'Reconstruct_Embedded';
        
        if ~isempty(ip.Results.img);
            imshow(-ip.Results.img,[]);
            hold on
            
            scatter(internalCandEPs(:,1),internalCandEPs(:,2),20,'b','filled'); % the endpoint to connect
            scatter(internalSeedEPs(:,1),internalSeedEPs(:,2),20,'r','filled');
            spy(labelMatRidgeCandEmbed>0,'b');
            spy(seedMask,'r')
            text(5,10,'Seeds','FontSize',10,'color','r');
            text(5,20,'Ridge Candidates Embedded','FontSize',10,'color','b');
            text(5,30,'With Local Geometry Vectors','FontSize',10,'color','k');
            text(5,40,'And Local Linkage Vectors','FontSize',10,'color','g');
            
            %
            quiver(internalCandEPs(:,1),internalCandEPs(:,2),...
                internalCandEPs(:,3),internalCandEPs(:,4),0.2,'color','b');
            quiver(internalSeedEPs(:,1),internalSeedEPs(:,2),internalSeedEPs(:,3),internalSeedEPs(:,4),...
                0.2,'color','r');
            
            arrayfun(@(x) plot([EPsCandPostKD(x,1),EPsSeedPostKD(x,1)],...
                [EPsCandPostKD(x,2),EPsSeedPostKD(x,2)],'color','g'),1:size(EPsCandPostKD,1));
            % plot the small vectors
            quiver(internalCandEPs(:,1),internalCandEPs(:,2),...
                internalCandEPs(:,3),internalCandEPs(:,4),0.2,'color','b');
            %
            quiver(EPsCandPostKD(:,1),EPsCandPostKD(:,2),vectCand2SeedNorm(:,1), vectCand2SeedNorm(:,2),...
                0.2,'color','g','linewidth',1);
            %
            quiver(EPsSeedPostKD(:,1),EPsSeedPostKD(:,2),vectSeedToCandNorm(:,1), vectSeedToCandNorm(:,2),...
                0.2,'color','g','linewidth',1);
            
            countFigs = countFigs+1;
        end
        
        %%
        
        %     vectConn = [deltX' deltY'];
        %     dConn = sqrt((deltX').^2+ (deltY').^2);
        %     costIntAndConn = arrayfun(@(i) dot(vectConn(i,:),vectInt(i,:))./dConn(i)./dInt(i),1:length(dConn));
        %     costSeedAndConn = arrayfun(@(i) dot(vectConn(i,:),vectSeed(i,:))./dConn(i)./dSeed(i),1:length(dConn));
        % costTotal =  abs(costIntAndConn) + abs(costSeedAndConn) + abs(costCandAndSeed);
        %costTotal  = abs(costIntAndConn) + abs(costSeedAndConn);
        % Sanity check
        %
        % could also use more of the local orientation instead of just the
        % ends.
        
        %% TSOverlays : Plot Linear Connections Color-Coded by Cost
        
        TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'off'); % reget the handle
        TSFigs(countFigs).name = 'KD_Results_By_Cost';
        TSFigs(countFigs).group = 'Reconstruct_Embedded';
        
        if ~isempty(ip.Results.img)
            
            imshow(-ip.Results.img,[]);
            hold on
        end
        
        if ~isempty(ip.Results.edgeMask);
            
            spy(ip.Results.edgeMask,'k');
            hold on
        end
        %Plot the full Candidate Pieces
        spy(labelMatRidgeCandEmbed,'k',5);
        hold on
        % Plot the full Seed Pieces
        spy(seedMask,'k',5);
        
        
        
        
        % Scatter the end points of the canidates considered (note only the
        % closest point  to the veil/stem estimate is currently considered
        
        scatter(internalCandEPs(:,1),internalCandEPs(:,2),5,'b','filled');
        scatter(internalSeedEPs(:,1),internalSeedEPs(:,2),5,'r','filled');
        cMapLength=128; cMap=jet(cMapLength);
        mapper=linspace(min(costTotal),max(costTotal),cMapLength)';
        
        % get closest colormap index for each feature
        D=createDistanceMatrix(costTotal,mapper);
        [sD,idxCMap]=sort(abs(D),2);
        
        for k=1:cMapLength
            idxCand = E(idxCMap(:,1) == k,2);
            idxSeed = E(idxCMap(:,1)==k,1);
            for iEdge = 1:length(idxCand) % some can have the same color
                plot([internalCandEPs(idxCand(iEdge),1),internalSeedEPs(idxSeed(iEdge),1)],...
                    [internalCandEPs(idxCand(iEdge),2),internalSeedEPs(idxSeed(iEdge),2)],'color',cMap(k,:),'Linewidth',2);
            end
        end
        text(5,5,'Color Paths By Cost ','FontSize',10,'Color','k');
        text(5,15,'Red : High : Stong Path' ,'FontSize',10,'Color','r');
        text(5,25,'Blue : Low : Weak Path', 'FontSize',10,'Color','b');
        text(5,35,['MaxRadius = ' num2str(ip.Results.maxRadiusLinkEmbedded) ' Pixels']);
        countFigs = countFigs+1;
        % c = colormap(lines(size(E,1)));
        % for i = 1:length(E(:,1))
        %     idxCand = E(i,2);
        %     idxSeed = E(i,1);
        %     plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)],'color',c(i,:));
        %     text(internalCandEPs(idxCand,1),internalCandEPs(idxCand,2),num2str(costCand1Path(idxCand),2),'color',c(i,:));
        %
        % end
    end % end ip.Results.TSOverlays
        %% Filter Based on a Geometric Threshold
       
        idxGood = (costGeo1>ip.Results.geoThreshEmbedded & costGeo2>ip.Results.geoThreshEmbedded);
        
        costTotal = costTotal(idxGood)';
        E = E(idxGood',:);
       
        %% TSOverlays: Histogram of Costs
        if ip.Results.TSOverlays == true
            TSFigs(countFigs).h = setAxis('off');
            TSFigs(countFigs).name = 'Cost Function Hist';
            TSFigs(countFigs).group = 'Reconstruct_Embedded'; 
            % plot the cost histogram
            subplot(2,1,1);
            hist(costTotal);
            ylabel('Number');
            xlabel('Cost');
            title('Cost Histogram After Geometric Filtering')
            
            % plot the cost scatter
            subplot(2,1,2);
            scatter(costGeo1(idxGood),costGeo2(idxGood),'b','filled');
            hold on
            scatter(costGeo1(~idxGood),costGeo2(~idxGood),'r','filled');
            legend('Connections Maintained (Blue) ','Connections Deleted (Red)' ...
                ,'Box','off','Location','BestOutside');
            
            
            line([ip.Results.geoThreshEmbedded,ip.Results.geoThreshEmbedded],[min(costGeo2),max(costGeo2)],'Linewidth',2,'color','r');
            line([min(costGeo1),max(costGeo1)],[ip.Results.geoThreshEmbedded,ip.Results.geoThreshEmbedded],'Linewidth',2,'color','r');
            xlabel('CostGeo1');
            ylabel('CostGeo2');
            countFigs = countFigs+1;
        end
        
        
        
        %% TSOverlays : Plot Linear Connections Color-Coded by Cost
        if ip.Results.TSOverlays == true;
            
             idxCMap = idxCMap(idxGood',:);
            
            TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'off'); % reget the handle
            TSFigs(countFigs).name = 'KD_Results_AFter_Filter_By_Geometry';
            TSFigs(countFigs).group = 'Reconstruct_Embedded'; 
            if ~isempty(ip.Results.img)
                
                imshow(-ip.Results.img,[]);
                hold on
            end
            
            if ~isempty(ip.Results.edgeMask);
                
                spy(ip.Results.edgeMask,'k');
            end
            %Plot the full Candidate Pieces
            spy(labelMatRidgeCandEmbed,'k',5);
            hold on
            % Plot the full Seed Pieces
            spy(seedMask,'k',5);
            
            
            
            
            % Scatter the end points of the canidates considered (note only the
            % closest point  to the veil/stem estimate is currently considered
            
            scatter(internalCandEPs(:,1),internalCandEPs(:,2),5,'b','filled');
            scatter(internalSeedEPs(:,1),internalSeedEPs(:,2),5,'r','filled');
            
            for k=1:cMapLength
                idxCand = E(idxCMap(:,1) == k,2);
                idxSeed = E(idxCMap(:,1)==k,1);
                for iEdge = 1:length(idxCand) % some can have the same color
                    plot([internalCandEPs(idxCand(iEdge),1),internalSeedEPs(idxSeed(iEdge),1)],[internalCandEPs(idxCand(iEdge),2),internalSeedEPs(idxSeed(iEdge),2)],'color',cMap(k,:),'Linewidth',1);
                end
            end
            text(5,5,'After Application of Geometry Threshold','FontSize',10);
            text(5,20,['of ' num2str(ip.Results.geoThreshEmbedded)],'FontSize',10);
            countFigs = countFigs+1;
            
            
        end % ip.Results.
        
        %% Use Graph Matching to choose quickly the best candidate among several competing edges
        
        [seedFiloNodes,~,nodeLabels] = unique(E(:,1),'stable'); % reason note: some of the filo will not be candidates as their endpoints are not within the given radius
        NNodeQuery = length(seedFiloNodes);
        [inputLinks,~,nodeLabelsInput] = unique(E(:,2),'stable'); % reason note: just in case two filo are competing over the same seed point
        nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
        EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
        numberNodes = length(inputLinks) + length(seedFiloNodes);
        
        
        
        M = maxWeightedMatching(numberNodes, EFinal, costTotal');
        E = E(M,:);
        paths=arrayfun(@(i) bresenham([internalSeedEPs(E(i,1),1) internalSeedEPs(E(i,1),2)], [internalCandEPs(E(i,2),1) internalCandEPs(E(i,2),2)]),...
            1:length(E(:,1)),'uniformoutput',0);
        linkMask = zeros(imSize);
        goodCands = zeros(imSize);
        % find labels of internal candidates to keep
        labels = arrayfun(@(i) labelMatRidgeCandEmbed(sub2ind(imSize,internalCandEPs(E(i,2),2),internalCandEPs(E(i,2),1))), 1:length(E(:,1)));
        % find the indexing of those labels
        idxCandKeep = arrayfun(@(i) find(labelMatRidgeCandEmbed == i),labels,'uniformoutput',0);
        goodCands(vertcat(idxCandKeep{:})) = 1;
        
        links = vertcat(paths{:});
        
        if ~isempty(links) % nothing that falls under this criteria
            % Add links to candidate mask
            idxLinks = sub2ind(imSize,links(:,2),links(:,1));
            
            linkMask(idxLinks) = 1;
            maskPostConnect = (linkMask|goodCands|seedMask);
            status = 1;
        else
            maskPostConnect = seedMask;
            status =0;
        end
        %% TS Figure : Show the final links
        if ip.Results.TSOverlays == true
            TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'off'); % reget the handle
            TSFigs(countFigs).name = 'Linking_Results';
            TSFigs(countFigs).group = 'Reconstruct_Embedded'; 
            
            if ~isempty(ip.Results.img)
                
                imshow(-ip.Results.img,[]);
                hold on
            end
            
            if ~isempty(ip.Results.edgeMask);
                
                spy(ip.Results.edgeMask,'k');
            end
            %Plot the full Candidate Pieces
            spy(labelMatRidgeCandEmbed,'k',5);
            hold on
            % Plot the full Seed Pieces
            spy(seedMask,'k',5);
            spy(linkMask,'y',5);
            
            text(10,10,'Yellow Marks Final Linking','FontSize',10,'Color','k');
            
            countFigs = countFigs+1;
            
        end % ip.Results.
        
        
        if ip.Results.TSOverlays == true 
            TSFigs(countFigs).h = setFigure(imSize(2),imSize(1),'off'); 
            TSFigs(countFigs).name = 'Final_Embedded';
                 
            
            
            if ~isempty(ip.Results.img)
                
                imshow(-ip.Results.img,[]);
                hold on
            end 
            
             if ~isempty(ip.Results.edgeMask);
                
                spy(ip.Results.edgeMask,'k');
            end
        
        spy(linkMask,'y'); 
        
        spy(goodCands,'b'); 
        
        spy(seedMask,'r'); 
        countFigs = countFigs +1; 
        
        end 
         
            
        
    else
        status = 0 ; % not links within range to even consider
        maskPostConnect = seedMask;
        linkMask = zeros(imSize);
    end
    
    
    
    
    
    
