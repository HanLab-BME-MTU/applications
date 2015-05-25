function [maskPostConnect,linkMask,status,TSFigs] = gcaConnectEmbeddedRidgeCandidates(internalCandEPs,internalSeedEPs,seedMask,labelMatRidgeCandEmbed,vectSeed,vectInt,dSeed,dInt,varargin)
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
%% PARAMS: 
% img : (OPTIONAL) : an rxc double array of the original image 
%    for plotting of troubleshoot overlays
% 
% veilStem : (OPTIONAL) : an 
% 
% 
% 'maxRadiusLinkEmbedded' (PARAM) : Scalar 
%    Only embedded ridge candidate end points that are within this max 
%    search radius around each seed ridge endpoint are considered for matching.
%    Default: 10 Pixels
%
% add linearity filter... currently 0.7 see line 189
%
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
ip.addRequired('vectSeed'); 
ip.addRequired('vectInt'); 
ip.addRequired('dSeed'); 
ip.addRequired('dInt'); 

% Optional for plotting
ip.addOptional('img',[]); 
ip.addOptional('edgeMask',[]); 


ip.addParameter('maxRadiusLinkEmbedded',10,@(x) isscalar(x));
ip.addParameter('TSOverlays',true); 

% Plotting Parameters 
ip.parse(internalCandEPs,internalSeedEPs,seedMask,labelMatRidgeCandEmbed,vectSeed,vectInt,dSeed,dInt,varargin{:});

p = ip.Results;

%% Initiate 
imSize = size(seedMask); 
countFigs = 1; 
%%
% find all query points (internal filo EP coordinates) within a search radius surrounding the seed points (external filo EP coords)
[idx,d] = KDTreeBallQuery(internalCandEPs, internalSeedEPs, ip.Results.maxRadiusLinkEmbedded);

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
    
    vectInt = vertcat(vectInt{idxCandAll});
    dInt = vertcat(dInt{idxCandAll});
    vectSeed = vertcat(vectSeed{idxSeedAll});
    dSeed = vertcat(dSeed{idxSeedAll});
    
    
    
   % costCandAndSeed = arrayfun(@(i) dot(vectInt(i,:),vectSeed(i,:))./dInt(i)./dSeed(i),1:length(dInt));
    
    EPidxSeed = sub2ind(imSize,internalSeedEPs(:,2),internalSeedEPs(:,1));
    EPidxCand = sub2ind(imSize,internalCandEPs(:,2),internalCandEPs(:,1));
    
    
    
    % need to likewise need to maintain linearity
    % among the connection and the two pieces
    % get orientation of the connecting piece
    deltX = arrayfun(@(i) internalSeedEPs(E(i,1),1)-internalCandEPs(E(i,2),1),1:length(E(:,1)));
    deltY= arrayfun(@(i) internalSeedEPs(E(i,1),2)-internalCandEPs(E(i,2),2),1:length(E(:,1)));
    vectConn = [deltX' deltY'];
    dConn = sqrt((deltX').^2+ (deltY').^2);
    costIntAndConn = arrayfun(@(i) dot(vectConn(i,:),vectInt(i,:))./dConn(i)./dInt(i),1:length(dConn));
    costSeedAndConn = arrayfun(@(i) dot(vectConn(i,:),vectSeed(i,:))./dConn(i)./dSeed(i),1:length(dConn));
    % costTotal =  abs(costIntAndConn) + abs(costSeedAndConn) + abs(costCandAndSeed);
    costTotal  = abs(costIntAndConn) + abs(costSeedAndConn); 
    % Sanity check
    %
    % could also use more of the local orientation instead of just the
    % ends.
   
%% TSOverlays : Plot Linear Connections Color-Coded by Cost 
if ip.Results.TSOverlays == true;
  
   
        
        TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'on'); % reget the handle 
        TSFigs(countFigs).name = 'KD Results By Cost';
        
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
                        D=createDistanceMatrix(costTotal',mapper);
                        [sD,idxCMap]=sort(abs(D),2);
%  %       
%  % for each Edge find the xy coords of the seed point and the candidate
%  % point and plot 
%  for i = 1:length(E(:,1))
%      % get the indexes of the edge relative to the original input EPs
%     idxCand = E(i,2);
%     idxSeed = E(i,1);
%     plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)]);
%     
% end
for k=1:cMapLength
    idxCand = E(idxCMap(:,1) == k,2);
    idxSeed = E(idxCMap(:,1)==k,1);
    for iEdge = 1:length(idxCand) % some can have the same color
        plot([internalCandEPs(idxCand(iEdge),1),internalSeedEPs(idxSeed(iEdge),1)],[internalCandEPs(idxCand(iEdge),2),internalSeedEPs(idxSeed(iEdge),2)],'color',cMap(k,:),'Linewidth',1);
    end
end
text(5,5,'Color Paths By Cost ','FontSize',10,'Color','k'); 
        text(5,15,'Red : High : Stong Path' ,'FontSize',10,'Color','r'); 
        text(5,25,'Blue : Low : Weak Path', 'FontSize',10,'Color','b'); 
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
    %idxGood = find(abs(costIntAndConn)>0.70 & abs(costSeedAndConn)>0.70 & abs(costCandAndSeed)>0.70); % quickest fix is to up the cost
    idxGood = find(abs(costIntAndConn)>0.7 & abs(costSeedAndConn)>0.7);  
    costTotal = costTotal(idxGood)';
    E = E(idxGood',:);   
    idxCMap = idxCMap(idxGood',:); 
%% TSOverlays: Filter  
%% TSOverlays : Plot Linear Connections Color-Coded by Cost 
if ip.Results.TSOverlays == true;
  
   
        
        TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'on'); % reget the handle 
        TSFigs(countFigs).name = 'KD Results AFter Filter By Geometry';
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
%      cMapLength=128; cMap=jet(cMapLength);
%                         mapper=linspace(min(costTotal),max(costTotal),cMapLength)';
%                        
%                         % get closest colormap index for each feature
%                         D=createDistanceMatrix(costTotal',mapper);
%                         [sD,idxCMap]=sort(abs(D),2);
%  %       
%  % for each Edge find the xy coords of the seed point and the candidate
%  % point and plot 
%  for i = 1:length(E(:,1))
%      % get the indexes of the edge relative to the original input EPs
%     idxCand = E(i,2);
%     idxSeed = E(i,1);
%     plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)]);
%     
% end
for k=1:cMapLength
    idxCand = E(idxCMap(:,1) == k,2);
    idxSeed = E(idxCMap(:,1)==k,1);
    for iEdge = 1:length(idxCand) % some can have the same color
        plot([internalCandEPs(idxCand(iEdge),1),internalSeedEPs(idxSeed(iEdge),1)],[internalCandEPs(idxCand(iEdge),2),internalSeedEPs(idxSeed(iEdge),2)],'color',cMap(k,:),'Linewidth',1);
    end
end
  countFigs = countFigs+1; 
  % c = colormap(lines(size(E,1)));
    % for i = 1:length(E(:,1))
    %     idxCand = E(i,2);
    %     idxSeed = E(i,1);
    %     plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)],'color',c(i,:));
    %     texend % end ip.Results.TSOverlays 

    %
    % end

end % ip.Results.

    %% Use Graph Matching to choose quickly the best candidate among several competing edges
  
    [seedFiloNodes,~,nodeLabels] = unique(E(:,1),'stable'); % reason note: some of the filo will not be candidates as their endpoints are not within the given radius
    NNodeQuery = length(seedFiloNodes);
    [inputLinks,~,nodeLabelsInput] = unique(E(:,2),'stable'); % reason note: just in case two filo are competing over the same seed point
    nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
    EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
    numberNodes = length(inputLinks) + length(seedFiloNodes);
    
    
    
    M = maxWeightedMatching(numberNodes, EFinal, costTotal);
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
    TSFigs(countFigs).h  = setFigure(imSize(2),imSize(1),'on'); % reget the handle 
      TSFigs(countFigs).name = 'Linking Results'; 
    
      
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
else
    status = 0 ; % not links within range to even consider
    maskPostConnect = seedMask;
    linkMask = zeros(imSize);
end






