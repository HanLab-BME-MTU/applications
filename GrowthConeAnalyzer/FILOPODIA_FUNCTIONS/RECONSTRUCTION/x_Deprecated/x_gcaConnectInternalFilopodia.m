function [maskPostConnect,linkMask,status] = gcaConnectInternalFilopodia(internalCandEPs,internalSeedEPs,maxTh,seedMask,searchRadius,labelMat,vectInt,vectSeed,dInt,dSeed)
%% SUMMARY: Connects Co-linear Candidate Internal Ridges (a filoTracker Filo Reconstruction Function)
% This function will take the internal filopodia endpoints and attempt to
% reconnect them to the user defined seed points (typically high-confidence external/internal filo)
% Candidates are filtered
% linking the two candidate ridges
% The cost between the candidate ridges and .
%% INPUT:

% internalCanEPs: an nx2 array of candidate endpoints: these will be input
%                 points (ie you are searching for points in this population
%                 that meet the specified dist criteria )
%
% internalSeedEPs: an nx2 array of candidate endpoints of high confidence
% ridge sigal (connected to cell edge-or added in a previous reconstruction step) that will serve as query pts
% for the KDTree Ball Query search

% maxTh: nxm matrix the size of the img
% gives the local orientation of the ridge at a given point so one need not recalculate
% used as matching criteria.
%
% searchRadius: scalar, only interalCanEPs that are within this max search radius
% around each internalSeedEP (inputPt) are considered for matching.
%%%% NOTE you are going to need a label matrix of the candidates to keep
%%%% ONLY those that you connected. otherwise the output is problematic..
%%%% the reson why many of these were broken in the first place was due to
%%%% junctions therefore putting them back together without the filtering just recreates the
%%%% mess. ( did you take care of that 2013_07-14
%
% OUTPUT:
% maskPostConnect: 2D double matrix the size of the image, binary mask of all ridges (ie filopdia) post connection
% linkMask: 2D double matrix the size of the image , binary mask of all viable links between ridge candidates made
% status: scalar, 1 if links were found, 0 if no viable links were found
%%
maxTh = maxTh + pi/2; % orientation is documented by ridge detector perpendicular to feature
% find all query points (internal filo EP coordinates) within a search radius surrounding the seed points (external filo EP coords)
[idx,d] = KDTreeBallQuery(internalCandEPs, internalSeedEPs, searchRadius);
% remove all with distance ==0
% idxInput{iQuery} = [idxInput1, idxInput2, idxIndput3...]);
% first col E = idx filo Cand, 2nd col = idx input point,

% tranfer data from cell to double:  column 1 = index of the seed while
% column 2 = the index of the candidate, rows =  n number of
% potential paths
E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:length(internalSeedEPs(:,1)), 'UniformOutput', false);
E = vertcat(E{:});
% sanity check
spy(labelMat,'b');
hold on
spy(seedMask,'r');
for i = 1:length(E(:,1))
    idxCand = E(i,2);
    idxSeed = E(i,1);
    plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)]);
    
end





%% sanity check
% plot all paths
examplePlot =1 ;
%  if examplePlot ==1 ;
% %      candFiloNum = 1;% the number of the candidate filo you would like to plot
% %      EPlot = E(E(:,1)==candFiloNum,:);
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
%%

% E(E(:,1)<=E(:,2),:)=[]; % remove redundancy and self-queries
% e is of the form col = 1  of inputPt seed pt and col
if ~isempty(E) % continue
    
    idxCandAll = E(:,2);
    idxSeedAll = E(:,1);
    
    vectInt = vertcat(vectInt{idxCandAll});
    dInt = vertcat(dInt{idxCandAll});
    vectSeed = vertcat(vectSeed{idxSeedAll});
    dSeed = vertcat(dSeed{idxSeedAll});
    
    
    
    costCandAndSeed = arrayfun(@(i) dot(vectInt(i,:),vectSeed(i,:))./dInt(i)./dSeed(i),1:length(dInt));
    
    EPidxSeed = sub2ind(size(maxTh),internalSeedEPs(:,2),internalSeedEPs(:,1));
    EPidxCand = sub2ind(size(maxTh),internalCandEPs(:,2),internalCandEPs(:,1));
    
    %% NOTE: 2013_07_22: Change from taking the orientation using the end points
    % to using the full disp vector
    
    % get
    % t1 = maxTh(EPidxSeed(E(:,1)));
    % t2 = maxTh(EPidxCand(E(:,2)));
    % % NOTE try taking the average. get pix indices of all labels
    % % test first to make sure linear otherwise get some artifacts at ends!!
    %
    %   a1 = abs(t1 - t2);
    %         a2 = abs(a1-pi);
    %         minAngle = min(a1,a2);
    %         cost = cos(minAngle);
    
    
    % get theta for those matches
    %   t1 = theta(endpointIdx(unmatchedIdx(E(:,1))));
    %         t2 = theta(endpointIdx(unmatchedIdx(E(:,2))));
    %         a1 = abs(t1 - t2);
    %         a2 = abs(a1-pi);
    %         minAngle = min(a1,a2);
    %         cost1 = cos(minAngle);
    %
    
    % calculate overall displacement vector
    %      arrayfun(@(x) plot(internalSeedEPs(E(i,:),1)),internalSeedEPs(E(i,:),2),1:length(E(:,1)));
    %      arrayfun(@(x) text(internalSeedEPs(E(i,1),1)),internalSeedEPs(
    
    % need to likewise need to maintain linearity
    % among the connection and the two pieces
    % get orientation of the connecting piece
    deltX = arrayfun(@(i) internalSeedEPs(E(i,1),1)-internalCandEPs(E(i,2),1),1:length(E(:,1)));
    deltY= arrayfun(@(i) internalSeedEPs(E(i,1),2)-internalCandEPs(E(i,2),2),1:length(E(:,1)));
    vectConn = [deltX' deltY'];
    dConn = sqrt((deltX').^2+ (deltY').^2);
    costIntAndConn = arrayfun(@(i) dot(vectConn(i,:),vectInt(i,:))./dConn(i)./dInt(i),1:length(dConn));
    costSeedAndConn = arrayfun(@(i) dot(vectConn(i,:),vectSeed(i,:))./dConn(i)./dSeed(i),1:length(dConn));
    
    
    
    %         angleConn= atan2(deltY,deltX)';
    
    %         deltAngleWithCan1 = abs(angleConn-t1);
    %         test1 = abs(deltAngleWithCan1-pi);
    %         minAngle1 = min(deltAngleWithCan1,test1);
    %         deltAngleWithCand2 = abs(angleConn-t2);
    %         test2 = abs(deltAngleWithCand2-pi);
    %         minAngle2 = min(deltAngleWithCand2,test2);
    %         costCand1Path = cos(minAngle1);
    %         costCand2Path = cos(minAngle2);
    %
    
    
    
    % for now just filter out;
    idxGood = find(abs(costIntAndConn)>0.70 & abs(costSeedAndConn)>0.70 & abs(costCandAndSeed)>0.70); % quickest fix is to up the cost
    % Sanity check
    %
    % could also use more of the local orientation instead of just the
    % ends.
    costTotal =  abs(costIntAndConn) + abs(costSeedAndConn) + abs(costCandAndSeed);
    costTotal = costTotal(idxGood)';
    %        % sanity check
    %        figure;
    %         spy(labelMat,'b',10);
    % hold on
    % spy(seedMask,'r',10);
    % c = colormap(lines(size(E,1)));
    % for i = 1:length(E(:,1))
    %     idxCand = E(i,2);
    %     idxSeed = E(i,1);
    %     plot([internalCandEPs(idxCand,1),internalSeedEPs(idxSeed,1)],[internalCandEPs(idxCand,2),internalSeedEPs(idxSeed,2)],'color',c(i,:));
    %     text(internalCandEPs(idxCand,1),internalCandEPs(idxCand,2),num2str(costCand1Path(idxCand),2),'color',c(i,:));
    %
    % end
    
    E = E(idxGood',:);
    
    %% just in case put through graph matching
    % this is likely a bit overkill and maybe not the fastest way
    % what we are doing here is simply picking the maximum weight of
    % if there are two competing nodes so that only 1 attachment can be
    % made (though sometimes that assumption is less useful if have
    % these triangle converging internal filo... hmmm. maybe only allow
    % for bifurcation at the very edge and make it have to be very
    % robust... (or we just ignore at this point and wait for the
    % little guys to poke their head out a bit more...
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
    linkMask = zeros(size(maxTh));
    goodCands = zeros(size(maxTh));
    % find labels of internal candidates to keep
    labels = arrayfun(@(i) labelMat(sub2ind(size(maxTh),internalCandEPs(E(i,2),2),internalCandEPs(E(i,2),1))), 1:length(E(:,1)));
    % find the indexing of those labels
    idxCandKeep = arrayfun(@(i) find(labelMat == i),labels,'uniformoutput',0);
    goodCands(vertcat(idxCandKeep{:})) = 1;
    
    links = vertcat(paths{:});
    
    if ~isempty(links) % nothing that falls under this criteria
        % Add links to candidate mask
        idxLinks = sub2ind(size(maxTh),links(:,2),links(:,1));
        
        linkMask(idxLinks) = 1;
        maskPostConnect = (linkMask|goodCands|seedMask);
        status = 1;
    else
        maskPostConnect = seedMask;
        status =0;
    end
else
    status = 0 ; % not links within range to even consider
    maskPostConnect = seedMask;
    linkMask = zeros(size(maxTh));
end






