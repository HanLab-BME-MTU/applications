function [ candidateMaskNew,linkMask,EPCandidateSort, labelMatFill,status] = connectLinearStructures(EPCandidateSort,maxTh,candidateMask,labelMat,linearityThresh,radius)
%connectLinearStructure: small function to connect the ends of ridges that
%fall within a certain distance/certain linearity criteria. Used for
%connecting backbone (large scale ridges of neurite) or as a preclustering
%step to connect smaller segments of filopodia for further reconstruction 


% INPUT:
% PERSONAL NOTE: maybe change the input such that it just takes in the
% binary mask of filaments? not a priority though 

% EPCandidateSort: a cell array 
% linearityThresh : a vector of values for the colinearity threshold. 1)
% between candidate1 and path , 2) between candidate2 and path, and 3) between the two candidates. Default (0.5,0.5,0.5)   
% 

maxTh = maxTh + pi/2; 
endPoints = vertcat(EPCandidateSort{:}); % taking these out of a cell array so 
% create a repmat for indexing 
labels = arrayfun(@(i) repmat(i,2,1),1:numel(EPCandidateSort),'uniformoutput',0); 
labels = vertcat(labels{:}); 
[idx,d] = KDTreeBallQuery(endPoints, endPoints, radius); % originally 5 
% remove all with distance ==0 
 E = arrayfun(@(i) [repmat(i, [numel(idx{i}) 1]) idx{i}], 1:length(endPoints), 'UniformOutput', false);
E = vertcat(E{:}); 
d = vertcat(d{:});





 
E(E(:,1)<=E(:,2),:)=[]; % remove redundancy 
% remove self associations (ie connection with same label)
label1 = arrayfun(@(i) labelMat(sub2ind(size(maxTh),endPoints(E(i,1),2), endPoints(E(i,1),1))),1:length(E(:,1)))  ; 
label2 = arrayfun(@(i) labelMat(sub2ind(size(maxTh),endPoints(E(i,2),2),endPoints(E(i,2),1))),1:length(E(:,1))); 
E = E(label1~=label2,:); 
d = d(label1~=label2); 

sanityCheck = 1; 
if sanityCheck == 1
    figure;
    % plot edge labels
    imagesc(labelMat); 
    %hold on 
   
    
end 

d = max(d)-d;
% normalize 
d= d./max(d);



if ~isempty(E); % check if there are reasonable edges 
    

        
endpointIdx = sub2ind(size(maxTh),endPoints(:,2),endPoints(:,1));
% get 
t1 = maxTh(endpointIdx(E(:,1))); 
t2 = maxTh(endpointIdx(E(:,2))); 
  a1 = abs(t1 - t2);
        a2 = abs(a1-pi);
        minAngle = min(a1,a2);
        cost = cos(minAngle); % 


% get theta for those matches 
%   t1 = theta(endpointIdx(unmatchedIdx(E(:,1))));
%         t2 = theta(endpointIdx(unmatchedIdx(E(:,2))));
%         a1 = abs(t1 - t2);
%         a2 = abs(a1-pi);
%         minAngle = min(a1,a2);
%         cost1 = cos(minAngle);
%         
        % need to likewise need to maintain linearity 
        % among the connection and the two pieces
        % get orientation of the connecting piece
       
        deltX = arrayfun(@(i) endPoints(E(i,1),1)-endPoints(E(i,2),1),1:length(E(:,1))); 
        deltY= arrayfun(@(i) endPoints(E(i,1),2)-endPoints(E(i,2),2),1:length(E(:,1))); 
        angleConn= atan2(deltY,deltX)'; 
       
        deltAngleWithCan1 = abs(angleConn-t1); 
        
        
        
       test1 = abs(deltAngleWithCan1-pi);
      
        minAngle1 = min(deltAngleWithCan1,test1);
        deltAngleWithCand2 = abs(angleConn-t2);
        test2 = abs(deltAngleWithCand2-pi); 
       
        minAngle2 = min(deltAngleWithCand2,test2); 
        costCand1Path = cos(minAngle1); 
        costCand2Path = cos(minAngle2); 
        % FOR NOW JUST FILTER OUT BASED ON A HARD THRESHOLD ON ORIENTATION MATCHING 
        idxGood = find(costCand1Path>linearityThresh(1)& costCand2Path>linearityThresh(2) & cost>linearityThresh(3)); % 0.95,0.95,0.95 was .80/.80/.80
        
        
        costTotal = costCand1Path+costCand2Path+cost+d; % but use d as well when considering competing matches
        costTotal = costTotal(idxGood); 
        
       
        E = E(idxGood,:); 
        
        % just in case there are competing nodes 
%         [seedFiloNodes,~,nodeLabels] = unique(E(:,1),'stable'); % reason note: some of the filo will not be candidates as their endpoints are not within the given radius
% %         NNodeQuery = length(seedFiloNodes);
%         
%         [inputLinks,~,nodeLabelsInput] = unique(E(:,2),'stable'); % reason note: just in case two filo are competing over the same seed point
%         nodeLabelsInputFinal = nodeLabelsInput+NNodeQuery;
%         EFinal = [nodeLabels nodeLabelsInputFinal]; % put in independent node form
%         numberNodes = length(inputLinks) + length(seedFiloNodes);
%         
        [numberNodes,~,nodeLabels] =unique(E(:),'stable');
        numberNodes = length(numberNodes); 
        EFinal(:,1) = nodeLabels(1:length(E(:,1))); 
        EFinal(:,2) = nodeLabels(length(E(:,1))+1:end); 
        
        M = maxWeightedMatching(numberNodes, EFinal, costTotal);
        E = E(M,:);  
        
        
        
        paths=arrayfun(@(i) bresenham([endPoints(E(i,1),1) endPoints(E(i,1),2)], [endPoints(E(i,2),1) endPoints(E(i,2),2)]),...
            1:length(E(:,1)),'uniformoutput',0); 
        linkMask = zeros(size(maxTh)); 
       
        links = vertcat(paths{:});
        
else % no goo links 
    
         links = []; % since links is empty this will default to mask bad
          linkMask = zeros(size(maxTh)); 
end % isempty E 
        
        
        
         if ~isempty(links) % nothing that falls under this criteria
        % Add links to candidate mask 
        idxLinks = sub2ind(size(maxTh),links(:,2),links(:,1)); 
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
          labelMatFill = labelMat; % initiate a labelMatrix that connects the two pieces
       for i =  1:length(E(:,1))
        label1 = labelMatFill(sub2ind(size(maxTh),endPoints(E(i,1),2), endPoints(E(i,1),1)));  
        label2 = labelMatFill(sub2ind(size(maxTh),endPoints(E(i,2),2),endPoints(E(i,2),1))); 
        test = [label1 ;label2]; 
       
        labelKeep  = max(test); 
        labelSwitch = min(test); 
        labelMatFill(labelMatFill==labelSwitch) = labelKeep; % change the label. 
  
        labelMatFill(sub2ind(size(maxTh),paths{i}(:,2),paths{i}(:,1)))= labelKeep; 
        
        % make new mask of the new connected filo and get endpoints 
       
        testMask = double(labelMatFill==labelKeep);
        sumKernel = [1 1 1];
    % find endpoints of the floating candidates to attach (note in the 
    % future might want to prune so that the closest end to the 
    % body is the only one to be re-attatched: this will avoid double connections)
    newEPs = double((testMask.* (conv2(sumKernel, sumKernel', padarrayXT(testMask, [1 1]), 'valid')-1))==1); 
    [ye,xe] = find(newEPs~=0); 
%     coords(:,1) = xe; 
%     coords(:,2) = ye;
       EPCandidateSort{labelKeep} = [xe ye]; 
       EPCandidateSort{labelSwitch} = [NaN NaN ;NaN NaN]; % take out those endpoints that were previously considered 
        
       end 
     status =1; % there were links 
         else 
             candidateMaskNew = candidateMask; 
             labelMatFill = labelMat; 
             status = 0; % no links 
        


end

