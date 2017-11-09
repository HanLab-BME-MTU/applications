function [aggregationProbVec,receptor2cluster,cluster2receptor,clusterSize,...
    positionsNew,numPotColl,numColl,numPotColl_Assoc,numCollProbPairs,collProbStats] =...
    receptorAggregationAlg_maxWeightedMatching_potColl(...
    cluster2receptor,receptor2cluster,positionsOld,positionsNew,...
    aggregationProbVec,aggregationProb,associationDistVals,clusterSize)
%RECEPTAGGREGATIONALG_MAXWEIGHTEDMATCHING_POTCOLL performs receptor-to-receptor
%   and receptor-to-cluster associations for those cases where sure
%   associations did not occur. maxWeightedMatching is used to perform 
%   associations, limiting all receptors/clusters to associate with only 
%   one other receptor. Cluster-to-cluster associations are not allowed.
%
%   This is the function that implements the "path-based" association and
%   associations based on probability of collision. The latter is applied
%   only to those cases that do not collide.  Refer to summaries from
%   08/20/14, 11/11/14 and 12/09/14 for details.
%   
%   INPUT:
%       cluster2receptor:       2D array of clusters along the rows and
%                               member receptors along the columns
%       receptor2cluster:       1D array of current cluster labels for each
%                               receptor (rows = receptors)
%       positionsOld:           2D array of receptor positions at current
%                               iteration
%       positionsNew:           2D array of receptor positions for next
%                               iteration
%       aggregationProbVec:     1D array of receptor aggregation flags
%       aggregationProb:        1D array giving domain from which 
%                               aggregation probabilities will be drawn,
%                               for each cluster size (formation). Element
%                               #1 is assumed to be a default value, used
%                               when aggregationProb is a scalar or when
%                               its size < a given cluster size.
%       associationDistVals:    a struct with the following two fields:
%                               aggregationDist - the distance threshold 
%                               for association
%                               aggregationDistCorr - the corrected 
%                               distance threshold
%       clusterSize:            1D array of cluster sizes 
%
%
%   OUTPUT:
%       aggregationProbVec:     updated 1D array of receptor aggregation 
%                               flags
%       receptor2cluster:       1D array of cluster labels for each receptor
%       cluster2receptor:       2D array of clusters along the rows and
%                               member receptors along the columns
%       clusterSize:            1D array of cluster sizes
%       positionsNew:           updated 2D array of recpetor positions
%       numPotColl:             number of potential collision that occured
%       numColl:                number of collision that occured out of the
%                               potential collisions
%       numPotColl_Assoc:       number of actual associations that occured
%       numCollProbPairs:       number of pairs colliding via collision
%                               probability
%       collProbStats:          a struct with the following fields
%                               detailing collisions via collision prob.:
%                               collisionProb - the probablity value
%                               pwDist - pairwise distance
%                               primaryNodeRadius - radius for primary node
%                               partnerNodeRadii - set of radii for partner
%                                                  nodes
%
%   Robel Yirdaw, July 2014
%       Modified, November 2014
%

    %Initialize varabiles that will track potentail and probablity based
    %collisions. This was added inorder to probe the implementation.
    numPotColl = 0;
    numPotColl_Assoc = 0;
    numColl = 0;
    numCollProbPairs = 0;
    
    %Statistics for collision probabilities will also be returned.
    collProbStats = struct('collisionProb',NaN,'pwDist',NaN,...
    'primaryNodeRadius',NaN,'partnerNodeRadii',NaN);

    %Pull out default and corrected association distances.
    aggregationDist = associationDistVals.associationDist;
    aggregationDistCorr = associationDistVals.associationDistCorr;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Construct nodes, get pairwise distances and save those with potential
    %to associate.    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Current nodes list
    currNodeList = makeNodesFromC2R(cluster2receptor,positionsOld,aggregationProbVec);
    
    %Get ID of those able to associate
    assocNodeID_assocFlag = find([currNodeList.assocFlag] == 1);
    
    %Pull out coordinates for the above
    tempNodePositions = [([currNodeList.xPos]'), ([currNodeList.yPos]')];
    
    %Get pairwise distances
    assocNodePWDist = squareform(pdist(tempNodePositions(assocNodeID_assocFlag,:)));
    
    %Pick those that are within 2*rmsd but at or outside of aggregationDist
    assocNodeIndx_PWDist = find(any((assocNodePWDist >= aggregationDist) &...
        (assocNodePWDist < aggregationDistCorr) &...
        (assocNodePWDist ~= 0) ) );
    
    
    if (~isempty(assocNodeIndx_PWDist))
        %ID of nodes corresponding to assocNodesIndx_PWDist to be checked
        assocNodeID_check = assocNodeID_assocFlag(assocNodeIndx_PWDist);
        
        numNodePairs = 0.5*(sum(sum((assocNodePWDist >= aggregationDist) &...
            (assocNodePWDist < aggregationDistCorr) &...
            (assocNodePWDist ~= 0) ) ) );
        
        %Initialize list of potential node pairs and pair counter
        nodePairsToMerge = NaN(numNodePairs,2);
        mergeWeights = NaN(numNodePairs,1);
        pairCount = 0;
        
        for nodeIter=1:length(assocNodeID_check)
            %The current node ID
            primaryNodeID = assocNodeID_check(nodeIter);
            %Get potential association partners for primary node
            partnerNodeIndx = find((assocNodePWDist(:,assocNodeIndx_PWDist(nodeIter)) >= aggregationDist) &...
                (assocNodePWDist(:,assocNodeIndx_PWDist(nodeIter)) < aggregationDistCorr) &...
                (assocNodePWDist(:,assocNodeIndx_PWDist(nodeIter)) ~= 0) );
            %Trim upper values - using the lower triangle entries in PWDist
            partnerNodeIndx(partnerNodeIndx < assocNodeIndx_PWDist(nodeIter)) = [];
            %Finally get ID of partners
            partnerNodeID = assocNodeID_assocFlag(partnerNodeIndx);
            
            %Cluster to cluster associations are not allowed
            if ( (currNodeList(primaryNodeID).size > 1) &&...
                    any([currNodeList(partnerNodeID).size] > 1) )
                %Get boolean of parnter nodes with size > 1
                cluster2clusterAssoc = ([currNodeList(partnerNodeID).size] > 1);
                %Remove these from partner node list
                partnerNodeID(cluster2clusterAssoc) = [];
                
                clear cluster2clusterAssoc
            end
            
            if (~isempty(partnerNodeID))
                numPartners = length(partnerNodeID);
                
                %070914: sum the number of potential collisions (global)
                numPotColl = numPotColl + numPartners;
                
                %For each partner, determine separation distances during
                %displacement from old to new position
                
                %Initialize array of separation distance and time values 
                closestDist = NaN(numPartners,1);                
                timeOfCollision = NaN(numPartners,1);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %{
                %To use distances per incrmental steps, need to pass into
                %this function numDispPerUnitTime and dispTimeVec                
                %timeOfClosestDist = NaN(numPartners,1);
                %}
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for partnerIndx=1:numPartners
                    %Construct old and new coordinates for the current
                    %primary node and its partner to pass to function
                    
                    tempNode1.xi = currNodeList(primaryNodeID).xPos;
                    tempNode1.yi = currNodeList(primaryNodeID).yPos;
                    %New positions are in separate vector - positionsNew
                    tempNode1.xf = positionsNew(currNodeList(primaryNodeID).members(1,1),1);
                    tempNode1.yf = positionsNew(currNodeList(primaryNodeID).members(1,1),2);
                    
                    tempNode2.xi = currNodeList(partnerNodeID(partnerIndx)).xPos;
                    tempNode2.yi = currNodeList(partnerNodeID(partnerIndx)).yPos;
                    %New positions are in separate vector - positionsNew
                    tempNode2.xf = positionsNew(currNodeList(partnerNodeID(partnerIndx)).members(1,1),1);
                    tempNode2.yf = positionsNew(currNodeList(partnerNodeID(partnerIndx)).members(1,1),2);
                    
                    %{  %} 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Determine if a collision has occured between the two
                    %nodes (path-based association).
                    [timeOfCollision(partnerIndx,1),closestDist(partnerIndx,1)] =...
                        calcTimeOfCollision(tempNode1,tempNode2);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    clear tempNode1 tempNode2
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    %{ 
                    %To use distances per incrmental steps instead, need to
                    %pass into this function numDispPerUnitTime and dispTimeVec. 
                    %calcTimeOfCollision above has to be blocked.
                    
                    %Collect distances for each pairing
                    distValues = calcPairDistPerUnitTime(...
                        tempNode1,tempNode2,numDispPerUnitTime);
                        
                    %Save the closest distance per partner and its
                    %corresponding time point. NOTE: more than one min
                    %value possible, so take the first one in that case
                    
                    %closestDistVals(partnerIndx) = min(distValues);
                    %minDist_indx = (distValues == min(distValues));   
                    minDist_indx = find(distValues == min(distValues),1,'first');
                    
                    %Save the closest distance of approach
                    closestDist(partnerIndx,1) = distValues(minDist_indx);
                    %... and the corresponding time
                    timeOfClosestDist(partnerIndx,1) = dispTimeVec(minDist_indx);                                                          
                    
                    clear tempNode1 tempNode2 distValues minDist_indx                    
                    %}
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                end %for all partners of the current primary node
                    
                %Determine if at any point during displacement any two
                %pairs had distances that satisfied the required condition
                %for association. If so, save the pair.                                
                closestPartnerNode_indx = find(closestDist(:,1) < aggregationDist);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %081814
                %Adding collision probability for those not handled by the
                %above case - no detected collision. 
                collProbPartnerNode_indx = [];
                if (isempty(closestPartnerNode_indx))
                    collProbPartnerNode_indx(:,1) = 1:numPartners;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %070914: sum the number of collisions that occured (global)                
                numColl = numColl + length(closestPartnerNode_indx);                
                
                if (~isempty(closestPartnerNode_indx))
                    %Possibile that more than two nodes have 
                    %intersectiong paths or paths within association
                    %distance. In this case closestPartnerNode_indx will
                    %have more than one elment.
                    
                    numClosestPartners = length(closestPartnerNode_indx);
                    %Pick node IDs for the partners
                    closestPartnerNodeID = partnerNodeID(closestPartnerNode_indx);               
                           
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Prep aggregationProb
                    aggregationProbTemp = aggregationProb;
                    finalSizes = currNodeList(primaryNodeID).size + [currNodeList(closestPartnerNodeID).size];
                    
                    %If all sizes are not accounted for in aggregationProb,
                    %expand it and assign default value which is element 1.
                    if (length(aggregationProbTemp) < max(finalSizes))
                        aggregationProbTemp(length(aggregationProbTemp)+1:max(finalSizes),1) =...
                            aggregationProbTemp(1,1);
                    end
                    %Draw a random value for each cluster to be formed,
                    %using its size, and get association flag.
                    eventFlag = (rand(numClosestPartners,1) < aggregationProbTemp(finalSizes) );
                    
                    %Remove from the list partner nodes not associating
                    closestPartnerNodeID(eventFlag == 0) = [];
                    closestPartnerNode_indx(eventFlag == 0) = [];       
                    
                    clear eventFlag finalSizes aggregationProbTemp
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %081814
                %Continue with those nodes considered for collision
                %probability just like the collision case
                %
                if (~isempty(collProbPartnerNode_indx))

                    %Get number in this group
                    numCollProbPartners = length(collProbPartnerNode_indx);
                    %Pick node IDs for the partners
                    collProbPartnerNodeID = partnerNodeID(collProbPartnerNode_indx);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Prep collisionProb                    
                    
                    %Determine center coordinates for primary and partner
                    %nodes. The centers is the midpoint of the path taken.                    
                    primaryNode_c(1:numCollProbPartners,1) = currNodeList(primaryNodeID).xPos +...
                        0.5*( positionsNew(currNodeList(primaryNodeID).members(1,1),1) -...
                        currNodeList(primaryNodeID).xPos );
                    primaryNode_c(1:numCollProbPartners,2) = currNodeList(primaryNodeID).yPos +...
                        0.5*( positionsNew(currNodeList(primaryNodeID).members(1,1),2) -...
                        currNodeList(primaryNodeID).yPos);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Next we need the two radii, which are equal to half
                    %the step taken by each particle.                    
                    primaryNodeRadius = 0.5*sqrt( (positionsNew(currNodeList(primaryNodeID).members(1,1),1) -...
                        currNodeList(primaryNodeID).xPos)^2 + (positionsNew(currNodeList(primaryNodeID).members(1,1),2) -...
                        currNodeList(primaryNodeID).yPos)^2 );                    
                    primaryNodeRadiusVec = repmat(primaryNodeRadius,numCollProbPartners,1);
                    
                    %Initialize partner vectors
                    partnerNode_c = NaN(numCollProbPartners,2);
                    partnerNodeRadii = NaN(numCollProbPartners,1);
                    
                    for tempIndx=1:numCollProbPartners
                        %Get x and y center coordinates
                        partnerNode_c(tempIndx,1) = currNodeList(collProbPartnerNodeID(tempIndx)).xPos +...
                            0.5*( positionsNew(currNodeList(collProbPartnerNodeID(tempIndx)).members(1,1),1) -...
                            currNodeList(collProbPartnerNodeID(tempIndx)).xPos );

                        partnerNode_c(tempIndx,2) = currNodeList(collProbPartnerNodeID(tempIndx)).yPos +...
                            0.5*( positionsNew(currNodeList(collProbPartnerNodeID(tempIndx)).members(1,1),2) -...
                            currNodeList(collProbPartnerNodeID(tempIndx)).yPos ); 

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %Next we need the two radii, which are equal to half
                        %the step taken by each particle.                    
                        partnerNodeRadii(tempIndx,1) = 0.5*sqrt( (positionsNew(currNodeList(collProbPartnerNodeID(tempIndx)).members(1,1),1) -...
                            currNodeList(collProbPartnerNodeID(tempIndx)).xPos)^2 + (positionsNew(currNodeList(collProbPartnerNodeID(tempIndx)).members(1,1),2) -...
                            currNodeList(collProbPartnerNodeID(tempIndx)).yPos)^2 );
                    
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    %Calculate the pairwise distance between the center
                    %points of primary node's path and that of it's
                    %partners. This will be used to calculate the collision
                    %probability below.                    
                    pwDist = sqrt( (partnerNode_c(:,1) - primaryNode_c(:,1)).^2 +...
                        (partnerNode_c(:,2) - primaryNode_c(:,2)).^2 );
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Get colllision probability
                    collisionProb = arrayfun(@calcCollisionProb_circlesIntArea,...
                        pwDist,primaryNodeRadiusVec,partnerNodeRadii);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Prep aggregationProb
                    aggregationProbTemp = aggregationProb;
                    finalSizes = currNodeList(primaryNodeID).size + [currNodeList(collProbPartnerNodeID).size];
                    
                    %If all sizes are not accounted for in aggregationProb,
                    %expand it and assign default value which is element 1.
                    if (length(aggregationProbTemp) < max(finalSizes))
                        aggregationProbTemp(length(aggregationProbTemp)+1:max(finalSizes),1) =...
                            aggregationProbTemp(1,1);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    %Draw a random value for each cluster to be formed,
                    %and get association flag using the product of the
                    %association and collision probabilities
                    eventFlag = (rand(numCollProbPartners,1) <...
                        aggregationProbTemp(finalSizes).*collisionProb );

                    
                    %Remove from the list partner nodes not associating
                    collProbPartnerNodeID(eventFlag == 0) = [];
                    collProbPartnerNode_indx(eventFlag == 0) = [];  
                    %collisionProb will be used as weight
                    collisionProb(eventFlag == 0) = [];
                    
                    pwDist(eventFlag == 0) = [];
                    partnerNodeRadii(eventFlag == 0) = [];

                    clear eventFlag finalSizes aggregationProbTemp
                    %clear primaryNodeRadius primaryNodeRadiusVec partnerNodeRadii pwDist
                    clear primaryNode_c partnerNode_c tempIndx
                    
                end % collProbPartnerNode_indx exists 
                
                
                numClosestPartners = length(closestPartnerNode_indx);
                if (numClosestPartners > 0)
                    %Save the primary node and its partners for association
                    nodePairsToMerge(pairCount+1:(pairCount+numClosestPartners),1) = primaryNodeID;
                    nodePairsToMerge(pairCount+1:(pairCount+numClosestPartners),2) = closestPartnerNodeID; 
                    %Inverse of time will be used as weight
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    %{ 
                    %To use distances per incrmental steps, need to pass into
                    %this function numDispPerUnitTime and dispTimeVec                                    
                    %mergeWeights(pairCount+1:(pairCount+numClosestPartners),1) = 1./timeOfClosestDist(closestPartnerNode_indx); 
                    %}
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                    %mergeWeights(pairCount+1:(pairCount+numClosestPartners),1) = 1./timeOfCollision(closestPartnerNode_indx); 
                    
                    %082514 - use 1-time instead of 1/time
                    mergeWeights(pairCount+1:(pairCount+numClosestPartners),1) = 1 - timeOfCollision(closestPartnerNode_indx); 

                    %Increment pair counter
                    pairCount = pairCount + numClosestPartners;
                    
                    clear closestPartnerNodeID closestPartnerNode_indx numClosestPartners
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %081814
                %Continue with those nodes considered for collision
                %probability just like the collision case
                %                
                numCollProbPartners = length(collProbPartnerNode_indx);
                numCollProbPairs = numCollProbPairs + numCollProbPartners;
                if (numCollProbPartners > 0)
                    %Save the primary node and its partners for association
                    nodePairsToMerge(pairCount+1:(pairCount+numCollProbPartners),1) = primaryNodeID;
                    nodePairsToMerge(pairCount+1:(pairCount+numCollProbPartners),2) = collProbPartnerNodeID; 
                    
                    %Inverse of the collision probability will be used as weight
                    mergeWeights(pairCount+1:(pairCount+numCollProbPartners),1) = collisionProb; 
                    
                    %Increment pair counter
                    pairCount = pairCount + numCollProbPartners;
                    
                    %collProbStats = [];
                    collProbStats.collisionProb = collisionProb;
                    collProbStats.pwDist = pwDist;
                    collProbStats.primaryNodeRadius = primaryNodeRadius;
                    collProbStats.partnerNodeRadii = partnerNodeRadii;
                    
                    clear collProbPartnerNodeID collProbPartnerNode_indx numCollProbPartners
                    clear primaryNodeRadius primaryNodeRadiusVec partnerNodeRadii pwDist
                end
                
                
                clear numPartners
                
            end %if one or more partner nodes exist (~empty(partnerNodeID))
            
            clear primaryNodeID partnerNodeID partnerNodeIndx
            
        end %for each node with potential associations (primaryNodeIDs) 
        
        %If found node pairs that should associate, construct new nodes,
        %assign approprate positions and insert into node list
        nodePairsToMerge_Selected = [];
        if (pairCount == 1)
            %First trim empty rows in node pairs list
            nodePairsToMerge(pairCount+1:end,:) = [];
            nodePairsToMerge_Selected = nodePairsToMerge;
        elseif (pairCount > 1)
            %First trim empty rows in node pairs list
            nodePairsToMerge(pairCount+1:end,:) = [];
            mergeWeights(pairCount+1:end,:) = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Peform weighted matching to pick a single pairwise
            %interaction for every node               
            selectionFlags = maxWeightedMatching(max(max(nodePairsToMerge(:,1:2))),...
                nodePairsToMerge(:,1:2),mergeWeights);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Per returned flags above, pick the edges
            nodePairsToMerge_Selected = nodePairsToMerge(selectionFlags,:);
            
        end
            
        %Process the selected pairs and finalize    
        if (~isempty(nodePairsToMerge_Selected))
            %Insert the node pairs save above, as merged nodes, into the
            %current node list. This will also update receptPositions which
            %is returned with the updated node list.            

            %Insert
            [updatedNodeList,positionsNew] = insertNodes(nodePairsToMerge_Selected,currNodeList,positionsNew);
            
            %Next update receptor2cluster,cluster2receptor
            [receptor2cluster,cluster2receptor,clusterSize] =...
                updateR2C2RwithNodeInfo(updatedNodeList);            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %After performing associations based on node paths, if any, need to
            %disable aggregationProbVec for all those receptors that participated.
            %Note that since nodePairsToMerge has two columns, using linear
            %indexing to access current node list.    
            aggregationProbVec([currNodeList(nodePairsToMerge_Selected(:)).members]) = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %070914: sum the number of potential associations that resulted
            %in an actuall association (global)
            numPotColl_Assoc = numPotColl_Assoc + length(nodePairsToMerge_Selected(:,1));

        end %if found node pairs to associate
            
            
    end %if pairs within assoc distance exist  
  

end %function

