function [refinedForceMesh]=refineMeshAndBasis(forceMesh,nodesToRef,doPlot)
if nargin <3 || isempty(doPlot)
    doPlot        =0;
end
% find the Basis functions with these nodes.
for i=1:size(nodesToRef,1)
    found=false;
    j=1;
    while ~found || j<=size(forceMesh.basis,1)
        if compPts(nodesToRef(i,:),forceMesh.basis(j).node)
            refNode(i).node   =nodesToRef(i,:);
            refNode(i).basisID=j;
            refNode(i).nodeID =forceMesh.basis(j).nodeID;
            % find the neighboring nodes:
            class       = forceMesh.basis(j).class;
            neighRelPos = forceMesh.basisClass(class).neighPos;
            numNeigh    = forceMesh.basis(j).numNeigh;
            neigh       = forceMesh.basisClass(class).neighPos + repmat(forceMesh.basis(j).node,numNeigh,1);
            % assign the found neighbors:
            refNode(i).neigh= neigh;
            
            % so far we are done with this node:
            found=true;
        end
        j=j+1;
    end
    if found==false
        display('For node [',num2str(nodesToRef(i,1)),' ',num2str(nodesToRef(i,2)),'] no basis function has been found!');
        display('Nothing has been done with this node');
    end
end

newPts=[];
for i=1:length(refNode)
    % triangulize with the neighboring nodes:
    allPts=vertcat(refNode(i).node,refNode(i).neigh);
    dt=DelaunayTri(allPts(:,1),allPts(:,2));
    p=dt.X;
    t=dt.Triangulation;
    
    % all triangels per definition contain the center node.
    % determine the center of all connected triangles.
    for j=1:length(t)
        % get the point IDs for each triangle:
        ptsID=t(j,:);
        % read out the position of these pts:
        ptsPos=p(ptsID,:);
        % take the mean of the x and y pos to get the center of each
        % triangle. These are the new points to be added to the mesh
        center(j,:)=mean(ptsPos);        
    end
    newPts=vertcat(newPts,center);
end
% sort out double entries:
newPts=sortrows(newPts);
newPts=removeDoublePoints(newPts);

% append these points at the end of the original mesh points.
p_new    =vertcat(forceMesh.p,newPts);
%ID_newPts=(size(forceMesh.p,1)+1):size(p_new,1);

% use the extend point list to create a new force mesh with the same
% criterias as the original mesh:
keepBDPts     =forceMesh.keepBDPts;
basisClassInit=forceMesh.basisClass;
[refinedForceMesh]=createMeshAndBasisFastBEM(p_new(:,1),p_new(:,2),keepBDPts,basisClassInit,doPlot);

% now compare new and old force mesh to see which basis function have to be
% recalculated!
updateVec=ones(refinedForceMesh.numBasis,1);
for basisID=1:forceMesh.numBasis
    updateVec(basisID)=~compareBasisFunctions(forceMesh,refinedForceMesh,basisID);
end
refinedForceMesh.updateVec=updateVec;
refinedForceMesh.IDsOfBasisSolToBeUpdated=find(updateVec);
