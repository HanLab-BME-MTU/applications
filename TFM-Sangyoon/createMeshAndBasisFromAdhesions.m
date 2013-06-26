function [myMesh]=createMeshAndBasisFromAdhesions(x_vec,y_vec,paxImage,displField,pixelSize)
% createMeshAndBasisFromAdhesions creates mesh from adhesion information in
% paxImage. x_vec and y_vec serve as mask for mesh nodes.
% Sangyoon Han June 2013
disp('Creating force mesh using adhesion information in paxillin image ...')
paxImage = double(paxImage);
% Get Cell mask
pId = paxImage/max(paxImage(:));
% alpha = graythresh(pId);
%estimate the intensity level to use for thresholding the image
level1 = graythresh(pId); %Otsu
[~, level2] = cutFirstHistMode(pId,0); %Rosin
alpha = 0.1*level2 + 0.9*level1;

bwPI = im2bw(pId,alpha);
bwPI2 = bwmorph(bwPI,'clean');
bwPI3 = bwmorph(bwPI2,'erode',5);
bwPI3 = bwmorph(bwPI3,'dilate',4);
bwPI4 = bwmorph(bwPI3,'close',10);
bwPI4 = bwmorph(bwPI4,'dilate',1);
% bwPI5 = refineEdgeWithSteerableFilterGM(pId,bwPI4);
% In case that there is still islands, pick only the largest chunk
[labelPI,nChunk] = bwlabel(bwPI4);
if nChunk>1
    eachArea = zeros(nChunk,1);
    for k=1:nChunk
        currBWPI = labelPI==k;
        eachArea(k) = sum(currBWPI(:));
    end
    [~,indCellArea] = max(eachArea);
    bwPI4 = labelPI == indCellArea;
end
% Get FA mask
disp('detecting focal adhesions ...')
minSize = round((1000/pixelSize)*(200/pixelSize)); %adhesion limit=1um*.5um
maskFA = blobSegmentThreshold(paxImage,minSize,false,bwPI2);

% For each FA, get points in the area. 
[labelFA,nFA]=bwlabel(maskFA);
FANodeX=[];
FANodeY=[];
[x_mat, y_mat]=meshgrid(1:size(paxImage,2),1:size(paxImage,1));
skelMask = false(size(x_mat));
for k=1:nFA
    eachMaskFA = labelFA==k;
    % try skeleton
    eachMaskFAskel = bwmorph(eachMaskFA,'skel',Inf);
    idxEachFAskel = find(eachMaskFAskel);
    [eachFANodesY,eachFANodesX] = ind2sub(size(paxImage),idxEachFAskel(1:5:end));
    % from this, we can try to find if there is any points apart by 5 pixel
    % in the mask
    for p=1:length(eachFANodesX)
        skelMask = skelMask | sqrt((x_mat-eachFANodesY(p)).^2+(y_mat-eachFANodesX(p)).^2)<5;
    end
    marginMask = eachMaskFA & ~skelMask;
    margNodesY = [];
    margNodesX = [];
    if any(marginMask(:))
        [labelMargin,nMarg] = bwlabel(marginMask);
        for q=1:nMarg
            % get node at the margin skeleton
            eachMaskMarg = labelMargin == q;
            eachMargskel = bwmorph(eachMaskMarg,'skel',Inf);
            idxEachMargskel = find(eachMargskel);
            [eachMargNodesY,eachMargNodesX]=ind2sub(size(paxImage),idxEachMargskel(1:5:end));
            margNodesX = [margNodesX; eachMargNodesX];
            margNodesY = [margNodesY; eachMargNodesY];
        end
    end
    FANodeX = [FANodeX; eachFANodesX; margNodesX];
    FANodeY = [FANodeY; eachFANodesY; margNodesY];
end
% Get equi-spaced nodes on the cell boundary
B=bwboundaries(bwPI4,'noholes');
boundary = B{1};
% Get rid of points at the image boundary
indAtImBD = boundary(:,1)==1 | boundary(:,1)==size(paxImage,1) ...
            | boundary(:,2)==1 | boundary(:,2)==size(paxImage,2);
boundary(indAtImBD,:)=[];
bdNodesX = boundary(1:9:end,2);
bdNodesY = boundary(1:9:end,1);
% mask for band from edge
iMask = imcomplement(bwPI4);
distFromEdge = bwdist(iMask);
bandwidth = 8; % in um
bandwidth_pix = round(bandwidth*1000/pixelSize);
bandMask = distFromEdge <= bandwidth_pix;
bandMask = bandMask & bwPI4;
band_FA = bandMask & ~maskFA;
% find nascent adhesions from paxImage
disp('detecting nascent adhesions ...')
sigmaPSF_NA = 1.48;
[pstruct_NA,mask] = pointSourceDetection(paxImage, sigmaPSF_NA,...
    'alpha', 0.05, 'mask', band_FA);
NANodeX = ceil(pstruct_NA.x');
NANodeY = ceil(pstruct_NA.y');
NAs = [NANodeX NANodeY];
NAs = [NAs; bdNodesX bdNodesY];
% NAs + boundaryNodes are separated by 5 pixel 

idx = KDTreeBallQuery(NAs, NAs, 5);
valid = true(numel(idx),1);
for i = 1:numel(idx)
    if ~valid(i), continue; end
    neighbors_KD = idx{i}(idx{i}~=i);
    valid(neighbors_KD) = false;
end
NAs = NAs(valid, :);
NANodeX = NAs(:,1);
NANodeY = NAs(:,2);

% Get intNodes in interiorMask
interiorMask = bwPI4 & ~bandMask & ~maskFA;
[pstruct_intNA] = pointSourceDetection(paxImage, sigmaPSF_NA*1.5,...
    'alpha', 0.05, 'mask', interiorMask);
intNAs = [ceil(pstruct_intNA.x') ceil(pstruct_intNA.y')];
% Add hexagonal grid in the interior region
spacing = ceil(2000/pixelSize);
[intHexNAX,intHexNAY]=createHexGridInMask(spacing,interiorMask);
intNAs = [intNAs; intHexNAX,intHexNAY];
% intNAs are separated by 20 pixel 
idx = KDTreeBallQuery(intNAs, intNAs, 20);
valid = true(numel(idx),1);
for i = 1:numel(idx)
    if ~valid(i), continue; end
    neighbors_KD = idx{i}(idx{i}~=i);
    valid(neighbors_KD) = false;
end
intNAs = intNAs(valid, :);
intNANodeX = intNAs(:,1);
intNANodeY = intNAs(:,2);
% Get hexagonal nodes with 2000 nm spacing for bgdMask
bgdMask = ~bwPI4;
[bgdNodeX,bgdNodeY]=createHexGridInMask(spacing,bgdMask);
% For all nodes, get their undeformed postitions using displField
allNodesX=[FANodeX; NANodeX; intNANodeX; bgdNodeX];
allNodesY=[FANodeY; NANodeY; intNANodeY; bgdNodeY];
undNodesX=allNodesX;
undNodesY=allNodesY;

for ii=1:length(allNodesX)
    distToNode = sqrt(sum(([displField.pos(:,1) displField.pos(:,2)]- ...
                            ones(length(displField.pos(:,1)),1)*[allNodesX(ii) allNodesY(ii)]).^2,2));
                        
    [~,closest_ind] = min(distToNode);
    dispVec = [displField.vec(closest_ind,1) displField.vec(closest_ind,2)];
    undNodesX(ii) = allNodesX(ii) - dispVec(1);
    undNodesY(ii) = allNodesY(ii) - dispVec(2);
end

% finally, mask the undNodes with ROImask about which xvec and yvec has
% information
% get grid information for entire image
[entireNodeX,entireNodeY]=createHexGridInMask(spacing,true(size(bgdMask)));
% filter entireNodeX and Y first
xmin = min(x_vec);
xmax = max(x_vec);
ymin = min(y_vec);
ymax = max(y_vec);
nEN = length(entireNodeX);
idxEN = false(nEN,1);
for ii=1:nEN
    if entireNodeX(ii)>=xmin && entireNodeX(ii)<=xmax ...
            && entireNodeY(ii)>=ymin && entireNodeY(ii)<=ymax
        idxEN(ii) = true;
    end
end
filteredHNX=entireNodeX(idxEN);
filteredHNY=entireNodeY(idxEN);

nPoints = length(undNodesX);
idxROI = false(nPoints,1);
xmin = min(filteredHNX);
xmax = max(filteredHNX);
ymin = min(filteredHNY);
ymax = max(filteredHNY);
for ii=1:nPoints
    if undNodesX(ii)>=xmin && undNodesX(ii)<=xmax ...
            && undNodesY(ii)>=ymin && undNodesY(ii)<=ymax
        idxROI(ii) = true;
    end
end

xvec_d = allNodesX(idxROI);
yvec_d = allNodesY(idxROI);

% we use deformed (more straight at the boundaries) points for mesh
% generation to prevent many obtuse triangles.
%delaunay_mesh=delaunay(x_vec,y_vec);
disp('constructing mesh ...')
xyvec = [xvec_d yvec_d];
xyvec = unique(xyvec,'rows');
xvec_d = xyvec(:,1);
yvec_d = xyvec(:,2);
dt=DelaunayTri(xvec_d,yvec_d);
delaunay_mesh=dt.Triangulation;

% from dt.X, we get the undeformed position again
allNodesX = dt.X(:,1);
allNodesY = dt.X(:,2);
undNodesX = zeros(size(allNodesX));
undNodesY = zeros(size(allNodesY));
for ii=1:length(allNodesX)
    distToNode = sqrt(sum(([displField.pos(:,1) displField.pos(:,2)]- ...
                            ones(length(displField.pos(:,1)),1)*[allNodesX(ii) allNodesY(ii)]).^2,2));
                        
    [~,closest_ind] = min(distToNode);
    dispVec = [displField.vec(closest_ind,1) displField.vec(closest_ind,2)];
    undNodesX(ii) = allNodesX(ii) - dispVec(1);
    undNodesY(ii) = allNodesY(ii) - dispVec(2);
end

xvec=undNodesX;
yvec=undNodesY;
% k = convexHull(dt);
%for each node find all its neighbors:
for n=1:length(xvec)
    candidates=[];
    for k=1:length(delaunay_mesh(:,1))
        if delaunay_mesh(k,1)==n || delaunay_mesh(k,2)==n || delaunay_mesh(k,3)==n
            for m=1:length(delaunay_mesh(1,:))
                if length(candidates)==0
                    if delaunay_mesh(k,m)~=n
                        candidates=horzcat(candidates,delaunay_mesh(k,m));
                    end                    
                elseif delaunay_mesh(k,m)~=n && isfinite(sum(1./(candidates-delaunay_mesh(k,m))))
                    candidates=horzcat(candidates,delaunay_mesh(k,m));
                end                    
            end
        end
    end
    neighbors(n).cand=sort(candidates);
    
    % find the minimal rectangular region around the central node which includes
    % all neighboring nodes. This will be needed for integration:
    % store the bounds parameterized so that only triangulated area is
    % considered for integration
    % see if the node is on the convexhull
%     if ismember(n,k)
%         bounds(n).x=[min(xvec(candidates)) max(xvec(candidates))];
%         bounds(n).y=[min(yvec(candidates)) max(yvec(candidates))];   
%         % see if this min max rectangle exceeds the convex hull
%         if bounds(n).y(1)<yvec(n)
%             x1 = xvec(n);
%             y1 = yvec(n);
%             x2 = xvec(candidates(yvec(candidates)==bounds(n).y(1)));
%             y2 = yvec(candidates(yvec(candidates)==bounds(n).y(1)));
%             bounds(n).y(1) = @(x) (y1-y2)/(x1-x2)*x+(x1*y2-x2*y1)/(x1-x2);
%         end
%         [xnei,ixnei] = sort(xvec(candidates));
%         ynei = yvec(ixnei);
%         x4 = xnei(2); y4=ynei(2);
%         if x4>xvec(n) && y4>yvec(n)
%             x1 = xvec(n);
%             y1 = yvec(n);
%             bounds(n).x(1) = @(y) (x1-x4)/(y1-y4)*y+(y1*x4-y4*x1)/(y1-y4);
%         end
%         if bounds(n).x(2)>xvec(n)
%             
%         end
%     else
    bounds(n).x=[min(xvec(candidates)) max(xvec(candidates))];
    bounds(n).y=[min(yvec(candidates)) max(yvec(candidates))];   
%     end
end

myMesh.p=[xvec,yvec];
myMesh.dt=dt;  % DelaunayTri structure
myMesh.neighbors=neighbors;
myMesh.bounds=bounds;
myMesh.numNodes=length(myMesh.p(:,1));

base = struct('f_disc',zeros(myMesh.numNodes,2));
base = repmat(base,2*myMesh.numNodes,1);

% meshBase = struct('f_intp_x',@(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(1).f_disc(:,1),'linear'),...
%                                     'f_intp_y',@(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numNodes+1).f_disc(:,1),'linear'));
% myMesh.base = repmat(meshBase, 2*myMesh.numNodes,1);                        
myMesh.base(2*myMesh.numNodes).f_intp_x = @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numNodes).f_disc(:,1),'linear'); 
myMesh.base(2*myMesh.numNodes).f_intp_y = @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numNodes*2).f_disc(:,1),'linear'); 
%create the basis functions and interpolate them using the Delaunay Triangulation:
for j=1:myMesh.numNodes
    old_cputime = cputime;
    base(j).f_disc(j,1)=1;
    base(myMesh.numNodes+j).f_disc(j,2)=1;
    
    disp(['Creating ' num2str(j) 'th force base ... (' num2str(cputime-old_cputime) ' sec passed)'])
end
for j=1:myMesh.numNodes
    old_cputime = cputime;
    myMesh.base(j).f_intp_x= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(j).f_disc(:,1),'linear');
    myMesh.base(j).f_intp_y= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(j).f_disc(:,2),'linear'); % only zeros
    myMesh.base(j).testNumber = @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(j).f_disc(:,2),'linear',j); % only zeros
    myMesh.base(myMesh.numNodes+j).f_intp_x= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numNodes+j).f_disc(:,1),'linear'); % only zeros
    myMesh.base(myMesh.numNodes+j).f_intp_y= @(x,y) nan2zeroTriScatteredInterp(x,y,myMesh.dt,base(myMesh.numNodes+j).f_disc(:,2),'linear'); 
    disp(['Creating ' num2str(j) 'th basis function ... (' num2str(cputime-old_cputime) ' sec passed)'])
end
    function vOut=nan2zeroTriScatteredInterp(x,y,dtIn,vIn,method,j)
        if nargin<6
            F=TriScatteredInterp(dtIn,vIn,method);
            vOut=F(x,y);
            checkVec=isnan(vOut);
            vOut(checkVec)=0;
        else
            vOut = j;
        end
    end
end
% % plot an example to see if it works correctly
% ind=80;
% if length(xvec)>ind-1
%     xmin=min(xvec);
%     ymin=min(yvec);
%     xmax=max(xvec);
%     ymax=max(yvec);
%     
%     pointsPerEdge=round(sqrt(length(xvec)));
%     [x_fine,y_fine]=meshgrid(linspace(xmin,xmax,10*pointsPerEdge) , linspace(xmin,xmax,10*pointsPerEdge));
% 
%     figure(10)
%     plot(myMesh.p(myMesh.neighbors(ind).cand,1),myMesh.p(myMesh.neighbors(ind).cand,2),'or')
%     hold on
%     plot(myMesh.p(ind,1),myMesh.p(ind,2),'ob')
%     triplot(myMesh.dt);
%     plot([myMesh.bounds(ind).x(1) myMesh.bounds(ind).x(1) myMesh.bounds(ind).x(2) myMesh.bounds(ind).x(2) myMesh.bounds(ind).x(1)],[myMesh.bounds(ind).y(1) myMesh.bounds(ind).y(2) myMesh.bounds(ind).y(2) myMesh.bounds(ind).y(1) myMesh.bounds(ind).y(1)],'k')
%     quiver(x_fine,y_fine,myMesh.base(ind).f_intp_x(x_fine,y_fine),myMesh.base(ind).f_intp_y(x_fine,y_fine),'r')
%     quiver(x_fine,y_fine,myMesh.base(myMesh.numNodes+ind).f_intp_x(x_fine,y_fine),myMesh.base(myMesh.numNodes+ind).f_intp_y(x_fine,y_fine),'g')
%     xlim([xmin xmax])
%     ylim([ymin ymax])
%     hold off
% end