function [myMesh]=createMeshAndBasisFromAdhesions(x_vec,y_vec,paxImage,pixelSize)
% createMeshAndBasisFromAdhesions creates mesh from adhesion information in
% paxImage. x_vec and y_vec serve as mask for mesh nodes.
% Sangyoon Han June 2013
disp('Creating force mesh with paxillin image ...')
paxImage = double(paxImage);
% Get Cell mask
pId = paxImage/max(paxImage(:));
alpha = graythresh(pId);
bwPI = im2bw(pId,alpha);
bwPI2 = bwmorph(bwPI,'clean');

% Get FA mask
minSize = round((1000/pixelSize)*(200/pixelSize)); %adhesion limit=1um*.5um
maskFA = blobSegmentThreshold(paxImage,minSize,false,bwPI2);

% For each FA, get points along long axis

% find adhesions from paxImage
sigmaPSF_NA = 1.48;
[pstruct_NA,mask] = pointSourceDetection(paxImage, sigmaPSF_NA, 'alpha', 0.05);
NAs = [ceil(pstruct_NA.x') ceil(pstruct_NA.y')];

% To get high-resolution information, subsample detected NAs ensuring 
% NAs are separated by 2 pixel 
disp('Subsampling detected NAs (high resolution)...')
idx = KDTreeBallQuery(NAs, NAs, 2);
valid = true(numel(idx),1);
for i = 1:numel(idx)
    if ~valid(i), continue; end
    neighbors = idx{i}(idx{i}~=i);
    valid(neighbors) = false;
end
NAs = NAs(valid, :);

% Add nodes in the sparse area

%delaunay_mesh=delaunay(x_vec,y_vec);
dt=DelaunayTri(x_vec,y_vec);
delaunay_mesh=dt.Triangulation;

%for each node find all its neighbors:
for n=1:length(x_vec)
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

    bounds(n).x=[min(x_vec(candidates)) max(x_vec(candidates))];
    bounds(n).y=[min(y_vec(candidates)) max(y_vec(candidates))];   
end

myMesh.p=[x_vec, y_vec];
myMesh.dt=dt;  % DelaunayTri structure
myMesh.neighbors=neighbors;
myMesh.bounds=bounds;
myMesh.numNodes=length(myMesh.p(:,1));

%create the basis functions and interpolate them using the Delaunay Triangulation:
for j=1:myMesh.numNodes
    base(j).f_disc=zeros(myMesh.numNodes,2);
    base(j).f_disc(j,1)=1;
    
    myMesh.base(j).f_intp_x= TriScatteredInterp(myMesh.dt,base(j).f_disc(:,1),'linear');
    myMesh.base(j).f_intp_y= TriScatteredInterp(myMesh.dt,base(j).f_disc(:,2),'linear'); % only zeros
    
    base(myMesh.numNodes+j).f_disc=zeros(myMesh.numNodes,2);
    base(myMesh.numNodes+j).f_disc(j,2)=1;
    
    myMesh.base(myMesh.numNodes+j).f_intp_x= TriScatteredInterp(myMesh.dt,base(myMesh.numNodes+j).f_disc(:,1),'linear'); % only zeros
    myMesh.base(myMesh.numNodes+j).f_intp_y= TriScatteredInterp(myMesh.dt,base(myMesh.numNodes+j).f_disc(:,2),'linear'); 
end

% plot an example to see if it works correctly
ind=80;
if length(x_vec)>ind-1
    xmin=min(x_vec);
    ymin=min(y_vec);
    xmax=max(x_vec);
    ymax=max(y_vec);
    
    pointsPerEdge=round(sqrt(length(x_vec)));
    [x_fine y_fine]=meshgrid(linspace(xmin,xmax,10*pointsPerEdge) , linspace(xmin,xmax,10*pointsPerEdge));

    figure(10)
    plot(myMesh.p(myMesh.neighbors(ind).cand,1),myMesh.p(myMesh.neighbors(ind).cand,2),'or')
    hold on
    plot(myMesh.p(ind,1),myMesh.p(ind,2),'ob')
    triplot(myMesh.dt);
    plot([myMesh.bounds(ind).x(1) myMesh.bounds(ind).x(1) myMesh.bounds(ind).x(2) myMesh.bounds(ind).x(2) myMesh.bounds(ind).x(1)],[myMesh.bounds(ind).y(1) myMesh.bounds(ind).y(2) myMesh.bounds(ind).y(2) myMesh.bounds(ind).y(1) myMesh.bounds(ind).y(1)],'k')
    quiver(x_fine,y_fine,myMesh.base(ind).f_intp_x(x_fine,y_fine),myMesh.base(ind).f_intp_y(x_fine,y_fine),'r')
    quiver(x_fine,y_fine,myMesh.base(myMesh.numNodes+ind).f_intp_x(x_fine,y_fine),myMesh.base(myMesh.numNodes+ind).f_intp_y(x_fine,y_fine),'g')
    xlim([xmin xmax])
    ylim([ymin ymax])
    hold off
end