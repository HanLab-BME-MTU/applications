function makiFitPlane(dataStruct)
%MAKIFITPLANE attempts to fit planes into kinetochore clusters
%
% SYNOPSIS: makiFitPlane(dataStruct)
%
% INPUT dataStruct (opt): maki data structure. If empty, it will be loaded
%                         via GUI
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 03-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Later: Give coordinates (in um) as input, rather than dataStruct, so that
% it is clear which kind of coordinates that should be read

%TEST input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

% get coordinates, dataProperties, etc
coordStruct = dataStruct.initCoord;
nTimepoints = length(coordStruct);

% loop through timepoints. Get covariance of point cloud, and the
% corresponding eigenvalues. Compare the two most similar eigenvalues to
% the third to find out how much of a plane there is
coordCov = zeros(3,3,nTimepoints);
eigenValues = zeros(nTimepoints, 3);
eigenVectors = zeros(3,3,nTimepoints);  %vectors in cols

eigenRatio = zeros(nTimepoints,8);

coordCellMu = cell(nTimepoints,1);
nSpots = cat(1,coordStruct.nSpots);

for t=1:nTimepoints
    
    % leaveCoordCellMu as a cell to avoid rewriting everything
    coordCellMu{t} = coordStruct(t).allCoord(:,1:3);
    
    % get covariance
    coordCov(:,:,t) = cov(coordCellMu{t});
    
    % eigenVectors/eigenValues
    [eigenVectors(:,:,t),eVals] = eig(coordCov(:,:,t));
    eigenValues(t,:) = diag(eVals);
    
    % compare eigenValues
    diffs = pdist(eigenValues(t,:)');
    % indices to understand lowest
    [u,v] = find(tril(ones(3),-1));
    [dummy,idx] = min(diffs);
    % find two close, one far
    closestIdx = [u(idx),v(idx)];
    farIdx = setdiff(1:3,closestIdx);
    
    % take ratio
    closeAvg = mean(eigenValues(t,closestIdx));
    % since we want the rest for comparison for later, we don't need to
    % divide yet
    eigenRatio(t,:) = [0,eigenValues(t,farIdx),closeAvg,eigenValues(t,closestIdx),farIdx,closestIdx];
    
end

eigenRatio(:,1) = eigenRatio(:,2)./eigenRatio(:,3);
figure,plot(1:nTimepoints,eigenRatio(:,1),'-ro',...
    1:nTimepoints,eigenRatio(:,2),'--m',...
    1:nTimepoints,eigenRatio(:,3),'--b+',...
    1:nTimepoints,eigenRatio(:,4),'.b',...
    1:nTimepoints,eigenRatio(:,5),'.b');

% goodFrames are those with eigenRatio<1 (for now)
goodFramesL = (eigenRatio(:,1)<1);
% often, the eigenRatio starts low and then peaks before going back down
% find this by using bwlabel and identify the largest group
goodFramesLb = bwlabel(goodFramesL);
[idx,cts] = countEntries(goodFramesLb);
goodLabel = idx(find(cts(2:end)==max(cts(2:end)))+1);
goodFrames = find(goodFramesLb==goodLabel);

% loop goodFrames, fit plane
planeFit(1:nTimepoints) = struct('plane',[],'goodFrames',[],'planeCoord',[]);
% store goodFrames in pf1
planeFit(1).goodFrames = goodFrames;

for t = goodFrames'
    % fit plane: ax+by+cz = d (or ax+by+cz-d=0)
    % [a,b,c] is the normal to the plane, d is a function of the distance
    % from zero
    % Use the eigenvector as an initial guess for abc, and the take the
    % average coordinate to get d as dot([abc][xmymzm]). After fitting,
    % normalize abc to 1.
    % fit using linearLeastMedianSquares
%     A = coordCellMu{t};
%     x0 = eigenVectors(:,eigenRatio(t,6),t);
%     B = repmat(sum(mean(A,1).*x0',2),nSpots(t),1);
%     [x,qxx,goodRows]=linearLeastMedianSquares(A,B,[],x0);
%     badRows = setdiff(1:nSpots(t),goodRows);
%     % norm
%     n = sqrt(sum(x.^2));
%     plane = [x'./n B(1)./n];
%     
% normal = plane(1:3);
% [xgrid,ygrid] = meshgrid(linspace(min(A(:,1)),max(A(:,1)),5), ...
%     linspace(min(A(:,2)),max(A(:,2)),5));
% zgrid1 = (1/normal(3)) .* (plane(4) - (xgrid.*normal(1) + ygrid.*normal(2)));

% do not fit. Somehow, the eigenVectors give a more stable solution. 

meanCoord = mean(coordCellMu{t},1);
normal = eigenVectors(:,eigenRatio(t,6),t);
planeFit(t).plane = [normal',meanCoord*normal];
% A = coordCellMu{t};
% [xgrid,ygrid] = meshgrid(linspace(min(A(:,1)),max(A(:,1)),5), ...
%     linspace(min(A(:,2)),max(A(:,2)),5));
% zgrid2 = (1/normal(3)) .* (meanCoord*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
% 
% figure,
% h = mesh(xgrid,ygrid,zgrid1,'EdgeColor',[0 1 0],'FaceAlpha',0);
% hold on
% h = mesh(xgrid,ygrid,zgrid2,'EdgeColor',[0 0 1],'FaceAlpha',0);
% plot3(A(:,1),A(:,2),A(:,3),'.r')
% plot3(A(badRows,1),A(badRows,2),A(badRows,3),'ok')

% transform all spots into plane coordinate system. In case we fit at some
% point, calculate coords from the normal only
e_plane = zeros(3);
e_plane(:,1) = normal;
e_plane(:,2) = [-normal(2),normal(1),0]./sqrt(sum(normal(1:2).^2));
e_plane(:,3) = cross(e_plane(:,1),e_plane(:,2));
% [d,xPlane,yPlane]'=[e_plane]^-1*(coord-origin)
% d: distance, xyPlane: coordinates in the plane
% origin: meanCoord (or [0,0,d/c])
planeFit(t).planeCoord = (inv(e_plane)*(coordCellMu{t}-repmat(meanCoord,nSpots(t),1))')';

end
1
% Plotting commands
% plot plane in matlab
% pc=planeFit(45).planeCoord;
% pos = pc(:,1)>0;
% neg = pc(:,1)<0;
% figure,plot3(pc(pos,1),pc(pos,2),pc(pos,3),'.k',pc(neg,1),pc(neg,2),pc(neg,3),'or')
% 
% c=-5.75:0.5:5.75;
% nBins = length(c);
% z = zeros(nBins,nTimepoints);
% for t=goodFrames',z(:,t)=hist(planeFit(t).planeCoord(:,1),c);end
% uiViewPanel,imshow(z,[]),colormap(cm)




