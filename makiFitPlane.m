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
coordCell = dataStruct.initCoord;
pixelSize = [dataStruct.dataProperties.PIXELSIZE_XY,...
    dataStruct.dataProperties.PIXELSIZE_XY,...
    dataStruct.dataProperties.PIXELSIZE_Z];
nTimepoints = length(coordCell);

% loop through timepoints. Get covariance of point cloud, and the
% corresponding eigenvalues. Compare the two most similar eigenvalues to
% the third to find out how much of a plane there is
coordCov = zeros(3,3,nTimepoints);
eigenValues = zeros(nTimepoints, 3);
eigenVectors = zeros(3,3,nTimepoints);  %vectors in cols

eigenRatio = zeros(nTimepoints,3);

coordCellMu = cell(nTimepoints,1);
nSpots = zeros(nTimepoints,1);

for t=1:nTimepoints
    
    % calculate coords in microns
    nSpots(t) = length(coordCell{t});
    coordCellMu{t} = coordCell{t}(:,1:3).*repmat(pixelSize,nSpots(t),1);
    
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
    eigenRatio(t,:) = [0,eigenValues(t,farIdx),closeAvg];
    
end

eigenRatio(:,1) = eigenRatio(:,2)./eigenRatio(:,3);
figure,plot(eigenRatio(:,1));
