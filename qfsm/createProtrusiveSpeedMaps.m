function protSpeedMap=createProtrusiveSpeedMaps(Md,mask,pixelSize,sampling,gridSize,varargin)
% fsmSpeedMaps creates speed maps from the flow maps
%
% createProtrusiveSpeedMaps goes through the whole M (or Md) and calculate
% V dot n where n is the normal vector of the edge of cell closest to
% individual vectors.
%   
% INPUT         Md          : the interpolated track matrix
%               maskProc    : mask process
%               pixelSize : pixel size in the image domain (nm)
%               imgSize   : pixel size in the image domain (pixels)
%               gridSize

% OUTPUT        protSpeedMap  : a cell array of matrices of size imgSize
%
%
% Sangyoon Han 2022


% Check input
ip = inputParser;
ip.addRequired('Md',@iscell);
ip.addRequired('mask',@islogical);
ip.addRequired('pixelSize',@isscalar);
ip.addRequired('sampling',@isscalar);
ip.addRequired('gridSize',@isscalar);
ip.parse(Md,mask,pixelSize,sampling,gridSize,varargin{:});

samBatch = false; %Wether to run Sam's function in batch mode or not.
splineTolerance = 30;
downSample = 50;

% We go with whatever frames where Md is non-zero.
startFrames = 1:numel(Md);
for i=startFrames
    
    % velocities 
    vel = Md{i}(:,4:-1:3)-Md{i}(:,2:-1:1);
    pos = Md{i}(:,2:-1:1);
    
    % surface normal vectors are too detailed. will make it smooth to remove
    % the details.
    curMask = mask(:,:,i);
    curMask = bwmorph(curMask,'close',20);
    curMask = bwmorph(curMask,'erode',5);
    smoothMask = filterGauss2D(curMask,10)>5e-1; %conv2(double(mask(:,:,i)), ones(5)/25, 'same'); %

    % Get closest edge - adopting from prSamProtrusion
    [currMask, isCurrClosed, currOutline] = ...
        getOutlineInfo(smoothMask, i);
    
    %Set up the inputs for sam's function
    samIn = struct();
    samIn.batch_processing = samBatch;
    samIn.TOLERANCE = splineTolerance;
    samIn.dl_rate = downSample;
    
    %Set up input for sam's function
    samIn.pixel_list_last = currOutline;
    samIn.pixel_list = currOutline;
    samIn.ISCLOSE = isCurrClosed;
    samIn.maskTM1 = currMask;
    samIn.time_idx = i;
    
    %----call sam's protrusion calc function----%
    samOut = prSamProtrusion(samIn);
        
    %Got the normal
    normals = samOut.xy_normal;
    magNormal = (normals(:,1).^2+normals(:,2).^2).^0.5;
    unitNormals = normals./magNormal;
    
    % Get the index of points at closest edge from vector positions
    [idx,~]=KDTreeClosestPoint(samOut.pixel_tm1_output,pos);
    corrNormVecs = unitNormals(idx,:);  
    
    % Calculating V*n (V=vel, n=corrNormVecs)
    protSpeed = dot(vel,corrNormVecs,2);
    
    % Reshape
    iMap=i;
    [X,Y]=meshgrid(floor(min(pos(:,1))):floor(max(pos(:,1))),floor(min(pos(:,2))):floor(max(pos(:,2))));
%     [X,Y]=meshgrid(1:size(currMask,2),1:size(currMask,1));
%     curProtSpeedMap = griddata(pos(:,1),pos(:,2),protSpeed,X,Y,'cubic');
    F = scatteredInterpolant(pos(:,1),pos(:,2),protSpeed,'linear','linear');
    curProtSpeedMap = F(X,Y);
%     curProtSpeedMap(isnan(curProtSpeedMap)) = 0;
    curProtSpeedMapWhole = NaN(size(currMask));
    curProtSpeedMapWhole(floor(min(pos(:,2))):floor(max(pos(:,2))),floor(min(pos(:,1))):floor(max(pos(:,1)))) = curProtSpeedMap;
    
    % Transform to nm/min
    curProtSpeedMapWhole=curProtSpeedMapWhole*(60/sampling)*pixelSize;
    
    % Use union of masks
    nanMask = double(mask(:,:,i));
    nanMask(nanMask==0) = NaN;
    protSpeedMap{iMap} = curProtSpeedMapWhole .* nanMask;
end
