function [ varargout ] = Display3DVideoAndResultsInImaris( videoInput, varargin )

    % default color map -- only 6 channels are defined
    % if more than 6 channels then user should specify the displayColors 
    cMap(1,:) = [1 0 0];
    cMap(2,:) = [0 1 0];
    cMap(3,:) = [0 0 1];
    cMap(4,:) = [1 0 1];
    cMap(5,:) = [0 1 1];
    cMap(6,:) = [1 1 0];

    % get and validate input video argument
    p = inputParser;
    p.addRequired( 'videoInput', @(x) (isnumeric(x) && ismember(ndims(x), [4,5])) );
    p.parse( videoInput );

    if ndims(videoInput) == 4
       videoSize = size(videoInput);    
       videoInput = reshape(videoInput, [videoSize(1:3), 1, videoSize(4)]);
    end
    
    videoSize = size(videoInput);
    numTimePoints = videoSize(5);
    numChannels = videoSize(4);

    default_displayranges = zeros(numChannels,2);
    for i = 1:numChannels
       imCurChannel = videoInput(:,:,:,i,1);
       default_displayranges(i,:) = double( [ min(imCurChannel(:)) max(imCurChannel(:)) ] );
    end
    
    % get and validate other optional input arguments if provided
    p.addParamValue( 'spacing', ones(1,3), @(x) (isnumeric(x) && numel(x) == 3) );
    p.addParamValue( 'units', 'um', @(x) (ischar(x)) );
    p.addParamValue( 'voxelType', 'eTypeUInt16', @(x) (ischar(x) && ismember(x, {'eTypeUInt8', 'eTypeUInt16', 'eTypeFloat'})) );
    p.addParamValue( 'displayColors', cMap(1:numChannels,:), @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 3 && size(x,1) == numChannels ) );
    p.addParamValue( 'displayRanges', default_displayranges, @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 2 && size(x,1) == numChannels ) );
    p.addParamValue( 'spots', [], @(x) (iscell(x) || isstruct(x)) );
    p.addParamValue( 'surfaceObjectContainers', [], @(x) (iscell(x) || isstruct(x)) );
    p.parse( videoInput, varargin{:} );
    
    PARAMETERS = p.Results;
    spacing = PARAMETERS.spacing;
    displayColors = PARAMETERS.displayColors;
    displayRanges = PARAMETERS.displayRanges;
    units = PARAMETERS.units;
    voxelType = PARAMETERS.voxelType;
    spots = PARAMETERS.spots;
    surfaceObjectContainers = p.Results.surfaceObjectContainers;
    
    % initialization
    imarisApp = imarisStartNew(nargout==0);
    imarisScene = imarisApp.mFactory.CreateDataContainer;
    imarisApp.mSurpassScene = imarisScene;
    imarisScene.AddChild(imarisApp.mFactory.CreateLightSource); %add the light to the scene
    imarisScene.AddChild(imarisApp.mFactory.CreateFrame); %add the frame to the scene    

    % create a dataset for the multichannel 3D video
    volData = imarisApp.mFactory.CreateDataSet;
    volData.Create(voxelType, videoSize(2), videoSize(1),  videoSize(3), numChannels, numTimePoints);

    volData.mExtendMinX = 0;
    volData.mExtendMinY = 0;
    volData.mExtendMinZ = 0;
    volData.mExtendMaxX = volData.mSizeX * spacing(2);
    volData.mExtendMaxY = volData.mSizeY * spacing(1);
    volData.mExtendMaxZ = volData.mSizeZ * spacing(3);
    volData.mUnit = units;    
    
    for tid = 1:numTimePoints        
        for cid = 1:numChannels 

            volData.SetDataVolume( typeCastVolume( permute(flipdim(videoInput(:,:,:,cid,tid), 1), [2,1,3]), voxelType ), cid-1, tid-1);        
            volData.SetChannelColor( cid-1, displayColors(cid,1), displayColors(cid,2), displayColors(cid,3), 0.8 );
            volData.SetChannelRange( cid-1, displayRanges(cid,1), displayRanges(cid,2) );

        end        
    end
    
    imarisApp.mDataSet = volData;
    
    % Display the dataset in the newly created scene.
    imarisApp.mSurpassCamera.Fit;
    imarisScene.AddChild(imarisApp.mFactory.CreateVolume);    
    
    % add spots
    if ~isempty(spots)

        spotGroup = imarisApp.GetFactory.CreateDataContainer;
        spotGroup.SetName('Spots');
    
        if isstruct(spots) 
            spots = { spots };
        end
        
        for spid = 1:numel(spots)

            curSpotObject = spots{spid};
            
            if ~isstruct(curSpotObject) || ~all(isfield(curSpotObject, {'locations', 'timepoints'}))
                error( 'ERROR: each spot object must contain the fields -- locations, timepoints -- specifying the location of the spots in image space and time' );
            end
            
            % get spot locations
            if ~ismatrix(curSpotObject.locations) || size(curSpotObject.locations, 2) ~= 3
                error( 'ERROR: the locations field of the spot object must be an n x 3 matrix where n is the number of spots and each row specified the (x,y,z) location of the spot in image space' );
            end
                
            spotLocations = curSpotObject.locations;
            numSpots = size(spotLocations, 1);
            spotLocations(:,2) = volSize(1) - spotLocations(:,2) + 1; % flip-y
            spotLocationsPhysp = (spotLocations - 1) .* repmat(spacing([2,1,3]), [numSpots, 1] );        

            % get spot times
            if isscalar(curSpotObject.timepoints)
                spotTimes = curSpotObject.timepoints + zeros(numSpots,1);
            else
                spotTimes = curSpotObject.timepoints;
                if ~isvector(spotTimes) || numel(spotTimes) ~= numSpots
                    error( 'ERROR: the timepoints field should either be a scalar or a vector of size equal to the number of spots in the locations field' );
                end
            end
            
            % get spot radii -- assumed to be in physical space
            if ~isfield(curSpotObject, 'radii')
                spotRadii = ones(numSpots, 1) * 3 * min(spacing);
            else
                spotRadii = curSpotObject.radii;
                if isscalar(spotRadii)
                   spotRadii = spotRadii + zeros(numSpots, 1); 
                end
                
                if ~isvector(spotRadii) || numel(spotRadii) ~= numSpots
                    error( 'ERROR: the optional radii field should either be a scalar or a vector of size equal to the number of spots in the locations field' );
                end
            end
            
            % track edges
            if isfield(curSpotObject, 'trackEdges')
                
                if ischar(curSpotObject.trackEdges, 'consecutive')
                    trackEdges = ([0:numSpots-1; 1:numSpots])';
                else
                    trackEdges = curSpotObject.trackEdges;
                end
                
            end
            
            % get spot color
            if ~isfield(curSpotObject, 'color')
                spotColor = [1, 0, 0];
            else
                spotColor = curSpotObject.spotColor;
                if ~isvector(spotColor) || numel(spotColor) ~= 3
                    error( 'ERROR: the optional color field should a 3-element vector specifying the red, green and blue values' );
                end
            end
            
            % create imaris spot object
            imarisSpots = imarisApp.mFactory.CreateSpots;
            imarisSpots.SetColor(spotColor(1), spotColor(2), spotColor(3), 0 );
            imarisSpots.Set(spotLocationsPhysp, spotTimes, spotRadii);                
            imarisSpots.SetTrackEdges(trackEdges);
            
            % set name
            if isfield( curSpotObject, 'name' )
                imarisSpots.mName = curSpotObject.name;
            end
            
            % add to spot group
            spotGroup.AddChild(imarisSpots, -1);
            
        end
           
        % add spot objects container to imaris surpass scene
        imarisApp.GetSurpassScene.AddChild(spotGroup, -1);    
        
    end    
    
    % add surfaces
    if ~isempty( surfaceObjectContainers )
       
        surfaceGroup = imarisApp.GetFactory.CreateDataContainer;
        surfaceGroup.SetName('Surfaces');
        
        if isstruct( surfaceObjectContainers )
            surfaceObjectContainers = { surfaceObjectContainers };
        end
        
        for gid = 1:numel(surfaceObjectContainers)
            
            if ~isstruct(surfaceObjectContainers{gid}) || ~isfield(surfaceObjectContainers{gid}, 'surfaces')
                error( 'ERROR: each surface object must be a structure/structure-array contain a field surfaces specifying the geometry (vertices, faces, normals) for each surface belonging to the surface object' );
            end
            
            if sum(isfield( surfaceObjectContainers{gid}.surfaces, {'vertices', 'faces', 'normals', 'timepoint'})) ~= 4
                error( 'ERROR: each surface belonging to a suface object must contain must contain three fields -- vertices, faces, normals, timepoint --- which degine the geometry of the surface' );
            end
                        
            curImarisSurfaceObject = imarisApp.mFactory.CreateSurfaces;
            
            % set name
            if isfield( surfaceObjectContainers{gid}, 'name' )
                curImarisSurfaceObject.mName = surfaceObjectContainers{gid}.name;
            end
            
            % set color
            if isfield( surfaceObjectContainers{gid}, 'color' )
                if numel( surfaceObjectContainers{gid}.color ) <= 3
                    surfaceObjectContainers{gid}.color(4) = 0.0;
                end
                curImarisSurfaceObject.SetColor( surfaceObjectContainers{gid}.color(1), ...
                                                 surfaceObjectContainers{gid}.color(2), ...
                                                 surfaceObjectContainers{gid}.color(3), ...
                                                 surfaceObjectContainers{gid}.color(4) );
            end
           
            for sid = 1:numel( surfaceObjectContainers{gid} )

                curSurface = surfaceObjectContainers{gid}.surfaces(sid);
                
                % flip y-dim of vertices
                curSurface.vertices(:,2) = volSize(1) - curSurface.vertices(:,2) + 1;

                % convert vertex coordinates to physical space
                curSurface.vertices = (curSurface.vertices - 1) .* repmat(spacing([2,1,3]), [size(curSurface.vertices,1), 1] );        

                % convert face indices to 0-based index
                curSurface.faces = curSurface.faces - 1;
           
                % timepoint
                %curSurface.timepoint =  curSurface.timepoint - 1; % 0-based?
                
                % add surface
                curImarisSurfaceObject.AddSurface( curSurface.vertices, ...
                                                   curSurface.faces, ...
                                                   curSurface.normals, ...
                                                   curSurface.timepoint );
                
            end
            
            surfaceGroup.AddChild(curImarisSurfaceObject, -1);
            
        end

        imarisApp.mSurpassScene.AddChild(surfaceGroup, -1);
        
    end
    
    % return imaris app object if requested
    if nargout >= 1
        varargout{1} = imarisApp;
    end
    
end

function [vout] = typeCastVolume( vin, voxelType )

    switch voxelType
    
        case 'eTypeUInt8'            
            vout = uint8( vin );
            
        case 'eTypeUInt16'
            
            vout = uint16( vin );
            
        case 'eTypeFloat'
        
            vout = single( vin );
    end
    
end