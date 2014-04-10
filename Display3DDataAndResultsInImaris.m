function [ varargout ] = Display3DDataAndResultsInImaris( im, spacing, varargin )

    p = inputParser;
    p.addRequired( 'im', @(x) (isnumeric(x) && ndims(x) == 3) );
    p.addRequired( 'spacing', @(x) (isnumeric(x) && numel(x) == 3) );
    p.parse( im, spacing );
    
    p.addParamValue( 'rgbMask', [], @(x) (isnumeric(x) && ndims(x) == 4 && ~any(size(x)~=[size(im), 3])) );
    p.addParamValue( 'spotLocations', [], @(x) (isnumeric(x) && ndims(x) == 2 && size(x,2) == 3) );
    p.addParamValue( 'spotColor', [1, 0, 0], @(x) (isnumeric(x) && numel(x) == 3) );
    p.addParamValue( 'spotRadius', 3, @(x) (isnumeric(x) && isscalar(x)) );
    p.addParamValue( 'surfaceObjects', [], @(x) (iscell(x) || isstruct(x)) );
    p.parse( im, spacing, varargin{:} );
    
    imRGBMask = p.Results.rgbMask;
    spotLocations = p.Results.spotLocations;
    spotColor = p.Results.spotColor;
    spotRadius = p.Results.spotRadius;
    surfaceObjects = p.Results.surfaceObjects;
    
    % initialization
    imarisApp = imarisStartNew(nargout==0);
    imarisScene = imarisApp.mFactory.CreateDataContainer;
    imarisApp.mSurpassScene = imarisScene;
    imarisScene.AddChild(imarisApp.mFactory.CreateLightSource); %add the light to the scene
    imarisScene.AddChild(imarisApp.mFactory.CreateFrame); %add the frame to the scene    
    
    volSize = size( im );
    volData = imarisApp.mFactory.CreateDataSet;
    
    if ~isempty(imRGBMask)
        volData.Create( 'eTypeUint16', volSize(2), volSize(1),  volSize(3), 4, 1);
    else
        volData.Create( 'eTypeUint16', volSize(2), volSize(1),  volSize(3), 1, 1);
    end
    
    volData.mExtendMinX = 0;
    volData.mExtendMinY = 0;
    volData.mExtendMinZ = 0;
    volData.mExtendMaxX = volData.mSizeX * spacing(1);
    volData.mExtendMaxY = volData.mSizeY * spacing(2);
    volData.mExtendMaxZ = volData.mSizeZ * spacing(3);
    volData.mUnit = 'um';
    
    ImageIntensityRange = ComputeImageDynamicRange( im, 99.0 );
    volData.SetDataVolume( uint16( permute(flipdim(im, 1), [2,1,3]) ), 0, 0);
    volData.SetChannelColor( 0, 1, 1, 1, 0.8 );
    volData.SetChannelRange( 0, ImageIntensityRange(1), ImageIntensityRange(2) );
    
    % add rgb mask
    if ~isempty(imRGBMask)        
        channelColor = [1, 0, 0; 0 1 0; 0 0 1];
        imRGBMask = permute(flipdim(imRGBMask, 1), [2,1,3,4]);
        imRGBMask = uint16( mat2gray( imRGBMask ) * 255.0 );
        for i = 1:3
            volData.SetDataVolume( imRGBMask(:,:,:,i), i, 0);
            volData.SetChannelColor( i, channelColor(i,1), channelColor(i,2), channelColor(i,3), 0.8 );
            volData.SetChannelRange( i, 0, 255 );
        end
    end
    
    imarisApp.mDataSet = volData;
    imarisApp.mSurpassCamera.Fit;
    
    %Strangely enough, this is both necessary and sufficient to display the mDataSet in the newly created scene.
    imarisScene.AddChild(imarisApp.mFactory.CreateVolume);
    
    % add spots
    if ~isempty(spotLocations)
        
        numSpots = size(spotLocations, 1);
        
        spotLocations(:,2) = volSize(1) - spotLocations(:,2) + 1; % flip-y
        spotLocationsPhysp = (spotLocations - 1) .* repmat(spacing([2,1,3]), [numSpots, 1] );        
        spotRadii = ones(numSpots, 1) * spotRadius * spacing(1);
        spotTimes = zeros(numSpots,1);

        imarisSpots = imarisApp.mFactory.CreateSpots;
        imarisSpots.SetColor(spotColor(1), spotColor(2), spotColor(3), 0 );
        imarisSpots.Set(spotLocationsPhysp, spotTimes, spotRadii);                
        imarisApp.mSurpassScene.AddChild(imarisSpots);
                
    end
    
    % add surfaces
    if ~isempty( surfaceObjects )
       
        if isstruct( surfaceObjects )
            surfaceObjects = { surfaceObjects };
        end
        
        for obid = 1:numel(surfaceObjects)
            
            if ~isstruct(surfaceObjects{obid}) || ~isfield(surfaceObjects{obid}, 'surfaces')
                error( 'ERROR: each surface object must be a structure/structure-array contain a field surfaces specifying the geometry (vertices, faces, normals) for each surface belonging to the surface object' );
            end
            
            if sum(isfield( surfaceObjects{obid}.surfaces, {'vertices', 'faces', 'normals'})) ~= 3
                error( 'ERROR: each surface belonging to a suface object must contain must contain three fields -- vertices, faces, normals --- which degine the geometry of the surface' );
            end
                        
            curImarisSurfaceObject = imarisApp.mFactory.CreateSurfaces;
            
            % set name
            if isfield( surfaceObjects{obid}, 'name' )
                curImarisSurfaceObject.mName = surfaceObjects{obid}.name;
            end
            
            % set color
            if isfield( surfaceObjects{obid}, 'color' )
                if numel( surfaceObjects{obid}.color ) <= 3
                    surfaceObjects{obid}.color(4) = 0.0;
                end
                curImarisSurfaceObject.SetColor( surfaceObjects{obid}.color(1), ...
                                                 surfaceObjects{obid}.color(2), ...
                                                 surfaceObjects{obid}.color(3), ...
                                                 surfaceObjects{obid}.color(4) );
            end
           
            for sid = 1:numel( surfaceObjects{obid} )

                curSurface = surfaceObjects{obid}.surfaces(sid);
                
                % flip y-dim of vertices
                curSurface.vertices(:,2) = volSize(1) - curSurface.vertices(:,2) + 1;

                % convert vertex coordinates to physical space
                curSurface.vertices = (curSurface.vertices - 1) .* repmat(spacing([2,1,3]), [size(curSurface.vertices,1), 1] );        

                % convert face indices to 0-based index
                curSurface.faces = curSurface.faces - 1;
           
                % time index
                if ~isfield( curSurface, 'timepoint' )
                    curSurface.timepoint = 0;
                end
                
                % add surface
                curImarisSurfaceObject.AddSurface( curSurface.vertices, ...
                                                   curSurface.faces, ...
                                                   curSurface.normals, ...
                                                   curSurface.timepoint );
                
            end
            
            imarisApp.mSurpassScene.AddChild(curImarisSurfaceObject);
            
        end
        
    end
    
    % return imaris app object if requested
    if nargout >= 1
        varargout{1} = imarisApp;
    end
    
end
