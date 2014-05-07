classdef ImarisDataVisualizer

    methods (Access = public)
        
        function this = ImarisDataVisualizer( dataIn, varargin )

            % get and validate input arguments
            p = inputParser;
            p.CaseSensitive = false;
            p.addRequired( 'videoInput', @(x) (iscell(x) || (isnumeric(x) && ismember(ndims(x), [3,4,5]))) );
            p.addParamValue( 'format', '3DC', @(x) (ischar(x) && ismember(x, {'3DC', '3DT'})));
            p.addParamValue( 'spacing', ones(1,3), @(x) (isnumeric(x) && numel(x) == 3) );
            p.addParamValue( 'units', 'um', @(x) (ischar(x)) );
            p.addParamValue( 'voxelType', 'eTypeUInt16', @(x) (ischar(x) && ismember(x, {'eTypeUInt8', 'eTypeUInt16', 'eTypeFloat'})) );
            p.addParamValue( 'displayColors', [], @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 3) );
            p.addParamValue( 'displayRanges', [], @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 2) );
            p.parse( dataIn, varargin{:} );
            
            PARAMETERS = p.Results;
            
            if iscell(dataIn)                
                
                if isvector(dataIn)
                    % Input is provided as a cell array. 
                    % Assumed that each element of cell array corresponds to one time point
                    dataIn = cat(5, dataIn{:});
                else                   
                    % Input is provided as a cell matrix. 
                    % Assumed that each row of cell array corresponds to one time point
                    % Assumed that each column corresponds to a channel 
                    im = cell(size(dataIn, 1), 1);
                    for i = 1:size(dataIn, 1)
                        im{i} = cat(4, dataIn{i,:});
                    end
                    im = cat(5, im{:});
                    dataIn = im;
                    clear im;
                end

                if ndims(dataIn) == 4
                   videoSize = size(dataIn);    
                   dataIn = reshape(dataIn, [videoSize(1:3), 1, videoSize(4)]);
                end
            else
                if ndims(dataIn) == 4 && strcmpi(PARAMETERS.format, '3DT')
                    videoSize = size(dataIn);    
                    dataIn = reshape(dataIn, [videoSize(1:3), 1, videoSize(4)]);
                end
            end

            videoSize = size(dataIn);
            videoSize(numel(videoSize)+1:5) = 1;

            numTimePoints = videoSize(5);
            numChannels = videoSize(4);

            spacing = PARAMETERS.spacing;
            units = PARAMETERS.units;
            voxelType = PARAMETERS.voxelType;
            
            if isempty(PARAMETERS.displayRanges)
                default_displayranges = zeros(numChannels,2);
                for i = 1:numChannels
                   imCurChannel = dataIn(:,:,:,i,1);
                   default_displayranges(i,:) = double( [ min(imCurChannel(:)) max(imCurChannel(:)) ] );
                end
                displayRanges = default_displayranges;
            else
                displayRanges = PARAMETERS.displayRanges;
            end

            if isempty(PARAMETERS.displayColors)
                
                % default color map -- only 6 channels are defined
                % if more than 6 channels then rest are assigned randomly
                cMap(1,:) = [1 0 0];
                cMap(2,:) = [0 1 0];
                cMap(3,:) = [0 0 1];
                cMap(4,:) = [1 0 1];
                cMap(5,:) = [0 1 1];
                cMap(6,:) = [1 1 0];
                
                if numChannels <= 6
                    displayColors = cMap(1:numChannels,:);
                else
                    displayColors = [cMap; random('Uniform', 0.00001, 1, [numChannels-6, 3])];
                end
                
            else
                displayColors = PARAMETERS.displayColors;
            end
            
            % initialization
            this.imarisApp = imarisStartNew(nargout==0);
            imarisScene = this.imarisApp.mFactory.CreateDataContainer;
            this.imarisApp.mSurpassScene = imarisScene;
            imarisScene.AddChild(this.imarisApp.mFactory.CreateLightSource); %add the light to the scene
            imarisScene.AddChild(this.imarisApp.mFactory.CreateFrame); %add the frame to the scene    
            
            % create a dataset for the multichannel 3D video
            volData = this.imarisApp.mFactory.CreateDataSet;
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

                    volData.SetDataVolume( this.typeCastVolume( permute(flipdim(dataIn(:,:,:,cid,tid), 1), [2,1,3]), voxelType ), cid-1, tid-1);        
                    volData.SetChannelColor( cid-1, displayColors(cid,1), displayColors(cid,2), displayColors(cid,3), 0.8 );
                    volData.SetChannelRange( cid-1, displayRanges(cid,1), displayRanges(cid,2) );

                end        
            end

            this.imarisApp.mDataSet = volData;

            % Display the dataset in the newly created scene.
            this.imarisApp.mSurpassCamera.Fit;
            imarisScene.AddChild(this.imarisApp.mFactory.CreateVolume);                
            
            % store
            this.metadata.spacing = spacing;
            this.metadata.numTimePoints = numTimePoints;
            this.metadata.numChannels = numChannels;
            this.metadata.volSize = size( dataIn(:,:,:,1,1) );
            
        end
        
        function hContainer = AddDataContainer(this, hParent )
            
            if ~exist( 'hContainer', 'var' )
                hParent = this.imarisApp.mSurpassScene;
            end
            
            hContainer = this.imarisApp.mFactory.CreateDataContainer;
            hParent.AddChild(hContainer);
            
        end
        
        function AddSpots(this, spotLocations, spotTimes, varargin )
            
            numSpots = size(spotLocations, 1);
            
            p = inputParser;
            p.addRequired( 'spotLocations', @(x) (size(x,2) == 3) );
            p.addRequired( 'spotTimes', @(x) (isscalar(x) | (isvector(x) && numel(x) == numSpots)) );
            p.addParamValue( 'name', [], @ischar );
            p.addParamValue( 'color', random_color(), @(x) (isvector(x) && ismember(numel(x), [3,4])) );
            p.addParamValue( 'radii', 3 * min(this.metadata.spacing), @(x) (isscalar(x) | (isvector(x) && numel(x) == numSpots)) );
            p.addParamValue( 'trackEdges', [], @(x) ( (ischar(x) && stricmpi('consecutive')) && (ismatrix(x) && all(size(x) == [numSpots, 2]) && x(:) >= 1 && x(:) <= numSpots)) );
            p.addParamValue( 'hContainer', this.imarisApp.mSurpassScene, @(x) ishandle(x));
            p.parse( spotLocations, spotTimes, varargin{:} );
            
            PARAMETERS = p.Results;
            hContainer = PARAMETERS.hContainer;
            
            spotColor = PARAMETERS.color;
            if numel(spotColor) == 3
                spotColor(4) = 0;
            end
            
            spacing = this.metadata.spacing;
            volSize = this.metadata.volSize;
            
            % convert spot locations to physical space
            spotLocations(:,2) = volSize(1) - spotLocations(:,2) + 1; % flip-y
            spotLocationsPhysp = (spotLocations - 1) .* repmat(spacing([2,1,3]), [numSpots, 1] );        
            
            % get spot times
            if isscalar(spotTimes)
                spotTimes = repmat(spotTimes, numSpots,1);
            end
            
            % get spot radii
            if isscalar(PARAMETERS.radii)
               spotRadii = repmat(PARAMETERS.radii, numSpots, 1); 
            else
               spotRadii = PARAMETERS.radii; 
            end
            
            % create imaris spot object
            imarisSpots = this.imarisApp.mFactory.CreateSpots;
            imarisSpots.SetColor(spotColor(1), spotColor(2), spotColor(3), spotColor(4));
            imarisSpots.Set(spotLocationsPhysp, spotTimes, spotRadii);                
            
            % add track edges
            if ~isempty(PARAMETERS.trackEdges)
                
               if ischar(PARAMETERS.trackEdges)
                  trackEdges = ([0:numSpots-2; 1:numSpots-1])';  
               else
                  trackEdges = PARAMETERS.trackEdges - 1; % 0-based index?
               end
                
               imarisSpots.SetTrackEdges(trackEdges);
               
            end
            
            % name
            if ~isempty(PARAMETERS.name)
                imarisSpots.mName = PARAMETERS.name;
            end
            
            % add to container
            hContainer.AddChild(imarisSpots, -1);
            
        end
        
        function AddSurfaces(this, surfaceList, hContainer, varargin )
            
            if ~exist( 'hContainer', 'var' )
                hContainer = this.imarisApp.mSurpassScene;
            end
            
            p = inputParser;
            p.addRequired( 'surfaces', @(x) (isstruct(x) && sum(isfield(x, {'vertices', 'faces', 'normals', 'timepoint'})) == 4) );
            p.parse( surfaceList );
            
            numSurfaces = numel( surfaceList );
            
            p.addParamValue( 'name', [], @ischar );
            p.addParamValue( 'color', random_color(), @(x) (isvector(x) && ismember(numel(x), [3,4])) );
            p.addParamValue( 'tracks', [], @(x) ( (ischar(x) && stricmpi('consecutive')) || (ismatrix(x) && size(x,2) == 2 && all(x(:) >= 1) && all(x(:) <= numSurfaces))) );
            p.parse( surfaceList, varargin{:} );
            
            PARAMETERS = p.Results;
            
            spotColor = PARAMETERS.color;
            if numel(spotColor) == 3
                spotColor(4) = 0;
            end
            
            % create imaris surface object
            imarisSurfaceObject = this.imarisApp.mFactory.CreateSurfaces;
            
            % name
            if ~isempty(PARAMETERS.name)
                imarisSurfaceObject.mName = PARAMETERS.name;
            end
            
            % color
            imarisSurfaceObject.SetColor(spotColor(1), spotColor(2), spotColor(3), spotColor(4) );
            
            % add surfaces
            spacing = this.metadata.spacing;
            volSize = this.metadata.volSize;
            
            for sid = 1:numSurfaces

                curSurface = surfaceList(sid);
                
                % flip y-dim of vertices
                curSurface.vertices(:,2) = volSize(1) - curSurface.vertices(:,2) + 1;

                % convert vertex coordinates to physical space
                curSurface.vertices = (curSurface.vertices - 1) .* repmat(spacing([2,1,3]), [size(curSurface.vertices,1), 1] );        

                % convert face indices to 0-based index
                curSurface.faces = curSurface.faces - 1;
           
                % timepoint
                curSurface.timepoint =  curSurface.timepoint - 1; % 0-based?
                
                % add surface
                imarisSurfaceObject.AddSurface( curSurface.vertices, ...
                                                curSurface.faces, ...
                                                curSurface.normals, ...
                                                curSurface.timepoint );
                
            end
            
            % add track edges
            if ~isempty(PARAMETERS.tracks)
                
               if ischar(PARAMETERS.tracks)
                  trackEdges = ([0:numSurfaces-2; 1:numSurfaces-1])';  
               else
                  trackEdges = PARAMETERS.tracks - 1; % 0-based index?
               end
                
               imarisSurfaceObject.SetTrackEdges(trackEdges);
               
            end
            
            % add to container
            hContainer.AddChild(imarisSurfaceObject);
            
        end
        
    end
    
    methods (Static)
    
        function [surfaceObj] = generateSurfaceFromMask(imObjectMask, varargin)

            p = inputParser;
            p.CaseSensitive = false;
            p.addParamValue('timepoint', 1.0, @(x) isscalar(x) && isnumeric(x) && x > 0);
            p.addParamValue('surfaceQuality', 1.0, @(x) isscalar(x) && isnumeric(x) && x > 0.0 && x <= 1.0);
            p.parse( varargin{:} );
            
            PARAMETERS = p.Results;
            
            % crop the mask for better performance 
            pixcoord = ind2submat( size(imObjectMask), find(imObjectMask) );
            minCoord = min(pixcoord);
            maxCoord = max(pixcoord);
            cropind = cell(1, numel(minCoord));
            for dim = 1:numel(minCoord)
                cropind{dim} = minCoord(dim):maxCoord(dim);
            end
            
            imObjectMaskCropped = imObjectMask( cropind{:} );
            
            % generate isosurface
            imMaskSmoothed = smooth3( padarray(imObjectMaskCropped, ones(1,3), 0) );            
            surfaceObj = reducepatch( isosurface(imMaskSmoothed, 0.5), PARAMETERS.surfaceQuality );
            surfaceObj.normals = isonormals( imMaskSmoothed, surfaceObj.vertices );

            % correct vertex positions by adding offset
            surfaceObj.vertices(:,1) = cropind{2}(1) - 1 + surfaceObj.vertices(:,1);
            surfaceObj.vertices(:,2) = cropind{1}(1) - 1 + surfaceObj.vertices(:,2);
            surfaceObj.vertices(:,3) = cropind{3}(1) - 1 + surfaceObj.vertices(:,3);
            
            % set time point
            surfaceObj.timepoint = PARAMETERS.timepoint;
            
        end
        
    end
    
    methods (Access = private)
        
        function [vout] = typeCastVolume(this, vin, voxelType )

            switch voxelType

                case 'eTypeUInt8'            
                    vout = uint8( vin );

                case 'eTypeUInt16'

                    vout = uint16( vin );

                case 'eTypeFloat'

                    vout = single( vin );
            end

        end        
        
    end
    
    properties (SetAccess = private)
       
        metadata
        imarisApp
        
    end
        
end % classdef