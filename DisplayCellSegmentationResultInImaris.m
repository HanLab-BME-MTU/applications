function [varargout] = DisplayCellSegmentationResultInImaris(imageData, imLabelCellSeg, CellSegColorMap, imCellSeedPoints, varargin)

    % default color map -- only 6 channels are defined
    % if more than 6 channels then user should specify the displayColors 
    cMap(1,:) = [1 0 0];
    cMap(2,:) = [0 1 0];
    cMap(3,:) = [0 0 1];
    cMap(4,:) = [1 0 1];
    cMap(5,:) = [0 1 1];
    cMap(6,:) = [1 1 0];

    if isnumeric(imageData) 
        imageData = {imageData};
    end
    
    % default display range
    numChannels = numel(imageData);
    default_displayranges = zeros(numChannels,2);
    volSize = size(imageData{1});
    for i = 1:numChannels
       curChannelVolSize = size(imageData{i});
       curChannelVolSize(3) = size(imageData{i}, 3);
       if ~isempty( setxor( curChannelVolSize, volSize ) )
          error('Images corresponding to all channels must of the same size and dimension'); 
       end
       default_displayranges(i,:) = ComputeImageDynamicRange(imageData{i}, 98.0);
    end
    
    % parse arguments
    p = inputParser;
    p.addParamValue( 'spacing', ones(1,3), @(x) (isnumeric(x) && numel(x) == 3) );
    p.addParamValue( 'displayColors', cMap(1:numChannels,:), @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 3 && size(x,1) == numChannels ) );
    p.addParamValue( 'displayRanges', default_displayranges, @(x) ( isnumeric(x) && ndims(x) == 2 && size(x,2) == 2 && size(x,1) == numChannels ) );
    p.addParamValue( 'spotRadius', 3, @(x) (isnumeric(x) && isscalar(x)) );
    p.parse(varargin{:});
    
    PARAMETERS = p.Results;
    
    % compute isosurface geometry for each cell
    hStatusDlg = waitbar( 0, 'computing surface geometry for each cell' );
    cellSurfaceObjectList = {};
    
    cellStats = regionprops(imLabelCellSeg, {'PixelIdxList', 'BoundingBox'});
    
    for cid = 1:numel(cellStats)        
        
        waitbar( cid/numel( cellStats ), hStatusDlg );

        if numel(cellStats(cid).PixelIdxList) == 0 || cellStats(cid).BoundingBox(6) <= 1
            continue;
        end
        
        % name
        curCellIsoSurface.name = sprintf( 'CellSeg_%d', cid);
        
        % color
        curCellIsoSurface.color = CellSegColorMap(cid, :);
        
        % create crop indices
        curCellPixelCoord = ind2submat( size(imLabelCellSeg), cellStats(cid).PixelIdxList );
       
        subinds = cell(1,3);
        for i = 1:3
            subinds{i} = min(curCellPixelCoord(:,i)):max(curCellPixelCoord(:,i));
        end    
        
        % crop segmentation mask
        imCurCellSegCropped = padarray( double(imLabelCellSeg(subinds{:}) == cid), ones(1,3), 0 );                                        
        imCurCellSegSmoothed = smooth3( imCurCellSegCropped );            
        curCellSurfaceGeometry = isosurface( imCurCellSegSmoothed, 0.5 );
        curCellSurfaceGeometry.normals = isonormals( imCurCellSegSmoothed, curCellSurfaceGeometry.vertices );
        
        % correct vertex positions by adding offset
        curCellSurfaceGeometry.vertices(:,1) = subinds{2}(1) + curCellSurfaceGeometry.vertices(:,1) - 1;
        curCellSurfaceGeometry.vertices(:,2) = subinds{1}(1) + curCellSurfaceGeometry.vertices(:,2) - 1;
        curCellSurfaceGeometry.vertices(:,3) = subinds{3}(1) + curCellSurfaceGeometry.vertices(:,3) - 1;
        
        % add cell surface to cell pattern surface object
        curCellIsoSurface.surfaces = curCellSurfaceGeometry;
        cellSurfaceObjectList{end+1} = curCellIsoSurface;
        
    end
    closeStatusDialog( hStatusDlg );
    
    % get cell seed point locations
    seedStats = regionprops( bwlabeln( imCellSeedPoints ), 'Centroid' );
    cellSeedPointLocations = cat( 1, seedStats.Centroid );

    % Display everything in imaris
    imarisApp = DisplayMultichannel3DDataInImaris( imageData, ...
                                                   'surfaceObjects', cellSurfaceObjectList, ...
                                                   'spotLocations', cellSeedPointLocations, ...
                                                   'spotRadius', PARAMETERS.spotRadius, ...
                                                   'spacing', PARAMETERS.spacing, ...
                                                   'displayRanges', PARAMETERS.displayRanges, ...
                                                   'displayColors', PARAMETERS.displayColors );

    % return imaris app object if requested
    if nargout >= 1
        varargout{1} = imarisApp;
    end
                                                       
end