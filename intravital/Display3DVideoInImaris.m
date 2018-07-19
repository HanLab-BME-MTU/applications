function [ varargout ] = Display3DVideoInImaris( videoInput, varargin )

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
    p.parse( videoInput, varargin{:} );
    
    PARAMETERS = p.Results;
    spacing = PARAMETERS.spacing;
    displayColors = PARAMETERS.displayColors;
    displayRanges = PARAMETERS.displayRanges;
    units = PARAMETERS.units;
    voxelType = PARAMETERS.voxelType;
    
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