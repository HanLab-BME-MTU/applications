    clc
    clear all
    close all   
    
    stefanOldDataDir = 'C:\deepak\data\Stefan_June2012';
    stefanNewDataDir = 'C:\deepak\data\Stefan_July2012';
    ralphDataDir = 'C:\deepak\data\Ralph_June2012';
    adhesionDir = 'C:\deepak\data\adhesions';
    stefanDilutedDataDir = 'Z:\intravital\data\Stefan_September2012\goodones';
    
    % dataFileDir = ralphDataDir; 
    %dataFileDir = stefanDataDir; 
    dataFileDir = stefanDilutedDataDir; 
    
    if ~exist( 'twoPhotonDataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.oif' ), 'Select the two-photon data file' );   
        twoPhotonDataFilePath = fullfile( pathName, fileName )
    end

    if ~exist( 'confocalDataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.oif' ), 'Select the confocal data file' );   
        confocalDataFilePath = fullfile( pathName, fileName )
    end
    
    % load two photon data using bfopen from the Bioformats toolbox
    [ twoPhotonImageSeries ] = loadIntravitalDataset( twoPhotonDataFilePath );    
    [ confocalImageSeries ] = loadIntravitalDataset( confocalDataFilePath );
    
    imseriesshow_multichannel( {twoPhotonImageSeries.imageData{1,1}, confocalImageSeries.imageData{1,1}, confocalImageSeries.imageData{1,2}} );
    set( gcf, 'Name', 'Pre-Shift: Two Photon and Confocal Channel Overlay' );
        
    % perform rigid registration to correct shift    
    [mfilepath, name, ext] = fileparts( mfilename('fullpath') );    
    elastixParameters{1} = readElastixParameters( fullfile( mfilepath, 'elastixParamFile_rigid3D.txt' ) );

    mergedImageData{1} = twoPhotonImageSeries(1).imageData{1,1};
    
    for i = 1:2
        
        % setup fixed and moving image data
        dataFixed.im = confocalImageSeries(1).imageData{1,i};
        dataFixed.Spacing = confocalImageSeries(1).metadata.voxelSpacing;   
        
        dataMoving.im = twoPhotonImageSeries(1).imageData{1,1};
        dataMoving.Spacing = twoPhotonImageSeries(1).metadata.voxelSpacing;

        imseriesshow_multichannel( {dataMoving.im, dataFixed.im} );
        set( gcf, 'Name', sprintf('Overlay of Moving and Fixed Volume Before Registration - channel %d', i) );
        
        % compute fixed mask
        ImageIntensityRange = ComputeImageDynamicRange( dataFixed.im, 99.0 );
        imAdjusted = AdjustImageIntensityRange( dataFixed.im, ImageIntensityRange );
        imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );
        
        fixedMask = callMatitkFilter( 'FOtsuThreshold', { 256 } , ComputeImageLogTransform( imAdjusted ) );
        fixedMask = imdilate(fixedMask, ones(3,3,3));

        imseriesmaskshow( dataFixed.im, fixedMask );
        set( gcf, 'Name', sprintf( 'Fixed Mask - channel %d', i ) );
        
        % perform registration
        elastixParameters{1}.DefaultPixelValue = -1;
        [dataMovingR, elastixTransformParameters, itersInfo] = elastixRegister(dataFixed, dataMoving, elastixParameters, 'fixedMask', fixedMask);                         

        elastixParameters{1}.DefaultPixelValue = -1;
        elastixTransformParameters{1}.TransformParameters = -1 * elastixTransformParameters{1}.TransformParameters;
        dataFixedR = elastixTransform( dataFixed.im, dataFixed.Spacing, elastixTransformParameters );
        
        dataFixedR( dataFixedR < 0 ) = min(dataFixed.im(:));
        mergedImageData{i+1} = dataFixedR;
        
        imseriesshow_multichannel( {dataMoving.im, dataFixedR} );
        set( gcf, 'Name', sprintf( 'Overlay of Moving and Fixed Volume After Registration - channel %d', i ) );        
        
    end
    
    imseriesshow_multichannel( mergedImageData, ...
                               'spacing', twoPhotonImageSeries(1).metadata.voxelSpacing, ...
                               'displayColors', [0 0 1; 0 1 0; 1 0 0] );
    set( gcf, 'Name', 'Post-Shift: Two Photon and Confocal Channel Overlay' );
    
    colorWeights = ones(1,3)/3;
    mergedGrayImage = zeros( size(mergedImageData{1}) );
    for i = 1:3
        mergedGrayImage = mergedGrayImage + colorWeights(i) * mat2gray(mergedImageData{i});        
    end
    imseriesshow( mergedGrayImage );
    
    
    
% % perform rigid registration to correct shift    
% [mfilepath, name, ext] = fileparts( mfilename('fullpath') );    
% elastixParameters{1} = readElastixParameters( fullfile( mfilepath, 'elastixParamFile_rigid3D.txt' ) );
% 
%     % setup fixed and moving image data
%     dataFixed.im = 0.5 * (confocalImageSeries(1).imageData{1,1} + confocalImageSeries(1).imageData{1,2});
%     dataFixed.Spacing = confocalImageSeries(1).metadata.voxelSpacing;   
% 
%     dataMoving.im = twoPhotonImageSeries(1).imageData{1,1};
%     dataMoving.Spacing = twoPhotonImageSeries(1).metadata.voxelSpacing;
% 
%     imseriesshow_multichannel( {dataMoving.im, dataFixed.im} );
%     set( gcf, 'Name', 'Overlay of Moving and Fixed Volume Before Registration' );
% 
%     % compute fixed mask
%     ImageIntensityRange = ComputeImageDynamicRange( dataFixed.im, 99.0 );
%     imAdjusted = AdjustImageIntensityRange( dataFixed.im, ImageIntensityRange );
%     imAdjusted = matitk( 'FMEDIAN', [1,1,1], imAdjusted );
% 
%     fixedMask = callMatitkFilter( 'FOtsuThreshold', { 256 } , ComputeImageLogTransform( imAdjusted ) );
% 
%     for i = 1:size(fixedMask,3)
%         fixedMask(:,:,i) = imdilate( fixedMask(:,:,i), strel( 'disk', 3 ) );            
%     end
% 
%     imseriesmaskshow( dataFixed.im, fixedMask );
%     set( gcf, 'Name', 'Fixed Mask' );
% 
%     % perform registration
%     elastixParameters{1}.DefaultPixelValue = -1;
%     [dataMovingR, elastixTransformParameters, itersInfo] = elastixRegister(dataFixed, dataMoving, elastixParameters, 'fixedMask', fixedMask);                         
% 
%     imseriesshow_multichannel( {double(dataMovingR.im), dataFixed.im} );
%     set( gcf, 'Name', 'Overlay of Moving and Fixed Volume After Registration' );
% 
%     mergedImageData{1} = twoPhotonImageSeries(1).imageData{1,1};
%     elastixTransformParameters{1}.TransformParameters = -1 * elastixTransformParameters{1}.TransformParameters;
%     for i = 1:2
%         ImageIntensityRange = ComputeImageDynamicRange( confocalImageSeries(1).imageData{1,i}, 99.0 );
%         elastixTransformParameters{1}.DefaultPixelValue = -1;
%         mergedImageData{i+1} = elastixTransform( confocalImageSeries(1).imageData{1,i}, confocalImageSeries(1).metadata.voxelSpacing, elastixTransformParameters );
%     end

