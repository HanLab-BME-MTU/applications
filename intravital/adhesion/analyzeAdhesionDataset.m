%function analyzeIntravitalImageDataset( dataFilePath )

    clc
    clear 
    close all
    
    dataFileDir = 'C:\deepak\data\adhesions';    
    
    if ~exist( 'dataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.oif' ), 'Select the data file' );   
        dataFilePath = fullfile( pathName, fileName )
    end

    % load data using bfopen from the Bioformats toolbox
    data_bfopen = ( bfopen( dataFilePath ) )';    
    metadata_bfopen = data_bfopen{4};
    
    metadata.numChannels = metadata_bfopen.getPixelsSizeC(0).getValue;
    metadata.numTimePoints = metadata_bfopen.getPixelsSizeT(0).getValue;
    metadata.volSize = [ metadata_bfopen.getPixelsSizeX(0).getValue, metadata_bfopen.getPixelsSizeY(0).getValue metadata_bfopen.getPixelsSizeZ(0).getValue ];
    metadata.voxelSpacing = [ metadata_bfopen.getPixelsPhysicalSizeX(0).getValue metadata_bfopen.getPixelsPhysicalSizeY(0).getValue metadata_bfopen.getPixelsPhysicalSizeZ(0).getValue ];    
  
    metadata
    
    imageData = cell( metadata.numTimePoints, metadata.numChannels );  
    for t = 1:metadata.numTimePoints         
        for c = 1:metadata.numChannels 
            
            indx = (t-1) * metadata.volSize(3) * metadata.numChannels + [ c : metadata.numChannels : metadata.volSize(3) * metadata.numChannels ];
            imageData{t,c} = double( cat(3, data_bfopen{1}{ indx , 1 }) );
            
        end
    end   
        
    % generate MIP video stacks for each channel    
    imageDataMIP = cell( 1, metadata.numChannels );    
    for c = 1:metadata.numChannels         
        volSize = size( imageData{1,c} );
        imageDataMIP{c} = zeros( [ volSize(1:2), metadata.numTimePoints ] );
        for t = 1:metadata.numTimePoints
            imageDataMIP{c}(:,:,t) = max( imageData{t,c}, [], 3 );            
        end
    end         
      
    %% compare various global thresholding methods
    %experiment_global_thresholding( imageData{1,1}, metadata );
    
    %% apply global thresholding algorithms in each slice separately
    % experiment_global_thresholding_per_slice( imageData{1,2}, metadata );
    
    %% compare various local thresholding methods   
    % experiment_local_thresholding( imageData{1,1}, metadata );    
    
    % segment cell and detect spots
    % experiment_segment_adhesion_cell( imageData{1,2}, metadata, true );
    
    % segment collagen fibers
    experiment_segment_adhesion_fibers( imageData{1,1}, metadata ); 
    
%end
