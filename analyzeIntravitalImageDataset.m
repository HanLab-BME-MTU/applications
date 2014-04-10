%function analyzeIntravitalImageDataset( dataFilePath )

    clc
    clear 
    %close all
    
    stefanOldDataDir = 'C:\deepak\data\Stefan_June2012';
    stefanNewDataDir = 'C:\deepak\data\Stefan_July2012';
    stefanDilutedDataDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\goodones';
    ralphDataDir = 'C:\deepak\data\Ralph_June2012';
    adhesionDir = 'C:\deepak\data\adhesions';
    intravitalDataDir = 'Z:\intravital\data';
    milesDataDir = 'Z:\intravital\Miles-Data\0689_M08_PKPD_HT1080_53bp1truncappple_1_9_2013';
    
    dataFileDir = intravitalDataDir; 
    % dataFileDir = stefanDataDir; 
    % dataFileDir = stefanNewDataDir; 
    
    if ~exist( 'dataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.oif; *.oib; *.tif' ), 'Select the data file' );   
        dataFilePath = fullfile( pathName, fileName )
    end

    % load data using bfopen from the Bioformats toolbox
    data_bfopen = ( bfopen( dataFilePath ) )';    
    metadata_bfopen = data_bfopen{4};
    
    metadata.numChannels = metadata_bfopen.getPixelsSizeC(0).getValue;
    metadata.numTimePoints = metadata_bfopen.getPixelsSizeT(0).getValue;
    metadata.volSize = [ metadata_bfopen.getPixelsSizeX(0).getValue, metadata_bfopen.getPixelsSizeY(0).getValue metadata_bfopen.getPixelsSizeZ(0).getValue ];
    metadata.voxelSpacing = [ metadata_bfopen.getPixelsPhysicalSizeX(0).getValue metadata_bfopen.getPixelsPhysicalSizeY(0).getValue metadata_bfopen.getPixelsPhysicalSizeZ(0).getValue ];    
    
%     if (metadata.voxelSpacing(1)/ metadata.voxelSpacing(3)) >= 1
%         metadata.voxelSpacing = [0.5, 0.5, 2]; % incorrect spacing in metadata - use something meaningful
%     end
        
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
        imageDataMIP{c} = zeros( [metadata.volSize(1:2), metadata.numTimePoints ] );
        for t = 1:metadata.numTimePoints
            imageDataMIP{c}(:,:,t) = max( imageData{t,c}, [], 3 );            
        end
    end         
      
    %% segment cells using watershed with seed points detected as local intensity maxima   
    %[ imLabelCellSeg ] = segmentCellsInIntravitalData( imageData{1,2}, metadata.voxelSpacing, 'flagDisplayResultsInImaris', true, 'flagParallelize', true);
%     [L, imSegRGBMask] = experiment_watershed_segmentation( imageData{1,1}, metadata, true );
%     segMovie = MakeMovieFromRGBMask( imageData{1,1}, imSegRGBMask );
%     [pathstr, name, ext] = fileparts( dataFilePath );
%     movie2avi( segMovie, fullfile( 'C:\deepak\results\WatershedSeg', [name '.avi'] ), 'FPS', 5 );

    %% segment cells using watershed with seed points detected as sinks in diffused GVF
    % experiment_gvf_based_seed_detection( imageData{1,1}, metadata, true );

    %% segment cells using watershed with seed points detected as maxima of Adaptive LoG
    % experiment_watershed_cell_seg_autoscale( imageData{1,1}, metadata, true );

    %% segment cell using gradient weighted watershed 
    % experiment_gradient_weighted_watershed( imageData{1,1}, metadata, true );    
    
    %% compare various global thresholding methods
    %experiment_global_thresholding( imageData{1,1}, metadata );
    
    %% apply global thresholding algorithms in each slice separately
    %experiment_global_thresholding_per_slice( imageData{1,1}, metadata );
    
    %% compare various local thresholding methods   
    %experiment_local_thresholding( imageData{1,1}, metadata );    
    
    %% apply locall adaptive thresholding in each slice
    %experiment_local_otsu_thresholding( imageData{1,1}, metadata, true );
    
    %% try hough-style voting in gradient-vector direction to enhance the 
    % core of each nucleus and try getting a seed point in each nucleus
    
    %% compute minimum error thresholds for each slice and see how they vary    
    % experiment_check_slice_threshold_variation( imageData, imageDataMIP, metadata );
    
    %% check if slices are getting darker as we go deeper in z
    % experiment_check_slice_mean_intensity_variation( imageData, metadata );
    
     %% randomly shoot 10 rays through the 3D volume in the time axis and
     %  visualize them to see if background subtraction can help
     %  experiment_will_background_subtraction_help( imageData, metadata );
    
%end
