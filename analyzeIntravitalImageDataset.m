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
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, bfGetFileExtensions), 'Select the data file' );   
        dataFilePath = fullfile( pathName, fileName )
    end

    % load data using bfopen from the Bioformats toolbox
    imageSeries = loadIntravitalDataset( dataFilePath );
    
    metadata = imageSeries(1).metadata;
    imageData = imageSeries(1).imageData;
    
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
    chid = 3;    
    %experiment_global_thresholding( double(imageData{chid}), metadata );    

    imPreprocessed = matitk( 'FMEDIAN', 2 * [1, 1, 0], double(imageData{chid}) );
    %imPreprocessed = filterGaussND(double(imageData{chid}), min(metadata.voxelSpacing), 'spacing', metadata.voxelSpacing);

%     [ imThresholdSurface, imMask ] = thresholdVariationalMinMaxOpt( imPreprocessed );
%     imseriesmaskshow(imageData{chid}, imMask);
    
    krnlMax = streldisknd( round(0.5 * 40 ./ metadata.voxelSpacing(1:2)) );
    imLocalBackground = imopen(imPreprocessed, krnlMax);
    imSignalToBackgroundRatio = imPreprocessed ./ (eps + imLocalBackground);
    imseriesmaskshow(imageData{chid}, imSignalToBackgroundRatio > 2.0);    
    
    %% apply global thresholding algorithms in each slice separately
    %experiment_global_thresholding_per_slice( double(imageData{1,1}), metadata );
    
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
