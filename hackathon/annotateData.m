function annotateData( MovieListPath, datasetId )

    if ~exist( 'MovieListPath', 'var' )
        [fileName,pathName] = uigetfile( '*.mat', 'Select the Movie List File' );   
        MovieListPath = fullfile( pathName, fileName );
    end

    % load ML
    fprintf( '\n\nLoding Movie List ... \n\n' );
    
    load( MovieListPath ); 

    % run sanity check
    fprintf( '\n\nRunning Sanity Check ... \n\n' );

    [pathstr, name, ext] = fileparts( MovieListPath );
    ML.sanityCheck( pathstr );

    % get datasetId if not provided
    numDatasets = numel(ML.movieDataFile_);
    if ~exist( 'datasetId', 'var' )        
        
        strPrompt = sprintf( 'There are %d datsets\n', numDatasets );
        for i = 1:numDatasets        
            curMoviePath = ML.movieDataFile_{i};
            curDatasetName = curMoviePath( numel(ML.outputDirectory_)+2:end );
            strPrompt = sprintf( '%s\n%d/%d: %s', strPrompt, i, numDatasets, curDatasetName );
        end
        strPrompt = sprintf( '%s\n\nWhich one do you want to annotate?\n', strPrompt );
        datasetId = str2double( inputdlg( strPrompt, 'input dataset id', 1, {'1'} ) );
        
    end
    
    % load MD
    load( ML.movieDataFile_{datasetId} ); % creates MD

    % read image series
    imSeries = readImageSeries( MD.channels_.channelPath_, MD.channels_.getImageFileNames );

    % read rough mask
    maskFileNames = MD.processes_{1}.getOutMaskFileNames(1);
    maskFileNames = maskFileNames{1};
    imSeriesMask = readImageSeries( MD.processes_{1}.outFilePaths_{1}, maskFileNames );

    % ask user to nnotate
    [mask, isDone] = manualSegmentationTweakGUI( imSeries, imSeriesMask );

    % save annotation
    resultsDir = fullfile( MD.packages_{1}.outputDirectory_, 'annotated_masks', 'annotated_masks_for_channel_1' );
    if ~isdir(resultsDir)
       mkdir(resultsDir); 
    end
    prefix = sprintf( 'annotated_mask_%s', MD.outputDirectory_( numel(ML.outputDirectory_)+2:end ) );
    writeImageSeries( mask, resultsDir, prefix );
    save( fullfile(resultsDir, 'annotation.mat'), 'mask', 'isDone' );

end
