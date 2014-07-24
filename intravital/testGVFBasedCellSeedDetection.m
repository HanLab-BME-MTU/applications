
    clc
    clear all
    close all

    stefanOldDataDir = 'C:\deepak\data\Stefan_June2012';
    stefanNewDataDir = 'C:\deepak\data\Stefan_July2012';
    ralphDataDir = 'C:\deepak\data\Ralph_June2012';
    adhesionDir = 'C:\deepak\data\adhesions';

    dataFileDir = stefanNewDataDir;

    if ~exist( 'dataFilePath', 'var' )
        [fileName,pathName] = uigetfile( fullfile( dataFileDir, '*.mat' ), 'Select the data file' );   
        dataFilePath = fullfile( pathName, fileName )
    end

    dataFileContents = load( dataFilePath );
    imInput = dataFileContents.im;

    % adjust the intensity range of the image
    minObjectDiameter = [ 10, 10 ];
    ImageIntensityRange = ComputeImageDynamicRange( imInput, 99.0 );
    imAdjusted = AdjustImageIntensityRange( imInput, ImageIntensityRange );    

    % apply median filter to remove any spotty noise
    imAdjusted = medfilt2( imAdjusted, [3,3] );    

    % compute threshold
    imThresh = zeros( size( imAdjusted ) );
    numHistogramBins = 256; 
    thresholdingAlgorithm = 'FOtsuThreshold';
    imLOG = mat2gray( ComputeImageLogTransform( imAdjusted ) ); 
    imThresh = callMatitkFilter( thresholdingAlgorithm, { numHistogramBins } , imLOG );                        
    
    % detect cell seed points    
%     imCellSeedPoints = detect_cell_seeds_gvf_edge_map( imAdjusted, minObjectDiameter, ...
%                                                        'debugMode', true, ...
%                                                        'sigma', 0.5 );

    imCellSeedPoints = detect_cell_seeds_gvf_blurred_image( imerode( imAdjusted, strel('ball',3,3) ), ...
                                                            minObjectDiameter, ...
                                                            'debugMode', true );
    imCellSeedPoints( ~imThresh ) = 0;
    
    imseriesmaskshow( imAdjusted, imThresh );
    hold on;
    [yind,xind] = ind2sub( size(imAdjusted), find( imCellSeedPoints > 0 ) );
    plot( xind, yind, 'g+', 'MarkerSize', 10.0 );
    hold off;
    set( gcf, 'Name', 'Cell Seed Point Detection Result' );

    


