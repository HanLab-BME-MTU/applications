clc
clear
close all

% if ~exist( 'MovieListPath', 'var' )
    [fileName,pathName] = uigetfile( fullfile('C:\deepak\data\hackathon', '*.mat'), 'Select the Movie List File' );   
    MovieListPath = fullfile( pathName, fileName );
% end

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
imSeries = readImageSeries( MD.channels_(1).channelPath_, MD.channels_.getImageFileNames );

% read rough mask
maskFileNames = MD.processes_{1}.getOutMaskFileNames(1);
maskFileNames = maskFileNames{1};
imSeriesMask = readImageSeries( MD.processes_{1}.outFilePaths_{1}, maskFileNames );

% visualize
imseriesmaskshow( imSeries, imSeriesMask );

% % topo graph-cuts interactive
%imLabeledMask = segmentCellsInteractivelyUsingTopoGraphcut( imSeries );

% topo graph-cuts automatic
imLabeledMask = segmentCellsUsingTopoGraphcut( imSeries );

% % segment interactively using graph cuts
% switch ndims(imSeries)
% 
%     case 2
%     
%         imSegMask = segObject2D_InteractiveGraphCuts( imSeries );    
%         
%     case 3
%     
%         imSegMask = segObject3D_InteractiveGraphCuts( imSeries );
%         
% end

% % topo chan vese
% imSeries = filterGaussND( imSeries, 1.5 );
% 
% [fg, bg] = get_fgnd_bgnd_seeds_3d_points( imSeries );
% imSeries = imSeries(:,:,fg(1,3));
% imSeedMask = zeros( size(imSeries) );
% imSeedMask( fg(1,2), fg(1,1) ) = 1;
% imSeedMask( fg(2,2), fg(2,1) ) = 1;
% imSeedMask = imdilate( imSeedMask, strel('disk',10) );
% 
% seg1 = gc_chan_vese(imSeries, 1, 10, 10);
% seg2 = gc_chan_vese_tp(imSeries, imSeedMask > 0, 1, 10, 10);
% imseriesmaskshow( imSeries, {seg1, seg2, imSeedMask} );