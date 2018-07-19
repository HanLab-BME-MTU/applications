
clc
clear 
close all

% ************************************************************************
%                         PARAMETERS                        
% ************************************************************************

    defaultDataDir = 'C:\deepak\data\fromLindsay';    
    adhesionAreaThreshold = 150;    

% ************************************************************************

dataDir = uigetdir( defaultDataDir, 'Select Folder Contained Palm Data Files' );

outputDir = uigetdir( defaultDataDir, 'Select Results Folder' );

dataFileList = dir( fullfile( dataDir, '*.mat' ) );

for fid = 1:numel(dataFileList)
    
    fprintf( '\n\n%d/%d: Processing file %s ...\n\n', fid, numel(dataFileList), dataFileList(fid).name );
    
    [pathstr, curFileName, ext] = fileparts(dataFileList(fid).name); 
    dataFilePath = fullfile( dataDir, dataFileList(fid).name );    
    dataFileContents = load( dataFilePath );
    imInput = double( dataFileContents.imNMap );
    ptInd = find( imInput > 0 );
    [yind, xind] = ind2sub( size(imInput), ptInd );

    [ imAdhesionSegLabel, adhesionStats ] = segAdhesionsClustering( imInput, ...
                                                                    'minAdhesionArea', adhesionAreaThreshold, ...
                                                                    'minClusterDistance', 40, ...
                                                                    'sparsityFactor', 0.5, ...
                                                                    'flagDebugMode', false);

    [imSegRGBMask, labelColorMap] = label2rgbND( imAdhesionSegLabel );    
    imRGBMaskOverlay = genImageRGBMaskOverlay( mat2gray(ComputeImageLogTransform( imInput )), imSegRGBMask, 0.2 );    
    imwrite( imRGBMaskOverlay, fullfile( outputDir, [curFileName, '_SegResult.png'] ), 'png' );
    
    imshow( imRGBMaskOverlay );
    set(gca, 'Units', 'normalized');
    set(gca, 'Position', [0,0,1,1]);
    hold on;
    for i = 1:max(imAdhesionSegLabel(:))
       ptCenter = adhesionStats(i).Centroid;
       a = 0.5 * adhesionStats(i).MajorAxisLength;
       b = 0.5 * adhesionStats(i).MinorAxisLength;

       theta = adhesionStats(i).Orientation;
       ptRect = ([ -a, -b; ...
                    a, -b; ...
                    a, b; ...
                   -a, b; ...
                   -a, -b ])';
       ptRect(end+1,:) = 1;
       transMat = [ cosd(theta), sind(theta), ptCenter(1); ...
                   -sind(theta) cosd(theta), ptCenter(2); ...
                    0, 0, 1 ];
       ptRect = round( transMat * ptRect );
       plot( ptRect(1,:), ptRect(2,:), '-', 'Color', labelColorMap(i,:), 'LineWidth', 2.0 ); 

       plot( ptCenter(1), ptCenter(2), 'o', 'Color', [0,0,0], 'MarkerFaceColor', labelColorMap(i,:), 'LineWidth', 2.0 ); 
    end
    hold off;
    SaveFigure( gcf, fullfile( outputDir, [curFileName, '_SegResultBBox.png'] ) );
    
end