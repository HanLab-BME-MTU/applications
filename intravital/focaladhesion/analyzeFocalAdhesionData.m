
clc
clear 
close all

% ************************************************************************
%                         PARAMETERS                        
% ************************************************************************

    defaultDataDir = 'C:\deepak\data\fromLindsay';    
    adhesionAreaThreshold = 150;    

% ************************************************************************

if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile( fullfile( defaultDataDir, '*.mat' ), 'Select the data file' );   
    dataFilePath = fullfile( pathName, fileName )
end

dataFileContents = load( dataFilePath );
imInput = double( dataFileContents.imNMap );
ptInd = find( imInput > 0 );
[yind, xind] = ind2sub( size(imInput), ptInd );

imseriesshow( imInput );
set( gcf, 'Name', 'Input PALM Rendered Image' );

imZMap = dataFileContents.imZMap;
figure, plot3( xind, yind, imZMap(ptInd), 'ro' );

%[ imAdhesionSegLabel, adhesionStats ] = segAdhesionsBlobColoring( imInput, adhesionAreaThreshold );
[ imAdhesionSegLabel, adhesionStats ] = segAdhesionsClustering( imInput, ...
                                                                'minAdhesionArea', adhesionAreaThreshold, ...
                                                                'minClusterDistance', 30, ...
                                                                'flagShowBBox', true, ...
                                                                'sparsityFactor', 4.0);