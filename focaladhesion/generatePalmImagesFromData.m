
clc
clear 
close all

% ************************************************************************
%                         PARAMETERS                        
% ************************************************************************

    defaultDataDir = 'C:\deepak\data\fromLindsay\rawdata';    
    magnification = 1;

% ************************************************************************

rawdataRootDir = uigetdir( defaultDataDir, 'Select the directory containing the raw data files' );   
outputDir = uigetdir( defaultDataDir, 'Select the output directory' );   

dataFileList = dir( fullfile( rawdataRootDir, '*.txt' ) );
for i = 1:numel(dataFileList)
    
    fprintf( '\n\nProcessing data file %d/%d of size %d MB ... \n\n', i, numel(dataFileList), round(dataFileList(i).bytes/10^6) );
    
    [pathstr, name, ext] = fileparts( dataFileList(i).name );    
    
    tic
    [imPalm, imNMap, imZMap] = getPalmImage( fullfile(rawdataRootDir, dataFileList(i).name), 'magnification', magnification );
    timeElapsed = toc;
    
    fprintf( '\n\tTotal processing took %.2f seconds\n', timeElapsed );
    
    imwrite( imPalm, fullfile( outputDir, sprintf( '%s_palm_%dxMag.png', name, magnification ) ), 'png' );
    imwrite( mat2gray(ComputeImageLogTransform( imNMap )), fullfile( outputDir, sprintf( '%s_NMap_%dxMag.png', name, magnification ) ), 'png' );
    imwrite( mat2gray(ComputeImageLogTransform( imZMap )), fullfile( outputDir, sprintf( '%s_ZMap_%dxMag.png', name, magnification ) ), 'png' );
   
    save( fullfile( outputDir, sprintf( '%s_palm_%dxMag.mat', name, magnification ) ), 'imPalm', 'imNMap', 'imZMap' );
    
end


