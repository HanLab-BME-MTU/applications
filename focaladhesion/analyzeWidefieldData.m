clc
clear 
close all

% ************************************************************************
%                         PARAMETERS                        
% ************************************************************************

    defaultDataDir = 'C:\deepak\data\fromLindsay\widefieldlike';    

% ************************************************************************

if ~exist( 'dataFilePath', 'var' )
    [fileName,pathName] = uigetfile( fullfile( defaultDataDir, '*.tif' ), 'Select the focal adhesion widefield data' );   
    dataFilePath = fullfile( pathName, fileName )
end

im = double( imread( dataFilePath ) );
imseriesshow( im );
figure, histogram( im( im > 0 ) );

imseriesmaskshow( im, im > 19 );