function [ imSeries ] = readImageSeries( rootdir, filenames )
%%  reads a series of time lapse images
% 
%   [ imSeries ] = readImageSeries( rootdir, filenames )
% 
%   Input Arguments:
% 
%       rootdir: directory in which the images are located
%       filenames: a cell array of filenames of the image series in the
%                  correct order
% 
%   Output Arguments:       
%       
%       imSeries: A z-stack/3d-matrix of images stacked on top of each
%                 other
% 
%   Author: Deepak Roy Chittajallu
% 

    fprintf( '\n\nReading image series in folder %s ... \n\n', rootdir );

    imSeries = [];
    
    for i = 1:numel(filenames)
       
        fprintf( '\n%.3d/%.3d: reading image %s ...', i, numel(filenames), filenames{i} );
        
        imSeries = cat(3, imSeries, imread(fullfile(rootdir, filenames{i})) );
        
    end

    fprintf( '\n\n' );
    
end