function writeImageSeries( imSeries, outdir, prefix, fileformat )
%%  write a series of time lapse images to disk
% 
%   [ imSeries ] = writeImageSeries( outdir, filenames )
% 
%   Input Arguments:
% 
%       outdir: directory in which the images should be written
%       filenames: prefix to add before filename
% 
%   Author: Deepak Roy Chittajallu
% 

    if ~exist( 'fileformat', 'var' )
        fileformat = 'tif';
    end
    
    fprintf( '\n\nWriting image series to folder %s ... \n\n', outdir );

    numFrames = size(imSeries,3);
    
    for i = 1:numFrames
       
        curOutFileName = sprintf( '%s_t%.3d.%s', prefix, i, fileformat );
        fprintf( '\n%.3d/%.3d: writing image %s ...', i, numFrames, curOutFileName );        
        imwrite( imSeries(:,:,i), fullfile(outdir, curOutFileName), fileformat );
        
    end

    fprintf( '\n\n' );
    
end