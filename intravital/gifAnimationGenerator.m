

datadir = uigetdir( 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\dsp\goodones' );
dircontents = dir( fullfile( datadir, '*.tif' ) );

fprintf( '\nGenerating Gif file ...\n' );
figure;
for i = 1:numel(dircontents)
    fprintf( '\nAdding file %d/%d\n', i, numel(dircontents) );
    im = imread( fullfile( datadir, dircontents(i).name ) );
    imshow( rgb2gray(im) );    
    plotTimeStamp( sprintf( '%d/%d', i, numel(dircontents) ), 'Location', 'SouthWest', 'Color', [1 1 1] );
    im = getframe(gca);
    im = rgb2gray(im.cdata);
    if i == 1
        imwrite( im, fullfile(datadir, 'gif_animation.gif'), 'gif', 'Loopcount', inf );
    else
        imwrite( im, fullfile(datadir, 'gif_animation.gif'), 'gif', 'WriteMode', 'append' );
    end
end