movie = avifile ('background.avi');

for movieStep = 1 : 205
    
    fprintf (1, 'Frame: %d ', movieStep);
    
    filename = ['/lccb/projects/polytrack/bg_movie/bgImage_' num2str(movieStep) '.mat'];
    
    load (filename);
    
    fig = figure ('Position', [200 300 500 400]), mesh (bgImage);
    axis ([0 1500 0 1500 0.35 0.46])
    
    F = getframe (fig);
      
    movie = addframe (movie, F);
      
    close;
    
end

movie = close(movie);