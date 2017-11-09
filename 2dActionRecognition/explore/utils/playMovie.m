function playMovie(movie)

    figure();
    for i=1:size(movie,3)
        imshow(movie(:,:,i));
        input('press enter');
    end

end