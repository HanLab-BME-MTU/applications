function ML = movieData2movieList(movieObj)
%If you can't figure it out what it does, you shouldn't use it
    ML = MovieList(movieObj,movieObj.movieDataPath_);
    ML.setPath(movieObj.movieDataPath_)
    ML.setFilename('movieList.mat')
    ML.sanityCheck;

end

