function [ ML ] = toMovieList( in )
%lamins.util.toMovieList Utility function to convert different kinds of input into a
%MovieList

if(isa(in,'MovieList'))
    ML = in;
    if(isempty(ML.movies_))
        ML.sanityCheck;
    end
    return;
end

if(ischar(in))
    try
        ML = MovieList.load(in);
        return;
    catch err
        % Maybe not a MovieList
        if(err.identifier('lccb:movieObject:load'))
            in = MovieData.load(in);
        end      
    end
end

if(isa(in,'MovieData'))
    ML = MovieList(MD,pwd);
    return;
end


end

