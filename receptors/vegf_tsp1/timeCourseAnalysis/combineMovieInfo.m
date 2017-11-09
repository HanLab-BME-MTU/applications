function [ combinedMovieInfo, movieInfoArray ] = combineMovieInfo( ML , movieIndices)
%combineMovieInfo Combine several movieInfo structures together
%
% ML is a MovieList object containing the movies of interest
% movieIndices are the movies numbers within ML that are of interest
%
% mergedMovieInfo 
%
% Mark Kittisospikul for Sungsoo Lee

if(nargin < 2)
    movieIndices = 1:length(ML.movies);
end

movies = ML.movies_(movieIndices);
outFiles = cellfun(@(MD) MD.processes_{1}.outFilePaths_{1},movies,'UniformOutput',false);
S = cellfun(@load,outFiles,'Unif',false);
S = [S{:}];
% F x N array, where F are frames, and N is the number of movies combined
% N == length(movieIndices)
movieInfoArray = [S.movieInfo];

% merge movies together
nFrames = size(movieInfoArray,1);
% nMovies = size(movieInfoArray,2);

fields = fieldnames(movieInfoArray);

combinedMovieInfo(nFrames) = movieInfoArray(1);

for frame = 1:nFrames
    for fieldNumber = 1:length(fields)
        field = fields{fieldNumber};
        combinedMovieInfo(frame).(field) = vertcat(movieInfoArray(frame,:).(field));
    end
end

end

