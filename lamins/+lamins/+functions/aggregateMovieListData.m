function [ data ] = aggregateMovieListData( ML , filename)
%aggregateMovieListData Pulls together movieList data
% ML MovieList
% filename to load

if(nargin < 2)
    filename = 'skeletons_2015_06_10.mat';
end

if(isempty(ML.movies_))
    ML.sanityCheck();
end

directories = cellfun(@(MD) MD.outputDirectory_,ML.movies_,'Unif',false)';
data = cellfun(@(D) load([D filesep filename]),directories,'Unif',false);
data = [data{:}];

Sf = cellfun(@(S,tz) S{(tz(1)-1)*length(tz)+1},{data.S3},{data.tz},'Unif',false);
[data.Sf] = Sf{:};

if(length(data(1).tz) == 2)
Sf2 = cellfun(@(S,tz) S{(tz(2)-1)*length(tz)+2},{data.S3},{data.tz},'Unif',false,'ErrorHandler',@(a,b,c) []);
[data.Sf2] = Sf2{:};
end

end

