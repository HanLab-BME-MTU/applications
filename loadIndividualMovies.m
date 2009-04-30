function [exp] = loadIndividualMovies()
% runTracking engages the user in a dialog in order to choose a subset of
% movies under a given condition to analyze.
%
% SYNOPSIS [exp] = loadIndividualMovies(experiment))
%
% INPUT     experiment:  Structure containing source, date information, and
%                        frame rate for all movies under a given condition.
%
% OUTPUT    exp:  Structure containing the same information as the input
%                 structure but trimed down as specified by the user.
%
% REMARKS   The function lists all movies in experiment input structure
%           next to an identifier integer. It then asks the user to specify
%           whether to keep all movies, list movies to keep, or list movies
%           to take out.
%
% Daniel Nunez, March 19, 2008

%Load Data until user changes loadData value to 0
[experiment] = loadConditionData();
loadData = 1;
while loadData == 1
    loadData = input('Do you wish to load more condition folders?\nEnter 1 if yes, 0 if no.');
    if loadData == 1
    [exper] = loadConditionData();
    [experiment] = [experiment exper];
    end
end

% print out all movie paths with an identifying integer
for iMovie = 1:length(experiment);
    fprintf([int2str(iMovie) '\t' experiment(iMovie).source '\n'])
end

% ask user to choose in between tracking all movies, listing the movies
% that should be tracked, or listing the movies that should not be tracked
movies2Track = input('Enter:\n 1 to include all movies\n 2 to list the movies to include\n 3 to list the movies to not include\n');

%if user inputs 1 analyze all movies that need analysis
if movies2Track == 1
    exp = experiment;

%if user inputs 2 ask him to specify movies to analyzed (if they need to be
%analyzed)
elseif movies2Track == 2
    movies2Track = input('List identifying integers of movies as vector (inside square brackets separated by commas or spaces).\n');
    %delete all other movies from experiment structure
    for iMovie = 1:length(movies2Track)
        exp(iMovie) = experiment(movies2Track(iMovie));
    end %of for loop

%if user inputs 3 ask him to specify movies to leave out
elseif movies2Track == 3
    %ask user to input integers of movies to exclude
    movies2Track = input('List identifying integers of movies as vector (inside square brackets separated by commas or spaces).\n');
    %delete these movies from experiment structure if can find an entry in
    %movies2Track that matches the movie index number in the experiment
    %structure
    %for all movies
    count = 1; %initiate counter
    for iExperiment = 1:length(experiment)
        %if all entries in movies2Track fail to match experiment index of movie
        if isempty(find(movies2Track == iExperiment,1));
            %save movie 
            exp(count) = experiment(iExperiment);
            count = count + 1;
        end %of if loop
    end %of while loop
end


end %of function