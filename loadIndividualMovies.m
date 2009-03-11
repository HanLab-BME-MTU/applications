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

%list all movies and ask user to select movies to use in analysis
[selection, selectionList] = listSelectGUI({experiment(:).source},[],'move',[]);

exp = experiment(selection);

end %of function