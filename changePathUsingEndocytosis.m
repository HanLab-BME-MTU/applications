function experiment = changePathUsingEndocytosis(experiment)

% changePathUsingEndocytosis changes the sourc path of a given movie by
% looking for the word endocytosis in the current working directory and
% substituting everything up to this word
%
% INPUT     experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%
% OUTPUT
%
% REMARKS   
%
% Uses: none
%
% Daniel Nunez, October 23, 2009

%find endocytosis on cwd
od = pwd;
findE = strfind(od,'endocytosis');
directoryString = od(1:findE-1);

%for each experiment
for iexp = 1:length(experiment)
%find endocytosis on on source path
findE = strfind(experiment(iexp).source,'endocytosis');
%change begining of source path with current path up to endocytosis
experiment(iexp).source = strcat(directoryString,experiment(iexp).source(findE:end));


end %of for each experiment

end %of function