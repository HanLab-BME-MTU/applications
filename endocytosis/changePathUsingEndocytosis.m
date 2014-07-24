function experiment = changePathUsingEndocytosis(experiment)

% changePathUsingEndocytosis changes the source path of a given movie by
% looking for the word endocytosis in the current working directory and
% substituting everything up to this word on movie directory
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
%remember current working dir
od = pwd;
%find the word 'endocytosis' on current working dir
findE = strfind(od,'endocytosis');
%grab dir structure up to this word
upstreamString = od(1:findE-1);
%find directory separators in source
findSeparators = strfind(upstreamString,'/');
findSeparators = [findSeparators strfind(upstreamString,'\')];
%change separators to filesep to make system os independent
upstreamString(findSeparators) = filesep;

%for each experiment
for iexp = 1:length(experiment)
%find endocytosis on on source path
findE = strfind(experiment(iexp).source,'endocytosis');
%grab string that denotes directory structure downstream of endocytosis
downstreamString = experiment(iexp).source(findE:end);
%find directory separators in source
findSeparators = strfind(downstreamString,'/');
findSeparators = [findSeparators strfind(downstreamString,'\')];
%change separators to filesep to make system os independent
downstreamString(findSeparators) = filesep;
%put paths together
experiment(iexp).source = strcat(upstreamString,downstreamString);


end %of for each experiment

end %of function