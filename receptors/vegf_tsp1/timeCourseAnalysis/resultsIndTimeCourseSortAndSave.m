function [ caseTimeList,caseResSummary, s ] = resultsIndTimeCourseSortAndSave( resSummary, caseParam, saveFile )
%resultsIndTimeCourseSortAndSave Go over each case and put its results together
% INPUT
%       resSummary: Struct array,
%                   see timeCourseAnalysis.util.emptyResSummaryStruct
%       caseParam: Structure array with number of entries = number of
%                  timecourses among which the movies will be divided. For
%                  each time course, it contains the following fields:
%           .indx    : Indices of movies in ML belonging to this
%                      timecourse.
%           .timeList: The time point of each movie, with order following
%                      the movie order in ML.
%           .name    : Case name, to be used in naming output variables.
%           .indx0min: Index of movie to be considered at time 0, e.g. when
%                      a ligand is added. If > 1, movies before will have a
%                      negative relative time, movies after wil have a
%                      positive relative time. If 1, relative time and
%                      absolute time are the same.
%       saveFile: logical of whether to save the file or not
%
%OUTPUT caseTimeList  : 2-column vector indicating the movie times. Column 1
%                       shows absolute time, Column 2 shows relative time.
%       caseResSummary: Structure array storing various results for each
%                       movie, ordered based on their time (corresponding
%                       to caseTimeList). Contains the fields: 
%           .diffSummary         : Diffusion analysis summary, as output
%                                  by summarizeDiffAnRes.
%           .diffCoefMeanPerClass: Mean diffusion coefficient per motion
%                                  class. Order: Immobile, confined, free,
%                                  directed, undetermined.
%           .confRadMeanPerClass : Mean confinement radius per motion
%                                  class. Same order as above.
%           .ampMeanPerClass     : Mean particle amplitude per motion
%                                  class. Rows in same order as above.
%                                  Columns: first for "absolute" amplitude,
%                                  second for amplitude normalized by
%                                  monomer amplitude, as derived from modal
%                                  analysis of particles in last 20 frames
%                                  of each movie.
%           .ampStatsF20         : Amplitude statistics in first 20
%                                  frame of each movie. 
%                                  Order: mean amplitude, first mode
%                                  mean, first mode std, first mode
%                                  fraction, number of modes, mean
%                                  normalized by monomer amplitude (see
%                                  ampMeanPerClass for details).
%           .ampStatsL20         : Same as ampStatsF20, but for last 20
%                                  frames of movie.
%        s             : struct containing the fields according to
%                        caseParam
%           .timeList_[caseName]
%           .resSummary_[caseName]
%
%Khuloud Jaqaman, March 2015
%Mark Kittisopikul, January 2016: Made into independent function

if nargin < 2
    caseParam = [];
end
if nargin < 3
    saveFile = true;
end

if(isempty(caseParam))
    numCases = 1;
else
    numCases = length(caseParam);
end

%go over each case and put its results together
for iCase = 1 : numCases
    if(~isempty(caseParam))
        %get parameters from input
        caseTimeList = caseParam(iCase).timeList;
        caseIndx = caseParam(iCase).indx;
        caseName = caseParam(iCase).name;
        caseMin0 = caseParam(iCase).indx0min;

        %sort time and add column for relative time
        offset = caseTimeList(caseMin0);
        [caseTimeList,indxSort] = sort(caseTimeList);
        caseTimeList = [caseTimeList caseTimeList-offset]; %#ok<AGROW>

        %collect and sort results
        caseResSummary = resSummary(caseIndx);
        caseResSummary = caseResSummary(indxSort);
    else
        caseTimeList = [];
        caseResSummary = resSummary;
        caseName = 'default';
    end

    %Saves by default
    if saveFile

        % mkitti: use struct assignment and save rather eval

        %name variables properly for saving
        s.(['timeList_' caseName]) = caseTimeList;
        s.(['resSummary_' caseName]) = caseResSummary;

        %save results
        file2save = fullfile(dir2save, ...
            s.(['resSummary_' caseName]), ...
            ['timeList_' caseName], ...
            ['resSummary_' caseName]);
        save(file2save,'-struct','s');

    end

end

end

