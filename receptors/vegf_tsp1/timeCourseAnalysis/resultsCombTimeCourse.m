function [resSummaryComb,resSummaryInd] = resultsCombTimeCourse(dsName,dsResFile,timeListComb,timeAbsOrRel)
%RESULTSCOMBTIMECOURSE compiles the results of multiple movies from multiple timecourses into one timecourse
%
%SYNOPSIS [resSummaryComb,resSummaryInd] = resultsCombTimeCourse(dsName,dsResFile,timeListComb,timeAbsOrRel)
%
%INPUT  dsName      : Cell array of the different timecourse names, needed to
%                     handle the different timecourse variables.
%       dsResFile   : Cell array indicating the name (including full path) of
%                     the results file for each individual timecourse.
%       timeListComb: Column vector with timecourse time points for
%                     combined dataset.
%       timeAbsOrRel: 'abs' or 'rel' to use, respectively, absolute or
%                     relative time to align and group the different datasets.
%                     The absolute and relative times are stored in the 
%                     results file of each dataset.
%                     Optional. Default: 'abs'.
%    
%OUTPUT resSummaryComb: Structure with results summary for combined
%                       timecourse. Contains the fields: 
%           .numAbsClass: Absolute number of particles in the various
%                         motion classes.
%                         Row = time points in time course (see timeList).
%                         Columns = immobile, confined, free, directed,
%                         undetermined, determined, total.
%           .numNorm0Class: Normalized number of particles in the various
%                         motion classes, such that the first time = 1.
%                         Rows and columns as numAbsClass.
%           .probClass  : Probability of the various motion classes.
%                         Rows as above.
%                         Columns = immobile, confined, free, directed (all
%                         relative to determined); and determined relative
%                         to total.
%           .diffCoefClass: Mean diffusion coefficient in the various
%                         motion classes.
%                         Rows as above.
%                         Columns = immobile, confined, free, directed,
%                         undetermined.
%           .confRadClass: Mean confinement radius in the various
%                         motion classes.
%                         Rows and columns as above.
%           .ampClass   : Mean amplitude of particles in the various motion
%                         classes.
%                         Rows and columns as above.
%           .ampNormClass:Mean normalized amplitude of particles in the
%                         various motion classes.
%                         Rows and columns as above.
%           .ampStatsF20: Amplitude statistics for particles in first 20
%                         frames.
%                         Rows as above.
%                         Columns = mean, first mode mean, first mode std,
%                         first mode fraction, number of modes, normalized
%                         mean.
%           .ampStatsL20: Amplitude statistics for particles in last 20
%                         frames.
%                         Rows and columns as above.
%           .rateMS     : Rate of merging and rate of splitting per
%                         feature. Columns: merging, splitting.
%           .timeList   : List of time points in combined time course. Same
%                         as input timeListComb.
%           .timeAbsOrRel:'abs' or 'rel' to indicate, respectively, whether
%                         absolute or relative time was used to align and
%                         combine datasets.
%       resSummaryInd : Structure array equivalent to resSummaryComb but
%                       for individual timecourses.
%
%Khuloud Jaqaman, March 2015

%% Input

if nargin < 3
    error('resultsCombTimeCourse: Too few input arguments')
end

if nargin < 4 || isempty(timeAbsOrRel)
    timeAbsOrRel = 'abs';
end

if strcmp(timeAbsOrRel,'abs')
    timeIndxComb = 1;
else
    timeIndxComb = 2;
end

%get number of datasets
numDS = length(dsName);

%% Individual dataset results

%reserve memory
[dsSummary,timeList,numAbsClassInd,numNorm0ClassInd,probClassInd,...
    diffCoefClassInd,confRadClassInd,ampClassInd,ampNormClassInd,...
    ampStatsF20Ind,ampStatsL20Ind,rateMSInd] = deal(cell(numDS,1));

%load all results and store in cell array
for iDS = 1 : numDS
    load(dsResFile{iDS});
    varName = ['resSummary_' dsName{iDS} ];
    eval(['dsSummary{iDS} = ' varName ';']);
    varName = ['timeList_' dsName{iDS} ];
    eval(['timeList{iDS} = ' varName ';']);
end

for iDS = 1 : numDS
    
    %absolute numbers of molecules in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet, 6 det, 7 tot
    %rows are for different timepoints, as listed in corresponding timeList
    diffSummary = vertcat(dsSummary{iDS}.diffSummary); %#ok<NASGU>
    numAbsClassCurrent = [...
        catStruct(1,'diffSummary.probMotionType(6,3)') ...
        catStruct(1,'diffSummary.probMotionType(7,3)') ...
        catStruct(1,'diffSummary.probMotionType(8,3)') ...
        catStruct(1,'diffSummary.probMotionType(9,3)') ...
        catStruct(1,'diffSummary.probMotionType(11,3)')];
    numAbsClassCurrent = [numAbsClassCurrent sum(numAbsClassCurrent(:,1:4),2) sum(numAbsClassCurrent(:,1:5),2)]; %#ok<AGROW>
    numAbsClassInd{iDS} = numAbsClassCurrent;
        
    %numbers in various motion classes normalized to 1 at absolute time 0
    %columns and rows same as absolute numbers
    numNorm0ClassInd{iDS} = numAbsClassCurrent ./ repmat(numAbsClassCurrent(1,:),length(timeList{iDS}),1);
    
    %absolute probabilities of various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir (all relative to det), 
    %5 det (relative to tot)
    probClassInd{iDS} = [numAbsClassCurrent(:,1:4)./repmat(numAbsClassCurrent(:,6),1,4) ...
        numAbsClassCurrent(:,6)./numAbsClassCurrent(:,7)];
    
    %diffusion coefficient in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    diffCoefClassInd{iDS} = (horzcat(dsSummary{iDS}.diffCoefMeanPerClass))';
    
    %confinement radius in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    confRadClassInd{iDS} = (horzcat(dsSummary{iDS}.confRadMeanPerClass))';

    %amplitude in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    tmpCollect = (horzcat(dsSummary{iDS}.ampMeanPerClass))';
    ampClassInd{iDS} = tmpCollect(1:2:end,:);
    ampNormClassInd{iDS} = tmpCollect(2:2:end,:);

    %amplitude statistics in first 20 frames
    %columns are: 1 mean, 2 first mode mean, 3 first mode std, 4 first mode fraction, 5 number of modes
    ampStatsF20Ind{iDS} = vertcat(dsSummary{iDS}.ampStatsF20);
    
    %amplitude statistics in last 20 frames
    %columns are: 1 mean, 2 first mode mean, 3 first mode std, 4 first mode fraction, 5 number of modes
    ampStatsL20Ind{iDS} = vertcat(dsSummary{iDS}.ampStatsL20);
    
    %rate of merging and rate of splitting
    tmp = vertcat(dsSummary{iDS}.statsMS);
    rateMSInd{iDS} = tmp(:,6:7);
    
end

%% Combine dataset results

%do direct combination for all variables except for
%normalized number of particles in various classes
[numAbsClassComb.msn,numAbsClassComb.sample] = combineWithTimeAlign(numAbsClassInd,timeList,timeListComb,timeIndxComb);
[probClassComb.msn,probClassComb.sample] = combineWithTimeAlign(probClassInd,timeList,timeListComb,timeIndxComb);
[diffCoefClassComb.msn,diffCoefClassComb.sample] = combineWithTimeAlign(diffCoefClassInd,timeList,timeListComb,timeIndxComb);
[confRadClassComb.msn,confRadClassComb.sample] = combineWithTimeAlign(confRadClassInd,timeList,timeListComb,timeIndxComb);
[ampClassComb.msn,ampClassComb.sample] = combineWithTimeAlign(ampClassInd,timeList,timeListComb,timeIndxComb);
[ampNormClassComb.msn,ampNormClassComb.sample] = combineWithTimeAlign(ampNormClassInd,timeList,timeListComb,timeIndxComb);
[ampStatsF20Comb.msn,ampStatsF20Comb.sample] = combineWithTimeAlign(ampStatsF20Ind,timeList,timeListComb,timeIndxComb);
[ampStatsL20Comb.msn,ampStatsL20Comb.sample] = combineWithTimeAlign(ampStatsL20Ind,timeList,timeListComb,timeIndxComb);
[rateMSComb.msn,rateMSComb.sample] = combineWithTimeAlign(rateMSInd,timeList,timeListComb,timeIndxComb);

%for normalized number of particles, normalize by mean value at time 0 for
%combined datasets
tmp = numAbsClassComb.sample;
tmp = tmp ./ repmat(numAbsClassComb.msn(1,:,1),[size(tmp,1) 1 size(tmp,3)]);
numNorm0ClassComb.sample = tmp;
numNorm0ClassComb.msn = numAbsClassComb.msn;
numNorm0ClassComb.msn(:,:,1) = nanmean(tmp,3);
numNorm0ClassComb.msn(:,:,2) = nanstd(tmp,[],3);

%% Output

resSummaryInd = struct('numAbsClass',numAbsClassInd,'numNorm0Class',numNorm0ClassInd,...
    'probClass',probClassInd,'diffCoefClass',diffCoefClassInd,...
    'confRadClass',confRadClassInd,'ampClass',ampClassInd,'ampNormClass',ampNormClassInd,...
    'ampStatsF20',ampStatsF20Ind,'ampStatsL20',ampStatsL20Ind,'rateMS',...
    rateMSInd,'timeList',timeList);

resSummaryComb = struct('numAbsClass',numAbsClassComb,'numNorm0Class',numNorm0ClassComb,...
    'probClass',probClassComb,'diffCoefClass',diffCoefClassComb,...
    'confRadClass',confRadClassComb,'ampClass',ampClassComb,'ampNormClass',ampNormClassComb,...
    'ampStatsF20',ampStatsF20Comb,'ampStatsL20',ampStatsL20Comb,...
    'rateMS',rateMSComb,'timeList',timeListComb,'timeAbsOrRel',timeAbsOrRel);


%% ~~~ the end ~~~
