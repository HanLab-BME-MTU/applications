function [resSummaryInd] = resultsCombTimeCourseMod(dsSummary)
%RESULTSCOMBTIMECOURSE compiles the results of multiple movies from multiple timecourses into one timecourse
%
%SYNOPSIS [resSummaryComb,resSummaryInd] = resultsCombTimeCourseFunction(dsSummary,timeList,timeListComb,timeAbsOrRel)
%
%INPUT  dsSummary   : Structure array storing various results for each
%                     movie, ordered based on their time (corresponding
%                     to caseTimeList). See 'resultsIndTimeCourse'.
%       timeList    : 2-column vector indicating the movie times. Column 1
%                     shows absolute time, Column 2 shows relative time.
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
%                         feature (per time unit per number of features).
%                         Columns: merging, splitting.
%           .timeList   : List of time points in combined time course. Same
%                         as input timeListComb.
%           .timeAbsOrRel:'abs' or 'rel' to indicate, respectively, whether
%                         absolute or relative time was used to align and
%                         combine datasets.
%       resSummaryInd : Structure array equivalent to resSummaryComb but
%                       for individual timecourses.
%
%Khuloud Jaqaman, March 2015
%modified from resultsCombTimeCourse
%Tae H Kim, June 2015

%% Input




%get number of datasets
%numDS = 1;

%% Individual dataset results

%reserve memory
%{
[numAbsClassInd,numNorm0ClassInd,probClassInd,...
    diffCoefClassInd,confRadClassInd,ampClassInd,ampNormClassInd,...
    ampStatsF20Ind,ampStatsL20Ind,rateMSInd] = deal(cell(numDS,1));
%}
for iDS = 1
    
    %absolute numbers of molecules in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet, 6 det, 7 tot, 8 imm+conf
    %rows are for different timepoints, as listed in corresponding timeList
    diffSummary = vertcat(dsSummary.diffSummary); %#ok<NASGU>
    numAbsClassCurrent = [...
        catStruct(1,'diffSummary.probMotionType(6,3)') ...
        catStruct(1,'diffSummary.probMotionType(7,3)') ...
        catStruct(1,'diffSummary.probMotionType(8,3)') ...
        catStruct(1,'diffSummary.probMotionType(9,3)') ...
        catStruct(1,'diffSummary.probMotionType(11,3)')];
    numAbsClassCurrent = [numAbsClassCurrent sum(numAbsClassCurrent(:,1:4),2) ...
        sum(numAbsClassCurrent(:,1:5),2) sum(numAbsClassCurrent(:,1:2),2)]; %#ok<AGROW>
    numAbsClassInd = numAbsClassCurrent;
        
    %numbers in various motion classes normalized to 1 at absolute time 0
    %columns and rows same as absolute numbers
    numNorm0ClassInd = numAbsClassCurrent ./ repmat(numAbsClassCurrent(1,:),length(dsSummary),1);
    
    %absolute probabilities of various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir (all relative to det), 
    %5 det (relative to tot), 6 imm+conf (relative to det)
    probClassInd = [numAbsClassCurrent(:,1:4)./repmat(numAbsClassCurrent(:,6),1,4) ...
        numAbsClassCurrent(:,6)./numAbsClassCurrent(:,7) ...
        numAbsClassCurrent(:,8)./numAbsClassCurrent(:,6)];
    
    %cell area in order to get densities
    cellAreaInd = repmat(vertcat(dsSummary.cellArea),1,8);
    
    %densities of molecules in various motion classes
    %columns and rows same as absolute numbers
    densityAbsClassInd = numAbsClassInd ./ cellAreaInd;
    
    %densities in various motion classes normalized to 1 at absolute time 0
    %columns and rows same as absolute numbers
    densityNorm0ClassInd = densityAbsClassInd ./ repmat(densityAbsClassInd(1,:),length(dsSummary),1);
    
    %diffusion coefficient in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    diffCoefClassInd = (horzcat(dsSummary.diffCoefMeanPerClass))';
    
    %confinement radius in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    confRadClassInd = (horzcat(dsSummary.confRadMeanPerClass))';

    %amplitude in various motion classes
    %columns are: 1 imm, 2 conf, 3 free, 4 dir, 5 undet
    tmpCollect = (horzcat(dsSummary.ampMeanPerClass))';
    ampClassInd = tmpCollect(1:2:end,:); %might want: ampClassInd = tmpCollect(1:5,:);
    ampNormClassInd = tmpCollect(2:2:end,:); %and: ampNormClassInd = tmpCollect(6:10,:); instead?

    %amplitude statistics in first 20 frames
    %columns are: 1 mean, 2 first mode mean, 3 first mode std, 4 first mode fraction, 5 number of modes
    ampStatsF20Ind = vertcat(dsSummary.ampStatsF20);
    
    %amplitude statistics in last 20 frames
    %columns are: 1 mean, 2 first mode mean, 3 first mode std, 4 first mode fraction, 5 number of modes
    ampStatsL20Ind = vertcat(dsSummary.ampStatsL20);
    
    %rate of merging and rate of splitting
    tmp = vertcat(dsSummary.statsMS);
    rateMSInd = tmp(:,8:9);
    
end

%% Combine dataset results
%{
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
%}
%% Output
%{
%Replace all nan with 0
numAbsClassInd(isnan(numAbsClassInd)) = 0;
numNorm0ClassInd(isnan(numNorm0ClassInd)) = 0;
probClassInd(isnan(probClassInd)) = 0;
diffCoefClassInd(isnan(diffCoefClassInd)) = 0;
confRadClassInd(isnan(confRadClassInd)) = 0;
ampClassInd(isnan(ampClassInd)) = 0;
ampNormClassInd(isnan(ampNormClassInd)) = 0;
ampStatsF20Ind(isnan(ampStatsF20Ind)) = 0;
ampStatsL20Ind(isnan(ampStatsL20Ind)) = 0;
%}

resSummaryInd = struct('numAbsClass',numAbsClassInd,'numNorm0Class',numNorm0ClassInd,...
    'densityAbsClass',densityAbsClassInd,'densityNorm0Class',densityNorm0ClassInd,...
    'probClass',probClassInd,'diffCoefClass',diffCoefClassInd,'confRadClass',confRadClassInd,...
    'ampClass',ampClassInd,'ampNormClass',ampNormClassInd,...
    'ampStatsF20',ampStatsF20Ind,'ampStatsL20',ampStatsL20Ind,'rateMS',rateMSInd);

%% ~~~ the end ~~~
