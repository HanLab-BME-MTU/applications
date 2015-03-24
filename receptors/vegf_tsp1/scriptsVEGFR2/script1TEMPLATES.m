
%% TEMPLATE 1: When both datasets (e.g. p/m VEGF) are in same movieList

%see 141030 and 141031 datasets as examples

%MAKE A DIRECTORY IN TOP DORECTORY CALLED analysisKJ

clear all

%first load new movieList
movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/....mat';
ML = MovieList.load(movieListFile);

%then run the below line in order to get the movie order in movie list
%you need this to fill in the variables below
for iM = 1 : length(ML.movieDataFile_), disp([num2str(iM) '  ' ML.movieDataFile_{iM}]); end

%then copy the commands below to the appropriate place in
%script1IndividualDatasets and fill the caseParam fields as described

%START COPY

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/....mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = ''; %name should be equivalent to: mVEGF_mAAL_HMVEC_150224. Simply change the p or m, HMVEC or HUVEC, and date in same format
caseParam(1).timeList = []'; %time list in same order as movies in movieList
caseParam(1).indx = (); %index of movies in movieList
caseParam(1).indx0min = ; %index of movie in movieList that is considered 0min after VEGF. If dataset has no VEGF, enter 1

caseParam(2).name = ''; %these entries are same as above, just for the other dataset
caseParam(2).timeList = []';
caseParam(2).indx = [];
caseParam(2).indx0min = ;

resultsIndTimeCourse(ML,caseParam);

%END COPY

%% TEMPLATE 2: When each dataset has its own movieList

%see 150302 and 150303 as examples

%MAKE A DIRECTORY IN TOP DORECTORY CALLED analysisKJ

clear all

%first load one movieList
movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150307_HMVECp6/!0307VEGFmovieList.mat';
ML = MovieList.load(movieListFile);

%then run the below line in order to get the movie order in movie list
%you need this to fill in the variables below
for iM = 1 : length(ML.movieDataFile_), disp([num2str(iM) '  ' ML.movieDataFile_{iM}]); end

%then copy the commands below to the appropriate place in
%script1IndividualDatasets and fill the caseParam fields as described

%START COPY

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150307_HMVECp6/!0307VEGFmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'pVEGF_mAAL_HMVEC_150307'; %see comments above
caseParam(1).timeList = [0 15 20 2 30 40 51 5 70 8]';
caseParam(1).indx = (1:10);
caseParam(1).indx0min = 10;

resultsIndTimeCourse(ML,caseParam);

%END COPY

%then repeat the above for the second movieList
