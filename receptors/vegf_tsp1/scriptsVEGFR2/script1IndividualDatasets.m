%% 141028

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2014/20141028_HMVECp6-R2-488/!AllmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_mAAL_HMVEC_141028';
caseParam(1).timeList = [0 11 116 21 31 41 51 61 64 6]';
caseParam(1).indx = (1:10);
caseParam(1).indx0min = 1;

resultsIndTimeCourse(ML,caseParam);

%% 141030

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2014/20141030_HMVECp7_R2_1ugML-AF488/AllmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_mAAL_HMVEC_141030';
caseParam(1).timeList = [0 2 12 22 32 42 7 61]';
caseParam(1).indx = (1:8);
caseParam(1).indx0min = 1;

caseParam(2).name = 'pVEGF_mAAL_HMVEC_141030';
caseParam(2).timeList = [0 2 13 12 22 32 42 52 6 60]';
caseParam(2).indx = (9:18);
caseParam(2).indx0min = 2;

resultsIndTimeCourse(ML,caseParam);

%% 141031

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2014/20141031_HMVECp7-R2_half-ugML-488_original/!All-movieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_mAAL_HMVEC_141031';
caseParam(1).timeList = [0 1 10 20 31 40 50 5 62 72]';
caseParam(1).indx = (10:19);
caseParam(1).indx0min = 1;

caseParam(2).name = 'pVEGF_mAAL_HMVEC_141031';
caseParam(2).timeList = [0  2 11 21 31 41 51 64 62]';
caseParam(2).indx = (1:9);
caseParam(2).indx0min = 2;

resultsIndTimeCourse(ML,caseParam);

%% 150128

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150128_NGSBlockedHMVECp6_EIC1ug_FabAF488/TotalmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'pVEGF_mAAL_HMVEC_150128';
caseParam(1).timeList = [0 16 19 27 40 47 57 67 32 87]';
caseParam(1).indx = (1:10);
caseParam(1).indx0min = 4;

caseParam(2).name = 'pVEGF_pAAL_HMVEC_150128';
caseParam(2).timeList = [0 4 7 9 14 18 24 26 34 44 55 65 84 89]';
caseParam(2).indx = (11:24);
caseParam(2).indx0min = 5;

resultsIndTimeCourse(ML,caseParam);

%% 150129

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150129_HMVECp6_EICFab0.2ug/!noAAL-movieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'pVEGF_mAAL_HMVEC_150129';
caseParam(1).timeList = [0 3 7 13 17 22 26 32 42 52 70]';
caseParam(1).indx = (1:11);
caseParam(1).indx0min = 4;

resultsIndTimeCourse(ML,caseParam);

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150129_HMVECp6_EICFab0.2ug/!wAAL-movieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'pVEGF_pAAL_HMVEC_150129';
caseParam(1).timeList = [0 2 4 8 12 16 26 28 35 45 62 64]';
caseParam(1).indx = (1:12);
caseParam(1).indx0min = 4;

resultsIndTimeCourse(ML,caseParam);

%% 150211 - HUVECs

%% 150212 - HUVECs

%% 150216A - HUVECs

%% 150216B - HUVECs

%% 150217A - HUVECs

%% 150217B - HUVECs

%% 150218 - HUVECs

%% 150224

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150224_HMVECp6_VEGF/!AllmovieList20150224.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_mAAL_HMVEC_150224';
caseParam(1).timeList = [0 10 20 31 40 50 5 59 70]';
caseParam(1).indx = (1:9);
caseParam(1).indx0min = 1;

caseParam(2).name = 'pVEGF_mAAL_HMVEC_150224';
caseParam(2).timeList = [0 10 14 20 1 30 39 49 59 69 6]';
caseParam(2).indx = [10:17 19:21];
caseParam(2).indx0min = 2;

resultsIndTimeCourse(ML,caseParam);

%% 150225

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150225_HMVECp6_AAL_VEGF/!AllmovieList20150225.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_pAAL_HMVEC_150225';
caseParam(1).timeList = [0 9 15 25 30 3 39 49 60 69 6]';
caseParam(1).indx = (1:11);
caseParam(1).indx0min = 1;

caseParam(2).name = 'pVEGF_pAAL_HMVEC_150225';
caseParam(2).timeList = [0 13 18 28 2 38 48 58 5 67 8]';
caseParam(2).indx = (12:22);
caseParam(2).indx0min = 11;

resultsIndTimeCourse(ML,caseParam);

%% 150302

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150302_HMVECp6_VEGF/!noVEGFAllmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_mAAL_HMVEC_150302';
caseParam(1).timeList = [0 13 16 21 26 2 36 46 3 56 76 78 5]';
caseParam(1).indx = (1:13);
caseParam(1).indx0min = 1;

resultsIndTimeCourse(ML,caseParam);

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150302_HMVECp6_VEGF/!VEGFAllmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'pVEGF_mAAL_HMVEC_150302';
caseParam(1).timeList = [0 11 16 26 2 36 4 46 56 66 7]';
caseParam(1).indx = (1:11);
caseParam(1).indx0min = 11;

resultsIndTimeCourse(ML,caseParam);

%% 150303

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150303_HMVECp6_AAL/!AAL_NoVEGFallmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'mVEGF_pAAL_HMVEC_150303';
caseParam(1).timeList = [0 10 15 21 30 41 5 51 7 70]';
caseParam(1).indx = (1:10);
caseParam(1).indx0min = 1;

resultsIndTimeCourse(ML,caseParam);

clear all

movieListFile = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/2015/20150303_HMVECp6_AAL/!AAL_VEGFallmovieList.mat';
ML = MovieList.load(movieListFile);

caseParam(1).name = 'pVEGF_pAAL_HMVEC_150303';
caseParam(1).timeList = [0 10 15 17 21 2 30 40 50 5 70]';
caseParam(1).indx = (1:11);
caseParam(1).indx0min = 2;

resultsIndTimeCourse(ML,caseParam);

