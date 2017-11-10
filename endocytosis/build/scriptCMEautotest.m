% Now using MATLAB-based testing/performance suite
% Andrew R. Jamieson 2016
% Test cmeAnalysis


% function scriptCMEautotest()

% preconditions
if strcmp(computer('arch'),'win64')
	data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
else
	data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
end	

% Cell mask -- get GUI program 
% cell mask draw override function (tweak)
% cellViewer(data(1), 'mode', 'mask') -- 

cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
cmeData_dir2 = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_mu_PIP2_mutant']);

epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);
epiTIRFData_dir2 = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_mu_PIP2_mutant']);


% - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sortTiffStacks; % helper function for organizing data structure

% Setup data to use loadConditionData
data = loadConditionData(cmeData_dir, {''}, {'EGFP'}); %[works]
data2 = loadConditionData(cmeData_dir2, {''}, {'EGFP'}); %[works]

%% RUN__cmeAnalysis__Data1
results = cmeAnalysis(data);  %[works]
results2 = cmeAnalysis(data2);  %[works]
plotCellArea({results.lftRes}, {'ctrl'}); % (included in cmeAnalysis)
plotInitiationDensity({results.lftRes, results2.lftRes}, {'ctrl', 'other'}) % error bar styling MATLAB version update

% ---
%  cmeAnalysis(data, 'Overwrite', [true false false true]) -- manually go
%  over tracks for semgnetation mask.
% ---

%% Run ccpSorter
ccpSorter(data); % xinxin [works]
ccpSorter(data2); % xinxin [works]

%% Run cmeDataViewer
cmeDataViewer(data(1)); % xinxin [works]
cmeDataViewer(data(2)); % xinxin [works]

%% analyzeBleaching
analyzeBleaching(data(1)); % one movie at a time.
analyzeBleaching(data(2)); % one movie at a time.

%% runLifetimeAnalysis__and__Plots
[lftRes, res] = runLifetimeAnalysis(data);  % works - error bar fixes
[lftRes2, res2] = runLifetimeAnalysis(data2);  % works - error bar fixes
plotLifetimes(lftRes);  % works (included in cmeAnalysis)
plotLifetimes(lftRes2);  % works (included in cmeAnalysis)
plotLifetimeComparison({lftRes, lftRes2}, {'con1','con2'}); % wrap with cell % ask Marcel? (need at least two data)

track_num = 1;
movie_num = 1;

track = loadTracks(data(movie_num)); % xinxin 
plotTrack(data(movie_num), track(track_num)); %  % (included in cmeAnalysis)
plotIntensityCohorts(data); % (included in cmeAnalysis)
plotMaxIntensityDistribution(data);  % works, but old % MATLAB update fixes
plotIntensityDistributions(data(movie_num)); % works





% %%%
% plotInitIntensityVsLifetime(data);  % ???? not used regularly ????
% %%%

%% epiTIRFAnalysis
dataEpi = loadConditionData(epiTIRFData_dir2, {'TIRF', 'EPI'}, {'EGFP', 'EGFP'});
epiTIRFAnalysis({dataEpi, dataEpi},  {'EGFP', 'EGFP'}, 9.5); % Still in development