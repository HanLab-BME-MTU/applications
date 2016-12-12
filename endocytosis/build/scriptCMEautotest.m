
if strcmp(computer('arch'),'win64')
	data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
else
	data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
end	

% Cell mask -- get GUI program 
% cell mask draw override function (tweak)
% cellViewer(data(1), 'mode', 'mask') -- 

cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sortTiffStacks; % helper function for organizing data structure


data = loadConditionData(cmeData_dir, {''}, {'EGFP'}); %[works]
cmeAnalysis(data);  %[works] 
%  cmeAnalysis(data, 'Overwrite', [true false false true]) -- manually go
%  over tracks for semgnetation mask.
ccpSorter(data); % xinxin [works]
cmeDataViewer(data(1)); % xinxin [works]
cmeDataViewer(data(2)); % xinxin [works]
analyzeBleaching(data(1)); % one movie at a time.
analyzeBleaching(data(2)); % one movie at a time.
[lftRes, res] = runLifetimeAnalysis(data);  % works - error bar fixes

plotLifetimes(lftRes);  % works

track_num = 1;
movie_num = 1;

track = loadTracks(data(movie_num)); % xinxin 
plotTrack(data(movie_num), track(track_num));
plotIntensityCohorts(data); % works (good)
% plotInitiationDensity({da  lftRes(1),lftRes(2)}); % ask Marcel 
plotInitiationDensity({results.lftRes}, {'ctrl'}) % error bar MATLAB version issue
plotCellArea({results.lftRes}, {'ctrl'}); % ask Marcel (sub-output in CME)
% plotInitIntensityVsLifetime(data);  % ???? not used regularly
plotIntensityDistributions(data(movie_num)); % works
plotLifetimeComparison({lftRes}); % wrap with cell % ask Marcel? (need at least two data)
plotMaxIntensityDistribution(data);  % work, but old % MATLAB update fixes

% get two condition test and two movie

data = loadConditionData(epiTIRFData_dir, {'TIRF', 'EPI'}, {'EGFP', 'EGFP'});
% master, slave analysis -- ask Marcel
% epiTIRFAnalysis(data); % not used by xinxin yet -- % Ask Marcel [need special dataset]
epiTIRFAnalysis({data, data},  {'EGFP', 'EGFP'}, 9.5); % not used by xinxin yet -- % Ask Marcel [need special dataset]
% ask about LabArchive -- documentation on CME analysis