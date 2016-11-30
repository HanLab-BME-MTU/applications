
if strcmp(computer('arch'),'win64')
	data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
else
	data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
end	

cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = loadConditionData(cmeData_dir, {''}, {'EGFP'}); %[works]
cmeAnalysis(data);  %[works] 
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
plotIntensityCohorts(data); % works
plotInitiationDensity(lftRes); % ask Marcel
plotCellArea(lftRes); % ask Marcel
plotInitIntensityVsLifetime(data);  % ask Marcel
plotIntensityDistributions(data(movie_num)); % works
plotLifetimeComparison({lftRes}); % wrap with cell % ask Marcel?
plotMaxIntensityDistribution(data);  % work, but old

sortTiffStacks; % ??? % Ask Marcel
data = loadConditionData(epiTIRFData_dir, {'TIRF', 'EPI'}, {'EGFP', 'EGFP'});
% master, slave analysis -- ask Marcel
epiTIRFAnalysis(data); % not used by xinxin yet -- % Ask Marcel

% ask about LabArchive -- documentation on CME analysis