
if strcmp(computer('arch'),'win64')
	data_root = 'C:\Users\Andrew\Data\raw\CME\xinxin\auto_test\';
else
	data_root = '/work/bioinformatics/s170480/Data/CME/xinxin/auto_test/';
end	

cmeData_dir = fullfile(data_root, ['single_channel' filesep 'Zuzana_ARPE_CLC_egfp_07032014_control']);
epiTIRFData_dir = fullfile(data_root, ['epi_tirf' filesep 'Zuzana_arpe_clc_egfp_14082015_control']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = loadConditionData(cmeData_dir, {''}, {'EGFP'});
cmeAnalysis(data);
ccpSorter(data);
cmeDataViewer(data(1));
cmeDataViewer(data(2));
analyzeBleaching(data);
lftRes = runLifetimeAnalysis(data);
track = loadTracks(data(1));


here is my change

plotLifetimes(lftRes);
plotTrack(data(1), track(1));
plotIntensityCohorts(data);
plotInitiationDensity(lftRes);
plotCellArea(lftRes);
plotInitIntensityVsLifetime(data);
plotIntensityDistributions(data);
plotLifetimeComparison(lftRes);
plotLifetimes(lftRes);
plotMaxIntensityDistribution(data);

sortTiffStacks; % ???
data = loadConditionData(epiTIRFData_dir, {'TIRF', 'EPI'}, {'EGFP', 'EGFP'});
epiTIRFAnalysis(data);