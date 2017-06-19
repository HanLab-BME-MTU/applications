% Instructions to run CellX or Dope

% shortcuts cheatsheet:

% Make sure figure is in foreground (i.e., click figure if not working)

% 't' : toggle time limit display 
% 'e' : export handle/data info  (handleCX, dataCX)
% 's' : toggle Stage drift correction selector
% 'f' : toggle frame rate selector
% 'l' : toggle display of MovieData information at top

% This includes ALL Gen2 and Gen 3, including cell lines and melanocytes (non-malignant)
% i.e., meteff = NaN,0,1
cellXploreDR('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat')
dopeAnnotator('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat')

% This set is a subset, including only "primary metastatic cells" meteff = 0 or 1
% i.e., meteff = 0,1

% tic
% load('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat');
% newcellDataSet = {};
% parfor i = 1:length(cellDataSet)

%     icell = cellDataSet{i};
%     if ~isnan(icell.metEff)
%         newcellDataSet = [newcellDataSet {icell}];
%     end
% end
% cellDataSet = newcellDataSet;
% matFileName = '/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May16_metEffOnly.mat';
% save(matFileName, 'cellDataSet', 'annotationSet');
% toc


cellXploreDR('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May16_metEffOnly.mat')
dopeAnnotator('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May16_metEffOnly.mat')


%% load with previous annotations
dopeAnnotator('/work/bioinformatics/shared/dope/data/OMETIFF/Gen2n3_May15_ALL.mat',...
	'/work/bioinformatics/shared/dope/data/OMETIFF/dopeCheckPoint_Anno09Jun2017_1037.mat')