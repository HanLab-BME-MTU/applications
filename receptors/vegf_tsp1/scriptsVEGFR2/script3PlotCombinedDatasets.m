
%% Common variables between all plots

%time list
timeList{1} = timeListComb_mVEGF;
timeList{2} = timeListComb_mVEGF;
timeList{3} = timeListComb_pVEGF;
timeList{4} = timeListComb_pVEGF;

%condition names
condName = {'-VEGF -AAL','-VEGF +AAL','+VEGF -AAL','+VEGF +AAL'};

%condition colors
condColor = {'k','b','r','g'};

%flag whether to shift negative time so that first time point for all
%datasets is 0
shiftNegTime = 1;

%directory for saving figures
dir2save = '/project/biophysics/jaqaman_lab/vegf_tsp1/slee/VEGFR2/analysisCombinedDatasets/150320_analysis/figures';

%% Absolute number of molecules in various motion classes

tmp{1} = resSummaryCombMM.numAbsClass.msn;
tmp{2} = resSummaryCombMP.numAbsClass.msn;
tmp{3} = resSummaryCombPM.numAbsClass.msn;
tmp{4} = resSummaryCombPP.numAbsClass.msn;

figNameList = {'numAbsImm','numAbsConf','numAbsFree','numAbsDir',...
    'numAbsUndet','numAbsDet','numAbsTot'};

yaxisUnits = '(molecules)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Normalized number of molecules in various motion classes

tmp{1} = resSummaryCombMM.numNorm0Class.msn;
tmp{2} = resSummaryCombMP.numNorm0Class.msn;
tmp{3} = resSummaryCombPM.numNorm0Class.msn;
tmp{4} = resSummaryCombPP.numNorm0Class.msn;

figNameList = {'numNorm0Imm','numNorm0Conf','numNorm0Free','numNorm0Dir',...
    'numNorm0Undet','numNorm0Det','numNorm0Tot'};

yaxisUnits = '(unitless)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Probability of various motion classes

tmp{1} = resSummaryCombMM.probClass.msn;
tmp{2} = resSummaryCombMP.probClass.msn;
tmp{3} = resSummaryCombPM.probClass.msn;
tmp{4} = resSummaryCombPP.probClass.msn;

figNameList = {'probImm','probConf','probFree','probDir','probDet'};

yaxisUnits = '(unitless)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Diffusion coefficient in various motion classes

tmp{1} = resSummaryCombMM.diffCoefClass.msn;
tmp{2} = resSummaryCombMP.diffCoefClass.msn;
tmp{3} = resSummaryCombPM.diffCoefClass.msn;
tmp{4} = resSummaryCombPP.diffCoefClass.msn;

figNameList = {'diffCoefImm','diffCoefConf','diffCoefFree','diffCoefDir','diffCoefUndet'};

yaxisUnits = '(pixels^2/frame)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% Confinement radius in various motion classes

tmp{1} = resSummaryCombMM.confRadClass.msn;
tmp{2} = resSummaryCombMP.confRadClass.msn;
tmp{3} = resSummaryCombPM.confRadClass.msn;
tmp{4} = resSummaryCombPP.confRadClass.msn;

figNameList = {'confRadImm','confRadConf','confRadFree','confRadDir','confRadUndet'};

yaxisUnits = '(pixels)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% amplitude in various motion classes

tmp{1} = resSummaryCombMM.ampClass.msn;
tmp{2} = resSummaryCombMP.ampClass.msn;
tmp{3} = resSummaryCombPM.ampClass.msn;
tmp{4} = resSummaryCombPP.ampClass.msn;

figNameList = {'ampImm','ampConf','ampFree','ampDir','ampUndet'};

yaxisUnits = '(arbitrary units)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% amplitude in first and last 20 frames

tmp{1} = resSummaryCombMM.ampFL20.msn;
tmp{2} = resSummaryCombMP.ampFL20.msn;
tmp{3} = resSummaryCombPM.ampFL20.msn;
tmp{4} = resSummaryCombPP.ampFL20.msn;

figNameList = {'ampF20','ampL20'};

yaxisUnits = '(arbitrary units)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

%% rates of merging and splitting

tmp{1} = resSummaryCombMM.rateMS.msn;
tmp{2} = resSummaryCombMP.rateMS.msn;
tmp{3} = resSummaryCombPM.rateMS.msn;
tmp{4} = resSummaryCombPP.rateMS.msn;

figNameList = {'rateMerge','rateSplit'};

yaxisUnits = '(per frame)';

plotResCombinedDataset(tmp,timeList,condName,condColor,figNameList,dir2save,yaxisUnits,shiftNegTime);

