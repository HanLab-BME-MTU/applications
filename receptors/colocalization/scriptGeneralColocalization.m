function [] = scriptGeneralColocalization()
MD = ColorStackMovieData(pwd);
MD.sanityCheck
imageDir = '/project/biophysics/jaqaman_lab/vegf_tsp1/touretLab/CtxB-CD36-Actin/NoTSP';
MD = MovieData.load([imageDir filesep 'colorStackMovieData.mat']);
%MD.setReader(LinearReader);
process = SubResolutionProcess(MD);
MD.addProcess(process);
p = MD.getProcess(1).getParameters();
p.ChannelIndex = [1];
p.detectionParam.psfSigma = 1;
p.detectionParam.calcMethod = 'c';
p.detectionParam.testAlpha = struct('alphaR',0.05,'alphaA',0.1,'alphaD',0.1,'alphaF',0);
p.detectionParam.alphaLocMax = 0.1;
 p.detectionParam.doMMF = 0;
MD.getProcess(1).setParameters(p);
MD.getProcess(1).run;

process = MultiThreshProcess(MD);
MD.addProcess(process);
p = MD.getProcess(2).getParameters();
p.ChannelIndex = 2;
p.GaussFilterSigma = 3;
p.MaxJump = 1;
MD.getProcess(2).setParameters(p);
MD.getProcess(2).run;

process = ColocalizationProcess(MD);
MD.addProcess(process);
p = MD.getProcess(3).getParameters();
p.ChannelRef = 2;
p.ChannelObs = 1;
p.ChannelMask = 2;
p.MethodIndx = 2;
MD.getProcess(3).setParameters(p);
MD.getProcess(3).run;
end