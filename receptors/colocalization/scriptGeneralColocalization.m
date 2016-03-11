function [] = scriptGeneralColocalization()
MD = ColorStackMovieData(pwd);
MD.sanityCheck
imageDir = '/project/biophysics/jaqaman_lab/vegf_tsp1/touretLab/CtxB-CD36-Actin/NoTSP';
MD = MovieData.load([imageDir filesep 'colorStackMovieData.mat']);
%MD.setReader(LinearReader);
process = SubResolutionProcess(MD);
MD.addProcess(process);
process = MultiThreshProcess(MD);
MD.addProcess(process);
process = ColocalizationProcess(MD);
MD.addProcess(process);

for j =1:4
    load(strcat(condition{j},'/colorStackMovieData.mat'));
p = MD.getProcess(1).getParameters();
p.ChannelIndex = [1];
p.detectionParam.psfSigma = 0.65;
p.detectionParam.calcMethod = 'g';
p.detectionParam.testAlpha = struct('alphaR',0.1,'alphaA',0.1,'alphaD',0.1,'alphaF',0);
p.detectionParam.alphaLocMax = 0.15;
 p.detectionParam.doMMF = 0;
MD.getProcess(1).setParameters(p);
MD.getProcess(1).run;
MD.sanityCheck
end
% process = MultiThreshProcess(MD);
% MD.addProcess(process);
p = MD.getProcess(2).getParameters();
p.ChannelIndex = 1;
p.GaussFilterSigma = 3;
p.MaxJump = 1;
MD.getProcess(2).setParameters(p);
MD.getProcess(2).run;

% process = ColocalizationProcess(MD);
% MD.addProcess(process);
for k = 1:8
    load(strcat(condition{k},'/colorStackMovieData.mat'))
    p = MD.getProcess(3).getParameters();
    p.ChannelRef = 2;
    p.ChannelObs = 3;
    p.ChannelMask = 3;
    p.SearchRadius = 2;
    MD.getProcess(3).setParameters(p);
    MD.getProcess(3).run;
end
end