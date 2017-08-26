% Current state of the art
function [] = pcAnalysis(mainDirname,pixelSize,timePerFrame)

[exps] = whArrangeData(mainDirname);
nexp = length(exps);

flags = -1;

allParams = cell(1,nexp);
allDirs = cell(1,nexp);
for i = 1 : nexp
     [allParams{i},allDirs{i}] = pcProcessTimeLapse(pixelSize, timePerFrame, mainDirname, exps{i}.name, flags);
end

%% Global analysis
pcCollectGlobalData(exps,allParams,allDirs); % todo:check that this works

end