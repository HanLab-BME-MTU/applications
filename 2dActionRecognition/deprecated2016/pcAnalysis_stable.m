function [] = pcAnalysis_stable(mainDirname,pixelSize,timePerFrame)

[exps] = whArrangeData(mainDirname);
nexp = length(exps);

flags = -1;

for i = 1 : nexp
    pcProcessTimeLapse_stable(pixelSize, timePerFrame, mainDirname, exps{i}.name, flags);
end


end