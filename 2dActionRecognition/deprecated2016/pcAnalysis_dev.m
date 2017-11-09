% Try-catch on whatever data is available
function [] = pcAnalysis_dev(mainDirname,pixelSize,timePerFrame)

[exps] = whArrangeData(mainDirname);
nexp = length(exps);

flags = -1;

for i = 1 : nexp
    pcProcessTimeLapse_dev(pixelSize, timePerFrame, mainDirname, exps{i}.name, flags);
end


end