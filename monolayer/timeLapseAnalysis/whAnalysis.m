function [] = whAnalysis(mainDirname,metaData,timePerFrame,initParams)

[exps] = whArrangeData(mainDirname);
nexp = length(exps);

for i = 1 : nexp
    ind = ismember(metaData.fnames,exps{i}.name);
    pixelSize = metaData.pixelSize{ind};
    if nargin == 4
        whProcessTimeLapse(pixelSize, timePerFrame, mainDirname, exps{i}, initParams);
    else
        whProcessTimeLapse(pixelSize, timePerFrame, mainDirname, exps{i});
    end
end

end