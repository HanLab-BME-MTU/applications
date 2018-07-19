function [] = pcDump2ImgFolder(expList)
% Outputs images with unique ID
% Assaf Zaritsky, October 2017

close all;

analysisDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/All';
outputDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/cellProfilerInput';

if nargin  == 0
    expList = getExperiments(analysisDir);
end

nExp = length(expList);

for iexp = 1 : nExp
    curExp = expList{iexp};
    expDir = [analysisDir filesep curExp filesep];
    
    filenamesTask = dir(expDir);
    nfilesTask = length(filenamesTask);
    
    for itask = 3 : nfilesTask
        filenameTask = filenamesTask(itask).name;
        
        [pathstr, nameTask, ext] = fileparts(filenameTask);
        taskDir = [expDir filesep nameTask];
        if exist(taskDir,'dir')
            mdFname = [taskDir filesep nameTask '.mat'];
            assert(logical(exist(mdFname,'file')));
            MD =  MovieData.load(mdFname);
            for t = 1 : MD.nFrames_
                I = MD.getChannel(1).loadImage(t);
                imgFname = [outputDir filesep nameTask '_' sprintf('%03d.tif',t)];
                imwrite(I,imgFname);
            end
            clear MD;
        end
    end
end
end

function expList = getExperiments(dataDir)
expList = {};

filenamesExp = dir(dataDir);
nfilesExp = length(filenamesExp);

for iexp = 3 : nfilesExp
    filenameExp = filenamesExp(iexp).name;
    
    [pathstr, nameExp, ext] = fileparts(filenameExp);
    expDir = [dataDir filesep nameExp];
    if exist(expDir,'dir')
        expList{length(expList) + 1} = nameExp;
    end
end
end