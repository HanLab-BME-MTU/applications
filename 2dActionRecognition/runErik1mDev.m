% No motion estimation / segmentation, ignore missing data

function [] = runErik1mDev()
pixelSize = 0.325;
timePerFrame = 1;
mainDirname = '/work/gdanuser/azaritsky/UTSW/Data/Erik/POC_1min/';
analysis(mainDirname,pixelSize,timePerFrame)
end



function [] = analysis(mainDirname,pixelSize,timePerFrame)

dirnames = dir(mainDirname);
ndirs = length(dirnames);

for d = 3 : ndirs
    name = dirnames(d).name;
    dirname = [mainDirname name '/'];
    try
        pcAnalysis_dev(dirname,pixelSize,timePerFrame);
    catch e
    end
    if exist([dirname '1/'],'dir')
        try
            pcAnalysis_dev([dirname '1/'],pixelSize,timePerFrame);
        catch e
        end
    end
    if exist([dirname '2/'],'dir')
        try
            pcAnalysis_dev([dirname '2/'],pixelSize,timePerFrame);
        catch e
        end
    end
end
end