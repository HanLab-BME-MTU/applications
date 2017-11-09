%% Returns list of experimetns and labeling
function [experimentsDirs,labels, labelsStr] = pcGetExperiments(baseDirname, rule, metaDataFname) % rule is a pointer to a function (e.g., fh = @sin, fh(1))

% [metaData] = pcMetaData();

load(metaDataFname); % metaData

dirnames = dir(baseDirname);
ndirs = length(dirnames);

iexp = 0;
experimentsDirs = {};
labels = [];
labelsStr = {};

for d = 3 : ndirs
    dayDirname = [baseDirname dirnames(d).name];
    
    if ~exist(dayDirname,'dir')
        continue;
    end
    
    filenames = dir(dayDirname);
    nfiles = length(filenames);
    
    for i = 1 : nfiles
        filename = filenames(i).name;
        [pathstr, name, ext] = fileparts(filename);
        if strcmp(ext, '.tif')
            ruleFH = str2func(rule);
            [curLabel, curLabelStr] = ruleFH(name,metaData); % nan if not stands in rule, otherwuse a number.
            if isnan(curLabel)
                continue;
            end
            iexp = iexp + 1;            
            experimentsDirs{iexp}.dir = dayDirname;
            experimentsDirs{iexp}.name = name;
            labels = [labels curLabel];
            labelsStr{iexp} = curLabelStr;
            assert(length(labelsStr) == length(labels));
        end
    end
end
end