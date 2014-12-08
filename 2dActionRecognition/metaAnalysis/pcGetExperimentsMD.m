%% Returns list of experimetns and labeling
function [experimentsDirs,labels, labelsStr] = pcGetExperimentsMD(dataDirname, analysisDirname, rule, metaDataFname) % rule is a pointer to a function (e.g., fh = @sin, fh(1))

% [metaData] = pcMetaData();

load(metaDataFname); % metaData

filenames = dir(dataDirname);
nfiles = length(filenames);

iexp = 0;
experimentsDirs = {};
labels = [];
labelsStr = {};

ruleFH = str2func(rule);

for i = 3 : nfiles    
    filename = filenames(i).name;
    [pathstr, name, ext] = fileparts(filename);
    if strcmp(ext, '.nd2')        
        for task = 1 : 20
            nameNew = [name '_s' pad(task,2)];
            [curLabel, curLabelStr] = ruleFH(nameNew,metaData); % nan if not stands in rule, otherwuse a number.
            if isnan(curLabel)
                continue;
            end
            iexp = iexp + 1;
            %             experimentsDirs{iexp}.dir = [analysisDirname name '/' nameNew '/'];
            experimentsDirs{iexp}.dir = [analysisDirname name '/'];
            experimentsDirs{iexp}.name = nameNew;
            labels = [labels curLabel];
            labelsStr{iexp} = curLabelStr;
            assert(length(labelsStr) == length(labels));
        end
    end
end
end