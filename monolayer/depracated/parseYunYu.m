function [] = parseYunYu()

dbFname = 'C:\Assaf\Yun-Yu\GTPasesScreenMetaData20140116.xlsx';
firstRow = 2;
lastRow = 18;
nrows = lastRow - firstRow + 1;

exps = getExperiments(dbFname,firstRow,lastRow);
exps.N = nrows;

exps.groupsByTreatments = clusterByTreatments(exps);

outFname = 'C:\Assaf\Yun-Yu\GTPasesScreenMetaData.mat';

save(outFname,'exps');
end

function [exps] = getExperiments(fname,firstRow,lastRow)
exps.treatment = {};
exps.fnames = {};
exps.pixelSize = {};

[numsAll,txtAll,rawAll] = xlsread(fname,sprintf('A%d:D%d',firstRow,lastRow));  

for line = firstRow : lastRow 
    i = line - firstRow + 1;
    exps.treatment{i} = txtAll{i,1};
    exps.fnames{i} = txtAll{i,2};    
    exps.pixelSize{i} = numsAll(i,1);    
end
end

% groupsByTreatments - cell array tha holds pairs <treatmentStr,[list of
% indices]>
function [groupsByTreatments] = clusterByTreatments(exps)
groupsByTreatments = {};
for r = 1 : exps.N
    t = 1;
    nGroups = length(groupsByTreatments);
    while t <= nGroups
        if strcmp(exps.treatment(r),groupsByTreatments{t}.treatment) == 1
            groupsByTreatments{t}.inds = [groupsByTreatments{t}.inds,r];            
            break;
        else
            t = t + 1;
        end
    end
    if t > nGroups
        groupsByTreatments{t}.treatment = exps.treatment(r);
        groupsByTreatments{t}.inds = [r];
    end
end
end

