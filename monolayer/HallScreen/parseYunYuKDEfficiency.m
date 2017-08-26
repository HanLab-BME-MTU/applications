
% dbFname - excel file
% kdEfficiencyTH - minimum % if KD
function [] = parseYunYuKDEfficiency(dbFname,kdEfficiencyTH)

% dbFname = '/project/cellbiology/gdanuser/collab/hall/Data20150928/GEFsScreenMetaData20150925.xlsx';
% dbFname = '/project/cellbiology/gdanuser/collab/hall/Data20150928/GEFsScreenMetaData20150925_new.xlsx';
% dbFname = '/project/cellbiology/gdanuser/collab/hall/allData/GEFsScreenMetaData20150424.xlsx';
% dbFname = '/project/cellbiology/gdanuser/collab/hall/allData/GEFsScreenMetaData20150424_Angeles_new.xlsx';
% dbFname = '/project/cellbiology/gdanuser/collab/hall/Data20150424/GEFsScreenMetaData20150424_Angeles_new.xlsx';
% dbFname = '/project/cellbiology/gdanuser/collab/hall/MetaAnalysis/201504/GEFsScreenMetaData20150415.xlsx';
% dbFname = '/project/cellbiology/gdanuser/collab/hall/Data201504/GEFsScreenMetaData20150415_new.xlsx';
% dbFname = 'C:\Users\CellBio1\Dropbox\ResearchInProgress\Hall\MetaData\GEFsScreenMetaData20140911_all.xlsx'
% dbFname = 'C:\Users\CellBio1\Dropbox\ResearchInProgress\Hall\MetaData\GEFsScreenMetaData20141013.xlsx';
% dbFname = 'C:\Users\CellBio1\Dropbox\ResearchInProgress\Hall\MetaData\GEFsScreenMetaData20141201_new.xlsx';
% dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20160201/YunYuValidationData20160202_AZ_0215.xlsx';
% kdEfficiencyTH = 50;
% dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/20151111a/GEFsScreenMetaData20151111_new_noARHGEF40_exclusions_duplicateControlWells.xlsx';
% kdEfficiencyTH = 50;

if nargin == 0
    dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612/ShefaliRevision201612_RhoC.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/Revision201612_RhoC3/ShefaliRevision201612_RhoC_Hairpin3.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data201612_shefali/ShefaliRevision201612All.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20161215_shefali/ShefaliRevision20161215.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd50time200_RHOA/GEFScreenFinal20160526_RHOA.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160523/GEFProjectAll20160523.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/20160523/YunYuFollowup20160523_AZ.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd50time200/GEFScreenFinal20160526.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20160523/YunYuFollowup20160523_AZ.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ProjectAll20160516/GEFProjectAll20160516.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/MetaAnalysis/ScreenFinal/kd0time200/GEFScreenFinal20160513.xlsx';
    %     dbFname = '/project/bioinformatics/Danuser_lab/GEFscreen/analysis/Data20160314/YunYuValidationData20160314.xlsx';
    
        kdEfficiencyTH = 0;
    %     kdEfficiencyTH = 50;
else
    if nargin < 2
        kdEfficiencyTH = -inf;
    end
end

[pathstr,name,ext] = fileparts(dbFname);

firstRow = 2;
% lastRow = 36;

% exps = getExperiments(dbFname,firstRow,lastRow);
[metaData,nExclude] = getExperiments(dbFname,firstRow,kdEfficiencyTH);
% metaData.N = nrows;

metaData.groupsByTreatments = clusterByTreatments(metaData);
metaData.groupsByDays = clusterByDays(metaData);

% outFname = 'C:\Assaf\Yun-Yu\Sos1MetaData20140404.mat';
outFname = [pathstr filesep name '_kd' num2str(kdEfficiencyTH) '.mat'];

save(outFname,'metaData');
end

% function [exps] = getExperiments(fname,firstRow,lastRow)
function [exps, nExclude] = getExperiments(fname,firstRow,kdEfficiencyTH)
exps.treatment = {};
exps.fnames = {};
exps.pixelSize = {};
exps.KD = {};
exps.dates = {};
exps.actualDates = {};

[numbers strings misc] = xlsread(fname);
lastRow = size(strings,1);

[numsAll,txtAll,rawAll] = xlsread(fname,sprintf('A%d:G%d',firstRow,lastRow));  

nExclude = 0;
nUnknownKD = 0;
n = 0;
for line = firstRow : lastRow 
    i = line - firstRow + 1;
    if isnan(numsAll(i,2)) && (length(rawAll{i,5}) == 3)
        %     if isnan(numsAll(i,4)) && strcmp(rawAll{i,5},'N/A')
        KD = -1;
        nUnknownKD = nUnknownKD + 1;        
    else
        KD = numsAll(i,2)*100;% numsAll(i,4)
    end
    
    if KD < kdEfficiencyTH && KD ~= 0 && KD ~= -1
        nExclude = nExclude + 1;
        continue;
    end
    
    n = n + 1;
    splitText = strsplit(txtAll{i,1},'_');
    if length(splitText) == 2
        exps.treatment{n} = [splitText{1} '_{' splitText{2} '}'];
    else if length(splitText) == 1
            exps.treatment{n} = splitText{1};
        else
            error('treatment name should be of the format X_Y, X - gene, Y - sh sequence');
        end
    end
    exps.fnames{n} = txtAll{i,2};    
    exps.pixelSize{n} = numsAll(i,1); %numsAll(i,3);    
    exps.KD{n} = KD;    
    exps.dates{n} = numsAll(i,3); %numsAll(i,5);
    exps.actualDates{n} = rawAll{i,7}; %numsAll(i,4); % 
end
exps.nExclude = nExclude;
exps.nUnknownKD = nUnknownKD;
exps.N = n;

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


% groupsByTreatments - cell array tha holds pairs <treatmentStr,[list of
% indices]>
function [groupsByDays] = clusterByDays(exps)
groupsByDays = {};
for r = 1 : exps.N
    t = 1;
    nGroups = length(groupsByDays);
    curDate = exps.dates(r); curDate = curDate{:}; curDate = num2str(curDate);
    while t <= nGroups
        if strcmp(curDate,groupsByDays{t}.dates) == 1
            groupsByDays{t}.inds = [groupsByDays{t}.inds,r];            
            break;
        else
            t = t + 1;
        end
    end
    if t > nGroups
        groupsByDays{t}.dates = curDate;
        groupsByDays{t}.inds = [r];
    end
end
end
