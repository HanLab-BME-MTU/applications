function targetStats(mData,exList)



doPlot=1;

% get output directory and excel sheet for writing
homeDir=pwd;
outputDir=uigetdir(pwd,'Select output directory.');
cd(outputDir)
filePath=[];
while exist(filePath,'file')==0
    [fileName,pathName]=uigetfile('*.xls','Select .xls file for data output.');
    filePath=[pathName fileName];
end
cd(homeDir)

% use doPlot variable for output directory path
if doPlot==1
    doPlot=outputDir;
else
    doPlot=[];
end

lastRow=0;

% get list of all target names
targetLabels=getlabels(mData.target);



for iTar=1:length(targetLabels);
    try

        %[movGroup]=splitByDate(mData,iTar,exList);
        [movGroup]=splitByMovie(mData,iTar,exList);

        discrimMatrices=[];
        prop2sample='growthSpeeds';
        [movGroup]=plusTipGetPooledData(mData,movGroup,prop2sample);
        [discrimMatrices]=plusTipSampleData(movGroup,prop2sample,discrimMatrices,doPlot);

        % extract the p-values from the ks test and t-test
        ksTestP=discrimMatrices.growthSpeeds     (find(triu(ones(size(discrimMatrices.growthSpeeds     ,1)),1)));
        tTestP =discrimMatrices.growthSpeedsMeans(find(triu(ones(size(discrimMatrices.growthSpeedsMeans,1)),1)));

        stats(iTar,1).kstest=ksTestP;
        stats(iTar,1).ttest = tTestP;

        %     prop2sample='shrinkSpeeds';
        %     [movGroup]=plusTipGetPooledData(mData,movGroup,prop2sample);
        %     [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);

        save([outputDir filesep 'targetStatsSplitByMovie ' num2str(iTar)],'movGroup','discrimMatrices')

        [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow);
    catch
        disp(['problem with' num2str(iTar)])
    end

end

save([outputDir filesep 'byMoviePvalues'],'stats')

for iTar=3:length(targetLabels);
    try
        [movGroup]=splitByOligo(mData,iTar,exList);

        discrimMatrices=[];
        prop2sample='growthSpeeds';
        [movGroup]=plusTipGetPooledData(mData,movGroup,prop2sample);
        [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);

        %     prop2sample='shrinkSpeeds';
        %     [movGroup]=plusTipGetPooledData(mData,movGroup,prop2sample);
        %     [discrimMatrices]=sampleData(movGroup,prop2sample,discrimMatrices,doPlot);

        % extract the p-values from the ks test and t-test
        ksTestP=discrimMatrices.growthSpeeds     (find(triu(ones(size(discrimMatrices.growthSpeeds     ,1)),1)));
        tTestP =discrimMatrices.growthSpeedsMeans(find(triu(ones(size(discrimMatrices.growthSpeedsMeans,1)),1)));

        stats(iTar,1).kstest=ksTestP;
        stats(iTar,1).ttest = tTestP;

        save([outputDir filesep 'targetStatsSplitByOligo ' num2str(iTar)],'movGroup','discrimMatrices')

        [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow);

    catch
        disp(['problem with' num2str(iTar)])
    end

end

save([outputDir filesep 'byOligoPvalues'],'stats')


%%%%% sub-functions %%%%%
function [lastRow]=writeData2Excel(filePath,discrimMatrices,lastRow)

letters={'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M' 'N'...
    'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'...
    'AA' 'AB' 'AC' 'AD' 'AE' 'AF' 'AG' 'AH' 'AI' 'AJ' 'AK' 'AL' 'AM' 'AN'...
    'AO' 'AP' 'AQ' 'AR' 'AS' 'AT' 'AU' 'AV' 'AW' 'AX' 'AY' 'AZ'...
    'BA' 'BB' 'BC' 'BD' 'BE' 'BF' 'BG' 'BH' 'BI' 'BJ' 'BK' 'BL' 'BM' 'BN'...
    'BO' 'BP' 'BQ' 'BR' 'BS' 'BT' 'BU' 'BV' 'BW' 'BX' 'BY' 'BZ'...
    'CA' 'CB' 'CC' 'CD' 'CE' 'CF' 'CG' 'CH' 'CI' 'CJ' 'CK' 'CL' 'CM' 'CN'...
    'CO' 'CP' 'CQ' 'CR' 'CS' 'CT' 'CU' 'CV' 'CW' 'CX' 'CY' 'CZ'...
    'DA' 'DB' 'DC' 'DD' 'DE' 'DF' 'DG' 'DH' 'DI' 'DJ' 'DK' 'DL' 'DM' 'DN'...
    'DO' 'DP' 'DQ' 'DR' 'DS' 'DT' 'DU' 'DV' 'DW' 'DX' 'DY' 'DZ'};

names = fieldnames(discrimMatrices);
for iName=1:length(names)
    tempMat=discrimMatrices.(names{iName});
    if iscell(tempMat)
        [r c]=size(tempMat);
        range=[letters{2} num2str(lastRow+3) ':' letters{1+c} num2str(lastRow+2+r)];
        xlswrite(filePath,tempMat,range);
        lastRow=lastRow+2+r;
    end

end






function [movGroup]=splitByOligo(mData,iTar,exList)
if nargin<3
    exList=[];
end
movGroup=[];

% get list of all target names
targetLabels=getlabels(mData.target);

% find all the movies corresponding to input target number iTar, and remove
% any that are to be excluded
movNums=find(ismember(mData.target,targetLabels(iTar)))';
movNums=setdiff(movNums,exList)';

% get list of all oligo names corresponding to iTar
iTarOligoLabels=getlabels(nominal(mData.oligo(movNums)));

% loop thru oligo names and find corresponding movie numbers, and remove
% any that are to be excluded
nOligos=length(iTarOligoLabels);
for iOligo=1:nOligos
    movGroup(iOligo,1).common2group=[targetLabels{iTar} '_split_by_oligo'];
    movGroup(iOligo,1).label=[targetLabels{iTar} '_' iTarOligoLabels{iOligo}];
    % movies are split by oligo number only - same oligos on different days
    % will be grouped together
    movNums=find(ismember(mData.oligo,iTarOligoLabels(iOligo)))';
    movGroup(iOligo,1).movNums=setdiff(movNums,exList)';
end


function [movGroup]=splitByMovie(mData,iTar,exList)
if nargin<3
    exList=[];
end
movGroup=[];

% get list of all target names
targetLabels=getlabels(mData.target);

% find all the movies corresponding to input target number iTar, and remove
% any that are to be excluded
movNums=find(ismember(mData.target,targetLabels(iTar)))';
movNums=setdiff(movNums,exList)';

for iMov=1:length(movNums)
    movGroup(iMov,1).common2group=[targetLabels{iTar} '_split_by_movie'];
    movGroup(iMov,1).label=[targetLabels{iTar} '_' mData.oligo{movNums(iMov)}...
        '_' mData.date{movNums(iMov)} '_mov_' num2str(movNums(iMov))];
    movGroup(iMov,1).movNums=movNums(iMov);

end


% this function mainly good for grouping control movies by date
function [movGroup]=splitByDate(mData,iTar,exList)
if nargin<3
    exList=[];
end

% get list of all target names
targetLabels=getlabels(mData.target);

% find all the movies corresponding to input target number iTar, and remove
% any that are to be excluded
movNums=find(ismember(mData.target,targetLabels(iTar)))';
movNums=setdiff(movNums,exList)';

% get list of all date names corresponding to iTar
iTarDateLabels=getlabels(nominal(mData.date(movNums)));

% loop thru date names and find corresponding movie numbers, and remove
% any that are to be excluded
nDates=length(iTarDateLabels);
for iDate=1:nDates
    movGroup(iDate,1).common2group=[targetLabels{iTar} '_split_by_date'];
    movGroup(iDate,1).label=[targetLabels{iTar} '_' iTarDateLabels{iDate}];
    % movies must be split by date and belong to target
    movNums=find(ismember(mData.date,iTarDateLabels(iDate)) & ismember(mData.target,targetLabels(iTar)))';
    movGroup(iDate,1).movNums=setdiff(movNums,exList)';
end
