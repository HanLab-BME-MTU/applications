function data = groupArma_loadData(options)
%GROUPARMA_LOADDATA is a loader routine for the grouping of ARMA descriptors
%
% options: need fields
%   multiply
%   type
% see groupArmaDescriptors for more info

if isunix
    oldDiag = getappdata(0,'UseNativeSystemDialogs');
    setappdata(0,'UseNativeSystemDialogs',false)
    [strainFile, dataPath] = uigetfile('strains_*.mat','Load strain list(s)!','MultiSelect','on');
    setappdata(0,'UseNativeSystemDialogs',oldDiag)
else
    [strainFile, dataPath] = uigetfile('strains_*.mat','Load strain list(s)!','MultiSelect','on');
end


if ~iscell(strainFile) && any(strainFile == 0)
    disp('--groupArmaDescriptors aborted')
    return
end

% make strainFile into a cell, so that we can use the code for
% multiSelections with one selection.
if ~iscell(strainFile)
    strainFile = {strainFile};
end

% read strainInfo so that user may select the files to be analyzed
% directly without having to wait for the data to be read from the
% server


% load data in loop. This is, unfortunately, sensitive to changes in the
% list and ordering of fields. One possible solution would be to
% define the key list of fields above, and to only read those

% load first strainInfo, then add the rest to the list, so that we
% automatically get the proper list of fieldnames

strainInfo = load(fullfile(dataPath,strainFile{1}));
% strainInfo is a structure in itself
fn = fieldnames(strainInfo);
strainInfo = strainInfo.(fn{1});

% loop for all the other files. If one file, the loop isn't executed
for iFile = 2:length(strainFile)
    strainInfoTmp = load(fullfile(dataPath,strainFile{iFile}));
    fnTmp = fieldnames(strainInfoTmp);
    strainInfoTmp = strainInfoTmp.(fnTmp{1});
    % reorder fields, so that we're at least insensitive to that
    strainInfoTmp = orderfields(strainInfoTmp,strainInfo);
    strainInfo = [strainInfo(:);strainInfoTmp(:)];
end





% let the user select the strains
[nameList{1:length(strainInfo)}] = deal(strainInfo.name);
selectionIdx = listSelectGUI(nameList);

if isempty(selectionIdx)
    disp('--groupArmaCoefficients aborted')
    return
end




% read armaData, lengthData. StrainInfo is a structure with
% the dataFile-names, while armaData and lengthData are collections of
% files that are labeled according to the dataFile-names.

% loop through potential list of files
armaData = struct;
lengthData = struct;

for iFile = 1:length(strainFile)

    dataSuffix = regexp(strainFile{iFile},'strains_([\w]+).mat','tokens');
    dataSuffix = char(dataSuffix{1});
    % armaData is a collection of many files, while strainInfo is already
    % the correct structure that contains the optimal model for each
    armaDataTmp = load(fullfile(dataPath,sprintf('resFitVelAndLen_%s',dataSuffix)));

    % this would be really annoying without a loop
    fn = fieldnames(armaDataTmp);
    for i=1:length(fn)
        armaData.(fn{i}) = armaDataTmp.(fn{i});
    end

    % try load length file (or velocity)
    % it will be called lengthData temporarily, but will be assigned
    % properly as length or vel in the data structure
    try
        switch options.type
            case 'Len'
        lengthDataTmp = load(fullfile(dataPath,sprintf('lengthSeries_%s',dataSuffix)));
            case 'Vel'
                lengthDataTmp = load(fullfile(dataPath,sprintf('velocitySeries_%s',dataSuffix)));
        end % switch
    catch
        lengthDataTmp = struct;
    end
    fn = fieldnames(lengthDataTmp);
    for i=1:length(fn)
        lengthData.(fn{i}) = lengthDataTmp.(fn{i});
    end

end



% get number of data
nData = length(selectionIdx);
orderType = ['order',options.type];

% read the correct model. If orderLen is [], skip.
for i=nData:-1:1,
    % construct fieldname - length or velocity?

    if isempty(strainInfo(selectionIdx(i)).(orderType))
        % remove only if there is data already
        if exist('data','var')
            data(i) = [];
            goodIdx(i) = [];
        end
    else
        % reorder data fields, in case there is a change in field-order in the future
        %             if exist('data','var')
        %             data = orderfields(armaData.(sprintf('fit%s%s',options.type,strainInfo(selectionIdx(i)).name))...
        %             (strainInfo(selectionIdx(i)).orderLen(1)+1,strainInfo(selectionIdx(i)).orderLen(2)+1),data);
        %             end

        % read the "correct" ARMA model according to orderVel/orderLen
        data(i)=...
            armaData.(sprintf('fit%s%s',options.type,strainInfo(selectionIdx(i)).name))...
            (strainInfo(selectionIdx(i)).(orderType)(1)+1,strainInfo(selectionIdx(i)).(orderType)(2)+1);
        goodIdx(i) = true;
    end
end
data = data(:);


% read strainInfo into data
[data.name] = deal(strainInfo(selectionIdx(goodIdx)).name);
% store type in data
[data.type] = deal(options.type);
[data.(orderType)] = deal(strainInfo(selectionIdx(goodIdx)).(orderType));

% read length into data (if there is any). We can't do this above,
% because filling up the structure requires the fields to be the same
% read either lengthSeries or velSeries (to remain consistent with Khuloud)
switch options.type
    case 'Len'
        seriesType = 'lenSeries';
        longType = 'length';
    case 'Vel'
        seriesType = 'velSeries';
        longType = 'vel';
    otherwise
        error('Type %s not supported yet',options.type)
end

if ~isempty(lengthData)
    % we already know what part of the selection is good
    goodIdx = find(goodIdx);
    for i=1:length(goodIdx)
        data(i).(seriesType) = ...
            lengthData.(sprintf('%s%s',longType,strainInfo(selectionIdx(goodIdx(i))).name));
        if options.multiply
            data(i).(seriesType) = repmat(data(i).lengthSeries,[1,options.multiply]);
            data(i).numObserve = data(i).numObserve * options.multiply;
        end
    end
end