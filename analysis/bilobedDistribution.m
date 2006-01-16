function bilobeData = bilobedDistribution
%BIOLOBEDDISTRIBUTION collects and plots tag positions along the spindle axis
%
% SYNOPSIS bilobeData = bilobedDistribution
%
% OUTPUT bilobeData:  structure array. Contains for every idlist used:
%                       .spindleLength  ntp-by-1. Spindle length in microns
%                       .cenPositions   ntp-by-2. Relative cen position
%                       .time           ntp-by-1. Timepoints
%                       .name           char. directory name
%
%
%
% MATLAB VERSION (originally written on): 7.1.0.246 (R14) Service Pack 3 Windows_NT
%
%
%
% USERNAME: Jonas Dorn
% DATE: 14-Jan-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%========================
% LOAD DATA
%========================

% let the user select top-directory, then let the user choose -data- files
% from listSelectGUI. Load lastResult. If idlist* and if no ? in
% labelcolor, write file into list

% select top-directory
topDir = uigetdir(cdBiodata(4));

% find all -data- files
fileList = searchFiles('-data-','log',topDir,1);

% launch listSelectGUI
selectIdx = listSelectGUI(fileList(:,1),[],'move');
% shorten fileList
fileList = fileList(selectIdx,:);

% loop through selected files
nFiles = length(fileList);
idlistList = cell(nFiles,2); % idlist, dirName
idlistCt = 1;
for iFile = 1:nFiles
    dataFileName = fullfile(fileList{iFile,2},fileList{iFile,1});
    % load lastResult
    load(dataFileName,'lastResult');
    % check if idlist
    if strfind(lastResult,'idlist')
        % if yes: load it
        tmp=load(dataFileName,lastResult);
        idlist = tmp.(lastResult);
        % check labelcolor: min. 3 tags? No '?'?
        if length(idlist(1).stats.labelcolor) > 3 && isempty(strmatch('?',idlist(1).stats.labelcolor))
            % store idlist
            idlistList{idlistCt,1} = idlist;
            % store directoryName
            dirName = fileList{iFile,2};
            fileSepIdx = strfind(dirName,filesep);
            idlistList{idlistCt,2} = dirName(fileSepIdx(end)+1:end);

            % update counter
            idlistCt = idlistCt + 1;
        end
    end
end
% remove empty entries
idlistList(idlistCt:end,:) = [];

%========================
% CALCULATE PROJECTION
%========================

% loop through idlists, extract positions. Store spindleLength,
% cenPositions, time (for plotting)
% this part of the code is not very clean because it has been copied from
% bilobeProjection.m
nIdlists = idlistCt - 1;
bilobeData(1:nIdlists) = struct('spindleLength',[],...
    'cenPositions',[],'time',[],'name',[]);

for iIdlist = 1:nIdlists

    idlist = idlistList{iIdlist,1};

    % find indices of spb, cen
    spb1idx = find(strcmpi(idlist(1).stats.labelcolor,'spb1'));
    spb2idx = find(strcmpi(idlist(1).stats.labelcolor,'spb2'));
    cen1idx = find(strcmpi(idlist(1).stats.labelcolor,'cen1'));
    cen2idx = find(strcmpi(idlist(1).stats.labelcolor,'cen2'));
    if isempty(cen2idx)
        cen2idx = cen1idx;
    end



    % loop through idlist, find spb-axis, spb-cen vectors
    tmax = length(idlist);
    spindleVector = zeros(tmax,3);
    s1c1Vector = zeros(tmax,3);
    s2c2Vector = zeros(tmax,3);
    s1s2c1c2int = zeros(tmax,4);
    goodTime = zeros(tmax,1);

    for t = 1:tmax
        if ~isempty(idlist(t).linklist)

            % find vector in direction of spindle
            spindleVector(t,:) =...
                diff(idlist(t).linklist([spb1idx, spb2idx],9:11));

            % spb - cen vectors
            s1c1Vector(t,:) = ...
                diff(idlist(t).linklist([spb1idx, cen1idx],9:11));
            s2c2Vector(t,:) = ...
                diff(idlist(t).linklist([spb2idx, cen2idx],9:11));

            % store intensities
            s1s2c1c2int(t,:) = ...
                idlist(t).linklist([spb1idx, spb2idx, cen1idx, cen2idx], 8)';

            % remember goodTime
            goodTime(t) = 1;

        else
            % do nothing
        end
    end % for t = 1:tmax

    % shorten Vectors
    spindleVector(~goodTime,:) = [];
    s1c1Vector(~goodTime,:) = [];
    s2c2Vector(~goodTime,:) = [];
    s1s2c1c2int(~goodTime,:) = [];

    % normalize spindleVector
    [n_spindleVector, e_spindleVector] = normList(spindleVector);

    % project spb-cen vectors. Distance from spb1
    cen1Dist = dot(s1c1Vector, e_spindleVector, 2);
    cen2Dist = n_spindleVector + dot(s2c2Vector, e_spindleVector, 2);

    % get normalized cen positions.
    cenPosNorm = [cen1Dist./n_spindleVector, cen2Dist./n_spindleVector];

    % store data
    bilobeData(iIdlist).spindleLength = n_spindleVector;
    bilobeData(iIdlist).cenPosition = cenPosNorm;
    bilobeData(iIdlist).time = find(goodTime);
    bilobeData(iIdlist).name = idlistList{iIdlist,2};

    %-------------------
    % individual plots
    %------------------

    %Normalized and absolute positions,
    % smooth histograms of true and 50% flipped positions

    % data plots
    x=find(goodTime);
    figure('Name',idlistList{iIdlist,2})

    % absolute projected position along the spindle
    subplot(2,2,2)
    plot(x,n_spindleVector,'-b.',x,cen1Dist,'-g.',x,cen2Dist,'-r.');
    title('absolute positions')

    % relative projected position along the spindle
    subplot(2,2,1)
    plot(x,ones(size(x)),'-b',x,cenPosNorm(:,1),'-g.',x,cenPosNorm(:,2),'-r.')
    title('relative positions')

    % histogram: cumulated positions
    ah = subplot(2,2,3);
    histogram(cenPosNorm,[],ah);
    xlim(ah,[0,1])
    title('cumulated positions')

    % histogram: 50% flip
    ah = subplot(2,2,4);
    % flip 50% cenPosNorm - it's already stored
    cenPosNorm(2:2:end,:) = 1-cenPosNorm(2:2:end,:);
    histogram(cenPosNorm,[],ah);
    xlim(ah,[0,1])
    title('cumulated positions - 50% flipped')

end % loop



%========================
% PLOT
%========================

% first plot individual: Normalized and absolute positions,
% smooth histograms of true and 50% flipped positions
% --> for simplicity moved into loop

% then plot all (norm, 50% flipped cumulated)
% potentially: group by spindle length with histogram
allCenPos = cat(1,bilobeData.cenPosition);
figure('Name',sprintf('%i movies cumulated',nIdlists))
ah = subplot(1,2,1);
histogram(allCenPos,[],ah);
xlim(ah,[0,1])
title('cumulated positions')
allCenPos(2:2:end,:) = 1-allCenPos(2:2:end,:);
ah = subplot(1,2,2);
histogram(allCenPos,[],ah);
xlim(ah,[0,1])
title('cumulated positions - 50% flipped')