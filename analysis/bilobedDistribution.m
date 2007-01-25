function [bilobeDataIdlist, bilobeDataSpindleLength, plotData] = bilobedDistribution(correct4Tag,spindleEdges)
%BIOLOBEDDISTRIBUTION collects and plots tag positions along the spindle axis
%
% SYNOPSIS bilobeDataIdlist = bilobedDistribution(kinetochcorrect4TagoreCorrection)
%
% INPUT  correct4Tag : (opt) correction for the distance between the
%                              kinetochore and the tag
%                              {0} No correction
%                               1  Correction by 0.1 um
%                               any other number: correction (in um)
%        spindleEdges (opt) : vector of bin edges for spindle length.
%                            Default : [1,1.2,1.4,1.6,2];
%                            Edge bins exclude extreme data
%
% OUTPUT bilobeDataXX:  structure array.
%                       XX indicates the grouping of the data - either
%                       according to the individual movies, or according to
%                       spindle length. Both arrays have the fields
%                       .spindleLength  ntp-by-1. Spindle length in microns
%                       .cenPositions   ntp-by-2. Relative cen position
%                       .time           ntp-by-1. Timepoints (movie# for
%                                       XX==spindleLength)
%                       .name           char. directory name

%
% MATLAB VERSION (originally written on): 7.1.0.246 (R14) Service Pack 3 Windows_NT
%
%
% USERNAME: Jonas Dorn
% DATE: 14-Jan-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if nargin < 1 || isempty(correct4Tag)
    correct4Tag = 0;
end
switch correct4Tag
    case 0
        tagCorrection = 0;
    case 1
        tagCorrection = 0.1;
    otherwise
        tagCorrection = correct4Tag;
end
if nargin < 2 || isempty(spindleEdges)
    spindleEdges = [1,1.2,1.4,1.6,2];
end

% if true, lower pole (at t=1) will be made spb1
randomAssignment = false;

% if true, plot data for individual distributions
plotIndividual = false;

%========================
% LOAD DATA
%========================

% let the user select top-directory, then let the user choose -data- files
% from listSelectGUI. Load lastResult. If idlist* and if no ? in
% labelcolor, write file into list

% condition: min. 3 tags, no '?'
idlistList = loadIdlistList(cdBiodata(4),...
    'length(idlist(1).stats.labelcolor) > 2 && isempty(strmatch(''?'',idlist(1).stats.labelcolor)) ');

%========================
% CALCULATE PROJECTION
%========================

% loop through idlists, extract positions. Store spindleLength,
% cenPositions, time (for plotting)
% this part of the code is not very clean because it has been copied from
% bilobeProjection.m
nIdlists = length(idlistList);
bilobeDataIdlist(1:nIdlists) = struct('spindleLength',[],...
    'cenPositions',[],'time',[],'name',[]);

for iIdlist = 1:nIdlists

    idlist = idlistList(iIdlist).idlist;

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

    % In case we want to correct for the centromere, we need to have the
    % sxcx vectors normed, too
    [n_s1c1Vector, e_s1c1Vector] = normList(s1c1Vector);
    [n_s2c2Vector, e_s2c2Vector] = normList(s2c2Vector);

    % correct - but don't make the length negative!
    s1c1Vector = s1c1Vector - repmat(min(tagCorrection,n_s1c1Vector-0.08),1,3) .* e_s1c1Vector;
    s2c2Vector = s2c2Vector - repmat(min(tagCorrection,n_s2c2Vector-0.08),1,3) .* e_s2c2Vector;


    % project spb-cen vectors. Distance from spb1
    cen1Dist = dot(s1c1Vector, e_spindleVector, 2);
    cen2Dist = n_spindleVector + dot(s2c2Vector, e_spindleVector, 2);

    % get normalized cen positions.
    cenPosNorm = [cen1Dist./n_spindleVector, cen2Dist./n_spindleVector];
    
    % potentially: make lower spb into spb1
    if randomAssignment
        lowerSpb2 = idlist(goodTime(1)).linklist(spb1idx,11) > ...
            idlist(goodTime(1)).linklist(spb1idx,11);
        if lowerSpb2
            cenPosNorm = cenPosNorm(:,[2,1]);
        end
    end

    % store data
    bilobeDataIdlist(iIdlist).spindleLength = n_spindleVector;
    bilobeDataIdlist(iIdlist).cenPosition = cenPosNorm;
    bilobeDataIdlist(iIdlist).time = find(goodTime);
    bilobeDataIdlist(iIdlist).name = idlistList(iIdlist).name;

    %-------------------
    % individual plots
    %------------------

    %Normalized and absolute positions,
    % smooth histograms of true and 50% flipped positions
    if plotIndividual
        % data plots
        x=find(goodTime);
        figure('Name',idlistList(iIdlist).name)

        % absolute projected position along the spindle
        subplot(2,2,2)
        plot(x,n_spindleVector,'-b.',x,cen1Dist,'-g.',x,cen2Dist,'-r.');
        hold on
        % plot lines for classification
        plot(x([1,end]),[1,1],'-k',x([1,end]),[1.2,1.2],'-k',...
            x([1,end]),[1.4,1.4],'-k',x([1,end]),[1.6,1.6],...
            '-k',x([1,end]),[2,2],'-k')
        ylim([0,2.5])
        title('absolute positions')

        % relative projected position along the spindle
        subplot(2,2,1)
        plot(x,ones(size(x)),'-b',x,cenPosNorm(:,1),'-g.',x,cenPosNorm(:,2),'-r.')
        ylim([0,1])
        title('relative positions')


        % histogram: cumulated positions
        ah = subplot(2,2,3);
        histogram(ah,cenPosNorm);
        xlim(ah,[0,1])
        title('cumulated positions')

        % histogram: 50% flip
        ah = subplot(2,2,4);
        % flip 50% cenPosNorm - it's already stored
        cenPosNorm(2:2:end,:) = 1-cenPosNorm(2:2:end,:);
        histogram(ah,cenPosNorm);
        xlim(ah,[0,1])
        title('cumulated positions - 50% flipped')
    end % if plotIndividual
end % loop

%========================
% BIN SPINDLE LENGTHS
%========================

% collect data
spindleLength = cat(1,bilobeDataIdlist.spindleLength);
allCenPos = cat(1,bilobeDataIdlist.cenPosition);
time = cat(1,bilobeDataIdlist.time);
tIdx = 0;
for i=nIdlists:-1:1
    nTime(i) = length(bilobeDataIdlist(i).time);
    time(tIdx+1:tIdx+nTime) = i;
    tIdx = tIdx + nTime;
end

% make spindle edges
spindleBins = (spindleEdges(2:end)+spindleEdges(1:end-1))/2;
nBins = length(spindleBins);
bilobeDataSpindleLength(1:nBins) = deal(bilobeDataIdlist(1));
for i=1:nBins
    spindleIdx = spindleLength > spindleEdges(i) & ...
        spindleLength < spindleEdges(i+1);
    bilobeDataSpindleLength(i).spindleLength = spindleLength(spindleIdx);
    bilobeDataSpindleLength(i).cenPosition = allCenPos(spindleIdx,:);
    bilobeDataSpindleLength(i).time = time(spindleIdx,:);
    bilobeDataSpindleLength(i).name = ...
        sprintf('Spindle length: %1.2f-%1.2f',...
        spindleEdges(i),spindleEdges(i+1));
end


%========================
% PLOT
%========================

% generate plotData. Make an array with
% - spindle Length
% - bin# along spindle axis
% - weight (1/#of datapoints of each movie)
% - index of movie

% there are two centromere positions per spindleLength
plotData = zeros(2*length(spindleLength),4);
% write spindle lengths, time
plotData(:,[1,4]) = repmat([spindleLength,time],2,1);
% write bins
cenPos = allCenPos(:);
boundaries = (-1/24:1/24:25/24);
for i=1:length(boundaries)-1
    plotData(cenPos>boundaries(i) & cenPos<=boundaries(i+1),2) = i;
end
% write weight
for t=1:nIdlists
    plotData(plotData(:,4)==t,3) = 1/nTime(t);
end

% call bilobe-plotting function
bilobePlot(plotData);





% % plot bilobes with shifting window of 0.2 um
% boundaries = [0.9:0.1:1.9;1.1:0.1:2.1];
% nBoundaries = size(boundaries,2);
% xall = zeros(26,nBoundaries);
% yall=zeros(26,nBoundaries);
% zall=zeros(26,nBoundaries);
% n = zeros(nBoundaries,1);
% yTickLabels = cell(nBoundaries,1);
% % figure, hold on
% for ct = 1:nBoundaries,
%     spindleIdx = (spindleLength>boundaries(1,ct) & spindleLength<boundaries(2,ct));
%     cenPosition = allCenPos(spindleIdx,:);
%     cenPosition(2:2:end,:) = 1-cenPosition(2:2:end,:);
%     [z,x]=histc(cenPosition(:),(-1/24:1/24:25/24));
%     x=(-1/48:1/24:49/48)';
%     z(end)=[];
%     % normalize histograms to equal total number of observations
%     n(ct) = numel(cenPosition);
%     z=z/n(ct);
%     %plot3(x,ct*ones(size(x)),z),
%     xall(:,ct)=x;yall(:,ct)=ct;zall(:,ct)=z;
%     yTickLabels{ct}=sprintf('%1.1f/%i/%i',mean(boundaries(:,ct)),...
%         nnz(spindleIdx),length(unique(time(spindleIdx))));
% end
% figure,surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
% axis tight
% set(gca,'yTickLabel',yTickLabels)
%  view([0 90])
% for ct = 1:nBoundaries,
%     spindleIdx = (spindleLength>boundaries(1,ct) & spindleLength<boundaries(2,ct));
%     cenPosition = allCenPos(spindleIdx,:);
%     cenPosition(2:2:end,:) = 1-cenPosition(2:2:end,:);
%     [z,x]=histc(cenPosition(:),(-1/24:1/24:25/24));
%     x=(-1/48:1/24:49/48)';
%     z(end)=[];
%     % normalize histograms to equal total number of observations
%     n(ct) = numel(cenPosition);
%     z=z/max(z);
%     %plot3(x,ct*ones(size(x)),z),
%     xall(:,ct)=x;yall(:,ct)=ct;zall(:,ct)=z;
%     yTickLabels{ct}=sprintf('%1.1f/%i/%i',mean(boundaries(:,ct)),...
%         nnz(spindleIdx),length(unique(time(spindleIdx))));
% end
% figure,surf(xall,yall,zall,'FaceColor','interp','FaceLighting','phong')
% axis tight
% set(gca,'yTickLabel',yTickLabels)
%  view([0 90])
% 
% 
% % first plot individual: Normalized and absolute positions,
% % smooth histograms of true and 50% flipped positions
% % --> for simplicity moved into loop
% 
% % plot bins
% for i=1:nBins
%     figure('Name',bilobeDataSpindleLength(i).name);
%     cenPos = bilobeDataSpindleLength(i).cenPosition;
% 
%     ah = subplot(1,2,1);
%     histogram(ah,cenPos);
%     xlim(ah,[0,1])
%     title(sprintf('%i normalized cen positions',numel(cenPos)));
%     cenPos(2:2:end,:) = 1-cenPos(2:2:end,:);
%     ah = subplot(1,2,2);
%     histogram(ah,cenPos);
%     xlim(ah,[0,1])
%     title('normalized cen positions - 50% flipped')
% end
% 
% % then plot all (norm, 50% flipped cumulated)
% % potentially: group by spindle length with histogram
% % allCenPos = cat(1,bilobeDataIdlist.cenPosition);
% figure('Name',sprintf('%i movies cumulated',nIdlists))
% ah = subplot(1,2,1);
% histogram(ah,allCenPos);
% xlim(ah,[0,1])
% title('cumulated positions')
% allCenPos(2:2:end,:) = 1-allCenPos(2:2:end,:);
% ah = subplot(1,2,2);
% histogram(ah,allCenPos);
% xlim(ah,[0,1])
% title('cumulated positions - 50% flipped')