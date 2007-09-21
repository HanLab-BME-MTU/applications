function [attached, detached] = fuzzyGroupArma(data, options)
%FUZZYGROUPARMA uses fuzzy clustering to group individual trajectories according to ARMA descriptors
%
% SYNOPSIS: fuzzyGroupArma(data,detachedBait, attachedBait, options)
%
% INPUT     data (opt): three additional fields: 
%                   - orderLen OR orderVel (number of fit parameters), 
%                   - type: 'Len' or 'Vel' depending on whether length or
%                     velocities were fitted.
%                   - name (name of the dataSet). 
%                 If empty, the code will ask for a file called
%                 strain_???.mat, an then look for a file called
%                 resFitVelAndLen_???.mat, and, possibly,
%                 lengthSeries_???.mat OR velocitySeries_???.mat
%           options (opt): structure with fields
%             .fishDetached : whether or not to try and remove detached
%                 chromosomes. If 1, the code will use hierarchical
%                 clustering to group the white noise variance of the data.
%                 To make sure that there are detached movies,
%                 fuzzyGroupArma needs a bait (see options.bait).
%             .bait : If you are trying to detect detached chromosomes,
%                 use this field to pass the ARMA descriptors of ndc10-1.
%                 If you don't, the code will ask for the bait if
%                 necessary.
%             .verbose : true if plots should be shown
%             .wnvThreshold : p-value threshold for white noise variance
%             .armaThreshold : p-value threshold for arma parameters
%             .type : 'Len' or 'Vel' - whether to group length or
%                     velocity series. Default: 'Len' 
%
% OUTPUT    attached : n-by-1 structure array with the same fields as data,
%                 characterizing the n significantly different groups. The
%                 name of each group is "attached_#", where # is the group
%                 number.
%           detached : structure with the same fields as data, plus a field
%                 movieList that lists the names of the data sets with
%                 detached chromosomes.
%
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 30-Oct-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================
%% Test Input / Load Data
%============================

% defaults

% fish for detached as default only if manually loading data
if nargin == 0 || isempty(data)
    def_options.fishDetached = true;
else
    def_options.fishDetached = false;
end
% default bait is empty, of course
def_options.bait = [];
% plot is good
def_options.verbose = true;
def_options.wnvThreshold = 1e-12;
def_options.armaThreshold = 5e-5;
def_options.recalcArma = 0; % 1: recalc till threshold. 2: also, try all combinations
def_options.recalcWnv = 0;
def_options.recalcCenters = 1; % can't have 0 here!
def_options.type = 'Len';

% set options now, because we might need them for loading data
if nargin < 2 || isempty(options)
    options = struct;
end
% assign default options
defaultOptions = fieldnames(def_options);
for fn = 1:length(defaultOptions)
    if ~isfield(options,defaultOptions{fn}) || isempty(options.(defaultOptions{fn}))
        options.(defaultOptions{fn}) = def_options.(defaultOptions{fn});
    end
end

% ensure correct spelling
options.type = lower(options.type);
options.type(1) = upper(options.type(1));


% % check verbose
% if nargin < 3 || isempty(verbose)
%     verbose = 1;
% end

if nargin == 0 || isempty(data)
    data = groupArma_loadData(options);
end % load data

% find common name from first data entry
out=regexp(data(1).name,'([\w]+)_movie','tokens');
commonName = out{1};
commonName = commonName{1};

% check for bait
if options.fishDetached
    if isempty(options.bait)
        h = helpdlg('Please select one bait strain');
        uiwait(h);

        bait = groupArma_loadData(options);
    else
        bait = options.bait;
    end
    if length(bait) > 1
        error('only bait of length 1 allowed')
    else
        data = [data;bait];
    end
end

%get number of data sets
nData = length(data);

% initialize output
dataFN = fieldnames(data);
attached = struct;
detached = struct;
for i=1:length(dataFN)
    attached.(dataFN{i}) = [];
    detached.(dataFN{i}) = [];
end



%============================




%============================
% CLUSTER DATA
%============================

% run groupArmaDescriptors. If we fish for detached chromosomes, allow
% WN-clustering, otherwise just run ARMA parameters.
% Always recalc until threshold is reached.

gad_options.wnv1_cutoff = -1 - options.fishDetached;
gad_options.wnv1_mode = [options.recalcWnv options.wnvThreshold];
gad_options.arma_cutoff = options.armaThreshold*5; % cutoff for starting centers is at higher p-values
gad_options.arma_mode = [options.recalcArma,options.armaThreshold];
gad_options.wnv2_cutoff = -1;
gad_options.wnv2_mode = [options.recalcWnv options.wnvThreshold];

[groupIdx,groupedData] = ...
    groupArmaDescriptors(data,gad_options,options.verbose);

% check clustering output:
% - remove detached chromosomes
% - find the number of potential groups
% - initialize fuzzy clustering

% remove detached chromosomes
if options.fishDetached
    % check in groupIdx for everything that shares a group with the bait
    baitIdx = groupIdx(end,1);
    detachedList = find(groupIdx(:,1) == baitIdx);
    attachedList = find(groupIdx(:,1) ~= baitIdx);

    % from the detachedList: calculate ARMA descriptors in place of the
    % bait
    % detachedList(end) is the bait
    %     data(detachedList(end)).lengthSeries = [];
    %     data = ...
    %         groupArma_recalcArma(data,...
    %         {detachedList(1:end-1),detachedList(end)},2,[],...
    %         sprintf('%s_detached',commonName));
    %     % store detached characteristics
    %     detached = data(end);
    %     detached.data = data(detachedList(1:end-1));
    %     detached.movieList = {data(detachedList(1:end-1)).name}';
    %     detached.movieIdxList = detachedList(1:end-1);

    % run fuzzy clustering on the detached data
    detOpt = options;
    detOpt.fishDetached = 0;
    detached = fuzzyGroupArma(data(detachedList(1:end-1)),detOpt);



    % remove detached from groupedData
    id1 = groupedData.collectedData(1).dataIdx;
    if any(id1 == nData)
        % bait is in first group - keep first
        groupedData = groupedData.collectedData(2).subGroups;
    else
        % bait is in second group
        groupedData = groupedData.collectedData(1).subGroups;
    end

    % remove bait entry
    data(detachedList) = [];
    nData = length(data);
else
    %detachedList = [];
    attachedList = (1:nData)';
    % remove first grouping level
    groupedData = groupedData.collectedData(1).subGroups;
end

% set attached
attached = data(1);
attached.data = data;
attached.movieList = {data.name}';
attached.movieIdxList = attachedList;

% check cluster tendency of attached data
fuzzyGroupArma_visualize(data);



%-------- FUZZY CLUSTERING --------------

% read number of groups, dataIndices from groupedData
% if nCenters == 1, find the farthest two insignificantly different groups

% Future: check groups for size - if X movies or less: remove

nCenters = length(groupedData.collectedData);

% collect indices
if nCenters == 1
    % find two most different groups
    clusterIdx = cluster(groupedData.collectedData(1).links,'MaxClust',2);
    for iCenter = 2:-1:1
        groupIdx = find(clusterIdx==iCenter);
        tmp = ...
            groupArma_recalcArma(groupedData.collectedData.data(groupIdx),...
            {1:length(groupIdx),...
            length(groupIdx)+1},options.recalcCenters,[],...
            sprintf('%s_attached_%i',commonName,iCenter));
        centers(iCenter) = tmp(end);
    end
else
    % calculate cluster prototype for every group. A potential center needs
    % at least 2 data points
    for iCenter = nCenters:-1:1
        nData = length(groupedData.collectedData(iCenter).dataIdx);
        if nData > 1
            tmp = ...
                groupArma_recalcArma(groupedData.collectedData(iCenter).data,...
                {1:nData,nData+1},...
                options.recalcCenters,[],...
                sprintf('%s_attached_%i',commonName,iCenter));
            centers(iCenter) = tmp(end);
        else
            centers(iCenter) = [];
        end
    end
end

% recount centers
nCenters = length(centers);

% loop until we only have significantly different centers. Even if we only
% find one cluster, we want to recalculate to remove outliers.
done = false;
while ~done

    % do fuzzy clustering
    eta = log10(options.armaThreshold); % negative eta -> no recalc
    m = 1.1;
    [centers, membership, distances, eta] = ...
        possibilisticClustering(data,centers,...
        'fuzzyGroupArma_centerFunction','fuzzyGroupArma_distanceFunction',...
        eta,m, ...
        struct('names',{{centers.name}},'recalc',options.recalcCenters), []);

    

        if nCenters > 1
            % for all the centers: p-values
            compMat = repmat(NaN,nCenters,nCenters);
            for iCenter=2:nCenters
                for jCenter=1:iCenter-1
                    compMat(iCenter,jCenter) = ...
                        armaxModelComp(centers(iCenter),centers(jCenter));
                end
            end

            % find the maximum p-value (=minimum -log10(p))
            [minLogP, minIdx] = nanmin(compMat(:));
        else
            minLogP = inf;
            removeCenters = false;
        end

        % check whether maxP is significant
        if minLogP > -log10(options.armaThreshold)
            % minLogP is above the threshold - all distances are significant
            done = true;

            % fill in attached
            for iCenter = nCenters:-1:1
                for fn=fieldnames(centers(iCenter))'
                    attached(iCenter).(fn{1}) = centers(iCenter).(fn{1});
                end
            end

            % fill in some more for attached(1)
            attached(1).membership = membership;
            attached(1).distances = distances;
            attached(1).eta = eta;

        else
            
            removeCenters = true;
            % loop to remove centers with p-values above 90%
            while removeCenters


                % find the two groups that should be fused
                [ig,jg] = ind2sub(size(compMat),minIdx);

                newMembership = membership(:,ig) + membership(:,jg);
                centers(min(ig,jg)) = fuzzyGroupArma_centerFunction(data,newMembership,[],...
                    struct('recalc',1,'names',{{sprintf('%s_attached_tmp',commonName)}}));

                % remove center
                centers(max(ig,jg)) = [];
                nCenters = nCenters - 1;

                % rename centers
                for iCenter = 1:nCenters
                    centers(iCenter).name = ...
                        sprintf('%s_attached_%i',commonName,iCenter);
                end

                % now check whether we want to remove additional centers
                if nCenters < 3
                    % always try to cluster again with two centers
                    % maybe try new starting values
                    removeCenters = false;
                else
                    % look at compMat
                    compMat = repmat(NaN,nCenters,nCenters);
                    for iCenter=2:nCenters
                        for jCenter=1:iCenter-1
                            compMat(iCenter,jCenter) = ...
                                armaxModelComp(centers(iCenter),centers(jCenter));
                        end
                    end
                    % find the maximum p-value (=minimum -log10(p))
                    [minLogP, minIdx] = nanmin(compMat(:));

                    if minLogP > -log10(0.9)
                        % the centers are relatively different. Stop removing
                        removeCenters = false;
                    else
                        % remove more
                    end
                end

            end % remove centers
        end % if minLogP


end


