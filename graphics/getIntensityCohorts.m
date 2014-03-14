
% Francois Aguet, 2013

function [res, cT] = getIntensityCohorts(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('lftData', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('CohortBounds_s', [10 20 40 60 80 100 120]);
% ip.addParamValue('ShowVariation', true, @islogical);
% ip.addParamValue('FillMode', 'SEM', @(x) any(strcmpi(x, {'SEM', 'pct'})));
% ip.addParamValue('FrontLayer', false, @islogical);
% ip.addParamValue('ShowBackground', false, @islogical);
ip.addParamValue('Rescale', true, @islogical);
ip.addParamValue('ExcludeVisitors', true, @islogical);
% ip.addParamValue('TrackIndex', []);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('Alpha', 0.05);
ip.addParamValue('SelIndex', [], @iscell); % selection index for each lftData
% ip.addParamValue('RemoveOutliers', false, @islogical);
% ip.addParamValue('AvgFun', @nanmean, @(x) isa(x, 'function_handle'));
ip.addParamValue('LftDataName', 'lifetimeData.mat');
ip.parse(data, varargin{:});
cohortBounds = ip.Results.CohortBounds_s;

lftData = ip.Results.lftData;
if isempty(lftData)
    lftData = getLifetimeData(data,...
        'LifetimeData', ip.Results.LftDataName, 'Scale', ip.Results.Rescale,...
        'Cutoff_f', ip.Results.Cutoff_f, 'ReturnValidOnly', true,...
        'ExcludeVisitors', ip.Results.ExcludeVisitors);
end

framerate = data(1).framerate;

% expand bounds to include longest lifetime
% cohortBounds = [cohortBounds max(vertcat(lftData.lifetime_s))+framerate];

nCh = numel(data(1).channels);
nd = numel(data);
nc = numel(cohortBounds)-1;

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1);

% selIndex = ip.Results.SelIndex;
% if isempty(selIndex)
%     selIndex = arrayfun(@(i) size(i.A,1), lftData, 'unif', 0);
% end

% # data points in cohort (including buffer frames)
b = size(lftData(1).sbA,2);
iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate) + 2*b, 1:nc);

% time vectors for cohorts
cT = arrayfun(@(i) (-b:i-b-1)*framerate, iLength, 'unif', 0);

res(1:nd) = struct('aInterp', [], 'rInterp', []);
cohortBounds(end) = cohortBounds(end)+framerate;
for i = 1:nd
    
    aInterp = NaN(size(lftData(i).A));
    rInterp = NaN(size(lftData(i).A));
    
    nt = numel(lftData(i).lifetime_s);
    
    % cohort index (index of tracks outside of selected cohorts is 0)
    cidx = zeros(nt,1);
    for c = 1:nc
        cidx(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<cohortBounds(c+1)) = c;
    end

    for ch = 1:nCh % channels
        % loop through tracks and interpolate to cohort time vector
        for t = 1:nt
            if cidx(t)~=0%~isnan(cidx(t))
                tidx = 1:lftData(i).trackLengths(t);
                A = [lftData(i).sbA(t,:,ch) lftData(i).A(t,tidx,ch) lftData(i).ebA(t,:,ch)];
                sigma_r = [lftData(i).sbSigma_r(t,:,ch) lftData(i).sigma_r(t,tidx,ch) lftData(i).ebSigma_r(t,:,ch)];
            
                % interpolate to mean length
                np = iLength(cidx(t));
                xi = linspace(1,lftData(i).trackLengths(t) + 2*b, np);
                aInterp(t,1:np,ch) = binterp(A, xi);
                rInterp(t,1:np,ch) = binterp(sigma_r, xi);
            end
        end
    end
    res(i).aInterp = aInterp;
    res(i).rInterp = kLevel*rInterp;
    res(i).cidx = cidx;
end
    
%     for ch = 1:nCh % channels
%         % interpolate tracks to mean cohort length
%         for c = 1:nc % cohorts
%             % tracks in current cohort (above threshold)
%             %cidx = find(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<cohortBounds(c+1));
%             nt = numel(cidx);
%             if nt>0
%                 
%                 % split as a function of slave channel signal
%                 if isfield(lftData(i), 'significantMaster')
%                     sigIdx = lftData(i).significantMaster(:,ch)==1;
%                     res(i).sigIdx{c}(:,ch) = sigIdx(cidx);
%                 else
%                     res(i).sigIdx{c}(:,ch) = ones(numel(cidx),1);
%                 end
%                 if ch==1
%                     res(i).trackIdx{c} = lftData(i).index(cidx);
%                 end
%             else
%                 res(i).interpTracks{ch,c} = NaN(1,iLength(c));
%                 res(i).interpSigLevel{ch,c} = NaN(1,iLength(c));
%                 res(i).sigIdx{c}(:,ch) = NaN;
%                 if ch==1
%                     res(i).trackIdx{c} = NaN;
%                 end
%             end
%         end
%     end

