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
% ip.addParamValue('RescalingReference', 'med', @(x) any(strcmpi(x, {'max', 'med'})));
% ip.addParamValue('ScaleSlaveChannel', true, @islogical);
% ip.addParamValue('ScalingFactor', ones(1,nCh));
ip.addParamValue('ExcludeVisitors', true, @islogical);
% ip.addParamValue('TrackIndex', []);
ip.addParamValue('Cutoff_f', 5);
ip.addParamValue('SelIndex', [], @iscell); % selection index for each lftData
% ip.addParamValue('RemoveOutliers', false, @islogical);
% ip.addParamValue('ShowLegend', false, @islogical);
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
cohortBounds = [cohortBounds max(vertcat(lftData.lifetime_s))+framerate];

mCh = 1;%find(strcmp(data(1).source, data(1).channels));
nCh = numel(data(1).channels);
nd = numel(data);
nc = numel(cohortBounds)-1;


selIndex = ip.Results.SelIndex;
if isempty(selIndex)
    selIndex = arrayfun(@(i) size(i.A,1), lftData, 'unif', 0);
end

% # data points in cohort (including buffer frames)
b = size(lftData(1).sbA,2);
iLength = arrayfun(@(c) floor(mean(cohortBounds([c c+1]))/framerate) + 2*b, 1:nc);

% time vectors for cohorts
cT = arrayfun(@(i) (-b:i-b-1)*framerate, iLength, 'unif', 0);

res(1:nd) = struct('aInterp', [], 'sInterp', []);
for i = 1:nd
    
    nt = numel(lftData(i).lifetime_s);
    
    % cohort index (lifetimes > last cohort have index 0)
    cidx = zeros(nt,1);
    for c = 1:nc
        cidx(cohortBounds(c)<=lftData(i).lifetime_s & lftData(i).lifetime_s<cohortBounds(c+1)) = c;
    end

    % loop through tracks and interpolate to cohort time vector
    AInterp = zeros(size(lftData(i).A));
    sInterp = zeros(size(lftData(i).A));
    for ch = 1:nCh
        for t = 1:nt
            % track index
            tidx = 1:lftData(i).trackLengths(t);
            
            A = [lftData(i).sbA(t,:,ch) lftData(i).A(t,tidx,ch) lftData(i).ebA(t,:,ch)];
            bgr = [lftData(i).sbSigma_r(t,:,ch) lftData(i).sigma_r(t,tidx,ch) lftData(i).ebSigma_r(t,:,ch)];
            
            % align to track start
            %w = min(numel(A),iLength);
            %interpTracks(t,1:w) = A(1:w);
            %sigma_r_Ia(t,1:w) = bgr(1:w);
            
            % interpolate to mean length
            tLength = lftData(i).trackLengths(t) + 2*b;
            np = iLength(cidx(t));
            xi = linspace(1,tLength, np);
            AInterp(t,1:np,ch) = binterp(A, xi);
            sInterp(t,1:np,ch) = binterp(bgr, xi);
            
            %interpTracks(t,:) = interp1(1:cLengths(t)+2*b, A, xi, 'cubic');
            %interpTracks(t,:) = binterp(A, xi);
            %sigma_rMat(t,:) = interp1(1:cLengths(t)+2*b, bgr, xi, 'cubic');
            %sigma_rMat(t,:) = binterp(bgr, xi);
        end
    end
    res(i).aInterp = AInterp;
    res(i).sInterp = sInterp;
    res(i).cidx = cidx;
end
