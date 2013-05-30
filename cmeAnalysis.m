

% Francois Aguet (last mod. 05/29/2013)

function cmeAnalysis(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('TrackingRadius', [3 6], @(x) numel(x)==2);
ip.addParamValue('TrackingGapLength', 2, @(x) numel(x)==1);
ip.parse(varargin{:});
data = ip.Results.data;
if isempty(data)
    data = loadConditionData();
end

opts = {'Overwrite', ip.Results.Overwrite};

runDetection(data, opts{:});

settingsloadTrackSettings('Radius', ip.Results.TrackingRadius, 'MaxGapLength', ip.Results.TrackingGapLength);
runTracking(data, settings, opts{:});
runTrackProcessing(data, opts{:});

runSlaveChannelClassification(data, opts{:});

lftRes = runLifetimeAnalysis(data, 'RemoveOutliers', true, opts{:});

% Additional graphical output
plotIntensityCohorts(data, 'MaxIntensityThreshold', lftRes.T, 'ShowBackground', false,...
    'DisplayMode', 'print', 'ScaleSlaveChannel', false, 'ShowLegend', false, 'ShowPct', false);
