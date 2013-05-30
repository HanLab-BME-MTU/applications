

% Francois Aguet (last mod. 05/29/2013)

function cmeAnalysis(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('TrackingRadius', [3 6], @(x) numel(x)==2);
ip.addParamValue('TrackingGapLength', 2, @(x) numel(x)==1);
ip.addParamValue('Parameters', [], @(x) numel(x)==3);
ip.parse(varargin{:});
data = ip.Results.data;

if isempty(data)
    parameters = ip.Results.Parameters;
    if isempty(parameters)
        parameters = zeros(1,3);
        parameters(1) = input('Enter the N.A. of the objective: ');
        parameters(2) = input('Enter the magnification of the objective: ');
        parameters(3) = 1e-6*input('Enter the camera pixel size, in [um]: ');
    end
    data = loadConditionData('Parameters', parameters);
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
