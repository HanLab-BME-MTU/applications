

% Francois Aguet (last mod. 05/29/2013)

function cmeAnalysis(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('SigmaSource', 'model', @(x) any(strcmpi(x, {'data', 'model'})));
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

% 'RemoveRedundant' inactivated on Windows as a temporary workaround for the problems with
% the KDTree MEX files for windows.
runDetection(data, 'SigmaSource', ip.Results.SigmaSource, 'RemoveRedundant', isunix, opts{:});

settings = loadTrackSettings('Radius', ip.Results.TrackingRadius, 'MaxGapLength', ip.Results.TrackingGapLength);
runTracking(data, settings, opts{:});
runTrackProcessing(data, opts{:});
if numel(data(1).channels)>1
    runSlaveChannelClassification(data, opts{:}, 'np', 5000);
end

lftRes = runLifetimeAnalysis(data, 'RemoveOutliers', true, 'Display', 'off', opts{:});

% Graphical output
plotLifetimes(lftRes, 'DisplayMode', 'print', 'ShowCargoDependent', true);

plotIntensityCohorts(data, 'MaxIntensityThreshold', lftRes.MaxIntensityThreshold,...
    'ShowBackground', false, 'DisplayMode', 'print', 'ScaleSlaveChannel', false,...
    'ShowLegend', false, 'ShowPct', false);
