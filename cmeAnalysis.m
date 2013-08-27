%[res, data] = cmeAnalysis(varargin) performs the analysis of clathrin-coated pit dynamics on a set of movies.
% The analysis comprises detection, tracking, and selection of bona fide CCP structures, and generates lifetime
% distribution and intensity cohort plots.
%
% Inputs (optional):
%           data : data structure returned by loadConditionData()
%
% Options:
%            'PlotAll' : true|[false] displays intermediary processing steps
%        'GaussianPSF' : ['model']|'data' toggles between a model-based or data-based estimation of the Gaussian PSF s.d.
%  'TrackingGapLength' : value defines the maximum number of consecutive missed frames in a trajectory. Default: 2
%     'TrackingRadius' : [minRadius maxRadius] search radii for frame-to-frame linking and gap closing. Default: [3 6]
%       'ChannelNames' : cell array of channel names (e.g., {'EGFP-CLCa'}). Default: fluorophore names
%          'Overwrite' : true|[false] overwrites results of previous analyses. 
%
% Outputs:
%            res : analysis results. Lifetime analysis in 'lftRes' field; intensity cohorts in 'cohorts' field
%           data : structure returned by loadConditionData()
%
% The function will ask for acquisition parameters and subsequently for the data location. The following acquisition
% parameters are required for the Gaussian point spread function (PSF) model used for CCP detection:
% numerical aperture (NA) and magnification of the objective used, and the physical pixel size of the camera (in µm).
%
%
% Notes:
%
% Gaussian PSF model parameters:
% ------------------------------
% The algorithm used for CCP detection estimates the Gaussian PSF s.d. based on an accurate model of a TIRFM PSF
% using the objective NA and magnification, camera pixel size, and fluorophores emission wavelength. In the case
% of misestimation of any of these parameters or non-ideal TIR conditions, the resulting PSF model may be sub-optimal
% for detection, and a data-derived parameterization of this model may be preferable.
%
% Tracking parameters:
% --------------------
% The two most important and sensitive parameters are the maximum number of consecutive missed detection, or ?gaps?,
% in a trajectory, and the search radii for frame-to-frame linking of the detections and for gap closing. These
% parameters are sensitive to the imaging frame rate and should be adjusted accordingly.
% The default values are recommended for data acquired at 0.5 frames/sec
%
% Comparing conditions:
% ---------------------
% cmeAnalysis() must be run separately on groups of movies from different experimental conditions
% (i.e., control vs. perturbation), since the automatic thresholds for identifying bona fide CCPs must be determined
% on control data. For such comparisons, first run, i.e.,
%  >> [resCtrl, dataCtrl] = cmeAnalysis;
% followed by
%  >> [resPert, dataPert] = cmeAnalysis('ControlData', resCtrl);
% In the first run, select the parent directory of the control data. In the second run, select the parent directory
% of the perturbation condition.
%
% Parallelization:
% ----------------
% This function takes advantage of Matlab?s parallelization capabilities. To enable this, enter
% >> matlabpool
% in the command prompt.

% Francois Aguet (last mod. 05/29/2013)

function [res, data] = cmeAnalysis(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('GaussianPSF', 'model', @(x) any(strcmpi(x, {'data', 'model'})));
ip.addParamValue('TrackingRadius', [3 6], @(x) numel(x)==2);
ip.addParamValue('TrackingGapLength', 2, @(x) numel(x)==1);
ip.addParamValue('Parameters', [], @(x) numel(x)==3);
ip.addParamValue('ControlData', [], @isstruct);
ip.addParamValue('PlotAll', false, @islogical);
ip.addParamValue('ChannelNames', []);
ip.addParamValue('FirstNFrames', [], @isposint);
ip.addParamValue('DisplayMode', 'screen', @(x) any(strcmpi(x, {'print', 'screen'})));
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

runDetection(data, 'SigmaSource', ip.Results.GaussianPSF, opts{:});

settings = loadTrackSettings('Radius', ip.Results.TrackingRadius, 'MaxGapLength', ip.Results.TrackingGapLength);
runTracking(data, settings, opts{:});
runTrackProcessing(data, opts{:});

if numel(data(1).channels)>1
    runSlaveChannelClassification(data, opts{:}, 'np', 5000);
end

chNames = ip.Results.ChannelNames;
if isempty(chNames)
    chNames = data(1).markers;
end

if ip.Results.PlotAll
    display = 'on';
else
    display = 'off';
end
if isempty(ip.Results.ControlData)
    res.lftRes = runLifetimeAnalysis(data, 'RemoveOutliers', true,...
        'Display', display, opts{:}, 'DisplayMode', ip.Results.DisplayMode, 'SlaveNames', chNames(2:end),...
        'FirstNFrames', ip.Results.FirstNFrames);
else
    res.lftRes = runLifetimeAnalysis(data, 'RemoveOutliers', true,...
        'Display', display, opts{:}, 'DisplayMode', ip.Results.DisplayMode, 'SlaveNames', chNames(2:end),...
        'MaxIntensityThreshold', ip.Results.ControlData.lftRes.MaxIntensityThreshold,...
        'FirstNFrames', ip.Results.FirstNFrames);
end

% Graphical output
%if ~ip.Results.PlotAll % otherwise this is generated in runLifetimeAnalysis()
    plotLifetimes(res.lftRes, 'DisplayMode', ip.Results.DisplayMode, 'PlotAll', false,...
        'SlaveNames', chNames(2:end));
%end

res.cohorts = plotIntensityCohorts(data, 'MaxIntensityThreshold', res.lftRes.MaxIntensityThreshold,...
    'ShowBackground', false, 'DisplayMode', 'print', 'ScaleSlaveChannel', false,...
    'ShowLegend', false, 'ShowPct', false, 'SlaveName', chNames(2:end));
