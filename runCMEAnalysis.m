function runCMEAnalysis(varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addOptional('data', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.parse(varargin{:});
data = ip.Results.data;
if isempty(data)
    data = loadConditionData();
end

runDetection(data, 'Overwrite', ip.Results.Overwrite);
runTracking(data, loadTrackSettings(), 'Overwrite', ip.Results.Overwrite);
runTrackProcessing(data, 'Overwrite', ip.Results.Overwrite);
