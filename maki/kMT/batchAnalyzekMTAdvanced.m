% Load movie list
ML = MovieList.load('C:\kjData\maki\newUnarchived\JulieSebastien\1306_cambridgeData\metaMovies\kMTAnalysis\kMTList.mat');

%% Recreate analysis processes

% Flag to reset all analysis
resetAnalysis = true;

if resetAnalysis
    
    arrayfun(@reset,ML); % Wipe out all existing analysis
    
    % Create process to correlate kEB signal with kinetochore dynamics
    ML.addProcess(CorrEBtoKinDynamicsProcess(ML));
        
end

%% Set analysis parameters

%NOTE: space units = pixels
%       time units = frames

% kEB-kDynamics correlation
%general
funParams = ML.getProcess(1).funParams_;
%function-specific
funParams.minDisp = 0.5; %pixels
%general
parseProcessParams(ML.getProcess(1), funParams);

%% Run analysis processes
cellfun(@run,ML.processes_);

%% Set visualization options
%nothing yet, but who knows
