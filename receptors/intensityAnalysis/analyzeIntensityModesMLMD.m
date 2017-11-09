function analyzeIntensityModesMLMD(MLMD,alpha,variableMean,variableStd,numModeMinMax,...
    plotResults,logData,modeParamIn,ampOrInt,ratioTol)
%analyzeIntensityModesMLMD does a modal analysis of particle intensity distributions in a list of movies
%
%SYNOPSIS  analyzeIntensityModesMLMD(MLMD,alpha,variableMean,variableStd,numModeMinMax,...
%    plotResults,logData,modeParamIn,ampOrInt,ratioTol)
%
%INPUT 
%   Mandatory
%       MLMD: MovieList or MovieData object for movie(s) to be analyzed
%   Optional
%       alpha        : Alpha-value for the statistical test that compares the
%                      fit of n+1 Gaussians to the fit of n Gaussians.
%                      Default: 0.05.
%       variableMean : Flag with multiple values:
%                      - 0 if assuming the fixed relationship
%                      (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                      - 1 if there is no relationship between the means of
%                      different Gaussians.
%                      - m > 1 if assuming the same fixed relationship as 0
%                      but that the first detected Gaussian is actually the
%                      mth Gaussian in the relationship.
%                      Optional. Default: 0.
%       variableStd  : Flag with multiple values:
%                      - 0 if assuming that all Gaussians have the same
%                      standard deviation. 
%                      - 1 if there is no relationship
%                      between the standard deviations of different
%                      Gaussians.
%                      - 2 if assuming the relationship
%                      (std of nth Gaussian) = sqrt(n) * (std of 1st Gaussian).
%                      This relationship is generalized if variableMean > 1.
%                      variableStd can equal 2 only if variableMean is not
%                      1.
%                      Optional. Default: 0.
%       numModeMinMax: Vector with minimum and maximum number of modes
%                      (Gaussian or log-normal) to fit. 
%                      Default: [1 9].
%                      If only one value is input, it will be taken as the
%                      maximum.
%       plotResults  : 1 if results are to be plotted, 0 otherwise.
%                      Default: 1.
%       logData      : 1 to fit intensity distributions with log-normal, 0
%                      to fit with normal.
%                      Default: 0.
%       modeParamIn  : Matrix with number of rows equal to number of
%                      modes and two columns indicating the mean/M
%                      (Gaussian/lognormal) and standard deviation/S
%                      (Gaussian/lognormal) of each mode. If input, the
%                      specified mode parameters are used, and only the mode
%                      amplitudes are determined by data fitting. In this
%                      case, the input alpha, variableMean, variableStd
%                      and numModeMinMax are not used.
%                      Default: [].
%       ampOrInt     : 1 to analyze Gaussian amplitudes, 2 to analyze raw
%                      intensities. If 2, the field intRawMinusBg must be
%                      present.
%                      Default: 1.
%       ratioTol     : Tolerance for ratio between mean/std of 1st Gaussian
%                      and mean/std of subsequent Gaussians.
%                      If 0, ratio is taken strictly.
%                      If > 0, ratio is allowed to wiggle by ratioTol about
%                      the theoretial ratio.
%                      If one entry, same value is used for both mean and
%                      std. If two entries, then first is for mean and
%                      second is for std.
%                      Example: If ratioTol = 0.1, then mean of 2nd Gaussian
%                      can vary between 1.9 and 2.1 of mean of 1st Gaussian
%                      (instead of exactly 2). Same holds for std.
%                      If ratioTol = [0.1 0.2], then mean can vary by 0.1
%                      and std by 0.2.
%                      Default: 0.
%                      Option currently implemented only for 3 cases:
%                      variableMean ~= 1 and variableStd = 1, 2 or 3.
%
%OUPUT Output is saved in directory IntModalAnalysis belonging to each
%      analyzed movie. See analyzeIntensityModes for detailed description
%      of saved output.
%
%Khuloud Jaqaman, March 2015

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--analyzeIntensityModesMLMD: Function needs at least 1 input argument!');
    return
end

%assign default values of optional input variables
alpha_def = 0.05;
variableMean_def = 0;
variableStd_def = 0;
numModeMinMax_def = [1 9];
plotResults_def = 1;
logData_def = 0;
modeParamIn_def = [];
ampOrInt_def = 1;
ratioTol_def = 0;

%check alpha
if nargin < 2 || isempty(alpha)
    alpha = alpha_def;
end

%check variableMean
if nargin < 3 || isempty(variableMean)
    variableMean = variableMean_def;
end

%check variableStd
if nargin < 4 || isempty(variableStd)
    variableStd = variableStd_def;
end

%check numModeMinMax
if nargin < 5 || isempty(numModeMinMax)
    numModeMinMax = numModeMinMax_def;
end

%check plotResults
if nargin < 6 || isempty(plotResults)
    plotResults = plotResults_def;
end

%check logData
if nargin < 7 || isempty(logData)
    logData = logData_def;
end
if logData && (variableMean==1&&variableStd~=1 || variableStd==1&&variableMean~=1)
    error('--analyzeIntensityModesMLMD: For log-normal fit,  mean and std must be either both variable or both constrained.')
end

%check input modes
if nargin < 8 || isempty(modeParamIn)
    modeParamIn = modeParamIn_def;
else
    numModeMinMax = size(modeParamIn,1)*[1 1];
end

%check what to analyze - amplitudes or intensities
if nargin < 9 || isempty(ampOrInt)
    ampOrInt = ampOrInt_def;
end

%check ratioTol
if nargin < 14 || isempty(ratioTol)
    ratioTol = ratioTol_def;
end

%% Analysis

%determine if input is a MovieList or MovieData object
if isa(MLMD,'MovieList') %if it is a movieDist
    
    listFlag = 1;
    
    %rename to ML
    ML = MLMD;
    clear MLMD
    
    %get number of movies
    numMovies = length(ML.movieDataFile_);
    
else %if it is a movieData
    
    listFlag = 0;
    
    %rename to MD
    MD = MLMD;
    clear MLMD
    
    %assign number of movies to 1
    numMovies = 1;
    
end

%go over all movies and run intensity modal analysis
for iM = 1 : numMovies
    
    %get movieData of current movie
    if listFlag == 1
        MD = MovieData.load(ML.movieDataFile_{iM});
    end
    
    %add intensity modal analysis process if never run before
    iProc = MD.getProcessIndex('IntModalAnalysisProcess',1,0);
    if isempty(iProc)
        iProc=numel(MD.processes_)+1;
        MD.addProcess(IntModalAnalysisProcess(MD));
    end
    %     if ~isempty(iProc)
    %         MD.deleteProcess(iProc);
    %         MD.addProcess(IntModalAnalysisProcess(MD));
    %     end
    
    %define function parameters
    
    %general
    funParams = MD.getProcess(iProc).funParams_;
    funParams.ChannelIndex = 1;
    
    %function-specific
    funParams.startFrame = 1;
    funParams.endFrame = [];
    funParams.alpha = alpha;
    funParams.variableMean = variableMean;
    funParams.variableStd = variableStd;
    funParams.numModeMinMax = numModeMinMax;
    funParams.plotResults = plotResults;
    funParams.logData = logData;
    funParams.modeParamIn = modeParamIn;
    funParams.ampOrInt = ampOrInt;
    funParams.ratioTol = ratioTol;
    
    %general again
    parseProcessParams(MD.getProcess(iProc),funParams);
    
    %% Run analysis processes
    cellfun(@run,MD.processes_(iProc));
    
end

