function correctMovieDisplacementField(movieData,varargin)
% correctMovieDisplacementField calculate the displacement field
%
% correctMovieDisplacementField 
%
% SYNOPSIS correctMovieDisplacementField(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Sep 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(DisplacementFieldCorrectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
displFieldCorrProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(displFieldCorrProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldCorrProc.getName());
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check displacement field process
iDisplFieldCalcProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,1);     
if isempty(iDisplFieldCalcProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])   
end

displFieldCalcProc=movieData.processes_{iDisplFieldCalcProc};
if ~displFieldCalcProc.checkChannelOutput
    error(['The channel must have a displacement field ! ' ...
        'Please calculate displacement field to all needed channels before '...
        'running force field calculation!'])
end

inFilePaths{1} = displFieldCalcProc.outFilePaths_{1};
displFieldCorrProc.setInFilePaths(inFilePaths);

% Set up the output directories
outputFile = {[p.OutputDirectory filesep 'displField.mat']};
mkClrDir(p.OutputDirectory);
displFieldCorrProc.setOutFilePaths(outputFile);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting correctiong displacement field...')
% Anonymous functions for reading input/output
displField=displFieldCalcProc.loadChannelOutput;


disp('Detecting and filtering vector field outliers...')
logMsg = 'Please wait, detecting and filtering vector field outliers';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

% %Parse input, store in parameter structure
% pd = parseProcessParams(displFieldCalcProc,paramsIn);

% Perform vector field outlier detection
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
for j= 1:nFrames
    % Outlier detection
    dispMat = [displField(j).pos displField(j).vec];
    % Take out duplicate points (Sangyoon)
    [dispMat,~,~] = unique(dispMat,'rows'); %dispMat2 = dispMat(idata,:),dispMat = dispMat2(iudata,:)
    displField(j).pos=dispMat(:,1:2);
    displField(j).vec=dispMat(:,3:4);

    if ~isempty(p.outlierThreshold)
        outlierIndex = detectVectorFieldOutliers(dispMat,p.outlierThreshold,1);
        %displField(j).pos(outlierIndex,:)=[];
        %displField(j).vec(outlierIndex,:)=[];
        dispMat(outlierIndex,3:4)=NaN;
        
        % I deleted this part for later gap-closing
        % Filter out NaN from the initial data (but keep the index for the
        % outliers)
%         ind= ~isnan(dispMat(:,3));
%         dispMat=dispMat(ind,:);
 
        displField(j).pos=dispMat(:,1:2);
        displField(j).vec=dispMat(:,3:4);

        % I deleted this part because artificially interpolated vector can
        % cause more error or false force. - Sangyoon June 2013
%         % Filling all NaNs with interpolated displacement vectors -
%         % We also calculate the interpolated displacements with a bigger correlation length.
%         % They are considered smoothed displacements at the data points. Sangyoon
%         dispMat = [dispMat(:,2:-1:1) dispMat(:,2:-1:1)+dispMat(:,4:-1:3)];
%         intDisp = vectorFieldSparseInterp(dispMat,...
%             displField(j).pos(:,2:-1:1),...
%             pd.minCorLength,pd.minCorLength,[],true);
%         displField(j).vec = intDisp(:,4:-1:3) - intDisp(:,2:-1:1);
    end
    
    % Update the waitbar
    if mod(j,5)==1 && feature('ShowFigureWindows')
        tj=toc;
        waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
    end
end

% Here, if nFrame>1, we do inter- and extrapolation of displacement vectors
% to prevent sudden, wrong force field change.
if nFrames>1
    disp('Performing displacement vector gap closing ...')
    % going through each point, see if there is NaN at each displacment
    % history and fill the gap
    logMsg = 'Performing displacement vector gap closing ...';

    nPoints = length(displField(1).pos(:,1));
    for k=1:nPoints
        % build each disp vector history
        curVecX = arrayfun(@(x) x.vec(k,1),displField);
        curVecY = arrayfun(@(x) x.vec(k,2),displField);
        if any(isnan(curVecX))
            t = 1:length(curVecX);
            t_nn = t(~isnan(curVecX));
            curVecX2 = interp1(t_nn,curVecX(~isnan(curVecX)),t,'linear');
            curVecY2 = interp1(t_nn,curVecY(~isnan(curVecX)),t,'linear');
            for ii=find(isnan(curVecX))
                displField(ii).vec(k,:) = [curVecX2(ii) curVecY2(ii)];
            end
        else
            continue
        end
        if mod(k,5)==1 && feature('ShowFigureWindows')
            tj=toc;
            waitbar(k/nPoints,wtBar,sprintf([logMsg timeMsg(tj*(nPoints-k)/k)]));
        end
    end
end

% Find rotational registration
if p.doRotReg, displField=perfRotReg(displField); end %#ok<NASGU>

save([p.OutputDirectory filesep 'displField.mat'],'displField');

% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished correcting displacement field!')