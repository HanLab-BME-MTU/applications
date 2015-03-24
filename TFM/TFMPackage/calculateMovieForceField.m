function calculateMovieForceField(movieData,varargin)
% calculateMovieForceField calculate the displacement field
%
% calculateMovieForceField 
%
% SYNOPSIS calculateMovieForceField(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Sep 2011
% Last updated by Sangyoon Han, Oct 2014

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('ForceFieldCalculationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ForceFieldCalculationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
forceFieldProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(forceFieldProc,paramsIn);
p.usePaxImg = false;
p.saveBEMparams = true;
p.lastToFirst = false;
p.LcurveFactor = 10;
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows'),
    wtBar = waitbar(0,'Initializing...','Name',forceFieldProc.getName());
    wtBarArgs={'wtBar',wtBar};
else
    wtBar=-1;
    wtBarArgs={};
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check optional process Displacement field correction first
iDisplFieldProc =movieData.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);     
if isempty(iDisplFieldProc)
    iDisplFieldProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
end

if isempty(iDisplFieldProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])
end

displFieldProc=movieData.processes_{iDisplFieldProc};

if ~displFieldProc.checkChannelOutput()
    error(['Missing displacement field ! Please apply displacement field '...
        'calculation/correction  before running force field calculation!'])
end

% define resolution depending on the grid information in displacementField
% step
iDisplFieldCalProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
displFieldCalProc=movieData.processes_{iDisplFieldCalProc};
pDisp = parseProcessParams(displFieldCalProc);
try
    pDisp.useGrid;
catch
    pDisp.useGrid = false;
end
if pDisp.useGrid
    p.highRes = false;
else
    p.highRes = true;
end
% Set up the input file
inFilePaths{1} = displFieldProc.outFilePaths_{1};
forceFieldProc.setInFilePaths(inFilePaths);

% Set up the output file
outputFile{1,1} = [p.OutputDirectory filesep 'forceField.mat'];
outputFile{2,1} = [p.OutputDirectory filesep 'tractionMaps.mat'];
outputFile{3,1} = [p.OutputDirectory filesep 'BEMParams.mat'];
% if  ~strcmp(p.solMethodBEM,'QR')
if p.useLcurve
    outputFile{4,1} = [p.OutputDirectory filesep 'Lcurve.fig'];
    outputFile{5,1} = [p.OutputDirectory filesep 'LcurveData.mat'];
end

% Add a recovery mechanism if process has been stopped in the middle of the
% computation to re-use previous results
firstFrame =1; % Set the strating fram eto 1 by default
if exist(outputFile{1},'file')
    % Check analyzed frames
    s=load(outputFile{1},'forceField');
    frameForceField=~arrayfun(@(x)isempty(x.pos),s.forceField);
    
    if ~all(frameForceField) && ~all(~frameForceField)
        % Look at the first non-analyzed frame
        firstFrame = find(~frameForceField,1);
        % Ask the user if display mode is active
        if ishandle(wtBar),
            recoverRun = questdlg(...
                ['A force field output has been dectected with ' ...
                num2str(firstFrame-1) ' analyzed frames. Do you' ...
                ' want to use these results and continue the analysis'],...
                'Recover previous run','Yes','No','Yes');
            if ~strcmpi(recoverRun,'Yes'), firstFrame=1; end
        end
    end
end
% asking if you want to reuse the fwdMap again (you have to make sure that
% you are solving the same problem with only different reg. param.) -SH
reuseFwdMap = 'No';
if strcmpi(p.method,'FastBEM') && exist(outputFile{3,1},'file')
    reuseFwdMap = questdlg(...
        ['BEM parameters were dectected. Do you' ...
        ' want to use these parameter and overwrite the results?'],...
        'Reuse Fwdmap','Yes','No','No');
end

% Backup the original vectors to backup folder
if firstFrame==1 && strcmpi(reuseFwdMap,'No')
    display('Backing up the original data')
    backupFolder = [p.OutputDirectory ' Backup']; % name]);
    if exist(p.OutputDirectory,'dir')
        ii = 1;
        while exist(backupFolder,'dir')
            backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
            ii=ii+1;
        end
        mkdir(backupFolder);
        copyfile(p.OutputDirectory, backupFolder,'f')
    end
    forceField(nFrames)=struct('pos','','vec','','par','');
    mkClrDir(p.OutputDirectory);
    M = [];
elseif strcmpi(reuseFwdMap,'Yes')
    fwdMapFile = load(outputFile{3,1},'M');
    M = fwdMapFile.M;
else
    % Load old force field structure 
    forceField=s.forceField;
    M = [];
end

forceFieldProc.setOutFilePaths(outputFile);

%% --------------- Force field calculation ---------------%%% 

disp('Starting calculating force  field...')
maskArray = movieData.getROIMask;
if min(min(maskArray(:,:,1))) == 0
    iStep2Proc = movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);
    step2Proc = movieData.processes_{iStep2Proc};
    pDisp = parseProcessParams(step2Proc,paramsIn);

    % Use mask of first frame to filter displacementfield
    iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(iSDCProc)
        SDCProc=movieData.processes_{iSDCProc};
        if ~SDCProc.checkChannelOutput(pDisp.ChannelIndex)
            error(['The channel must have been corrected ! ' ...
                'Please apply stage drift correction to all needed channels before '...
                'running displacement field calclation tracking!'])
        end
        %Parse input, store in parameter structure
        refFrame = double(imread(SDCProc.outFilePaths_{2,pDisp.ChannelIndex}));
        firstMask = false(size(refFrame));
        tempMask = maskArray(:,:,1);
        firstMask(1:size(tempMask,1),1:size(tempMask,2)) = tempMask;
        displFieldOriginal=displFieldProc.loadChannelOutput;
        displField = filterDisplacementField(displFieldOriginal,firstMask);
    else
        firstMask = maskArray(:,:,1);
        displFieldOriginal=displFieldProc.loadChannelOutput;
        displField = filterDisplacementField(displFieldOriginal,firstMask);
    end        
else
     displField=displFieldProc.loadChannelOutput;
     firstMask = maskArray(:,:,1);
end
   
% % Prepare displacement field for BEM
% if strcmpi(p.method,'fastBEM')
%     displField(end).par=0; % for compatibility with Achim parameter saving
%     displField=prepDisplForBEM(displField,'linear');
% end
% I don't think this is necessary anymore. - SH

% For Benedikt's software to work, the displacement field has to be
% interpolated on a rectangular grid, with an even number of grid points
% along each edge. Furthermore, one has to make sure that the noisy data 
% has not to be extrapolated. This may happen along the edges. To prevent
% this, extract the corner of the displacement grid, calculate how often 
% (even number) the optimal gridspacing fits into each dimension, then 
% place the regular grid centered to the orignal bounds. Thereby make sure 
% that the edges have been eroded to a certain extend. This is performed by
% the following function.
if ~p.highRes
    [reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField,1,0);
else
    [reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField,1.0,0); %no dense mesh in any case. It causes aliasing issue!
end
tMap = cell(1,nFrames);

disp('Calculating force field...')
logMsg = 'Please wait, calculating force field';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

if p.lastToFirst 
    frameSequence = nFrames:-1:firstFrame;
else
    frameSequence = firstFrame:nFrames;
end
jj = 0;
if strcmpi(p.method,'FastBEM')
    % if FastBEM, we calculate forward map and mesh only in the first frame
    % and then use parfor for the rest of the frames to calculate forces -
    % SH
    i=frameSequence(1); % For the first frame
    [grid_mat,~, ~,~] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
    
    % If grid_mat=[], then an optimal hexagonal force mesh is created
    % given the bead locations defined in displField:
    if p.usePaxImg && length(movieData.channels_)>1
        for i=frameSequence
            paxImage=movieData.channels_(2).loadImage(i);
            [pos_f, force, forceMesh, M, ~, ~, ~, ~]=...
                reg_FastBEM_TFM(grid_mat, displField, i, ...
                p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                'thickness',p.thickness/movieData.pixelSize_,'paxImg',paxImage,'pixelSize',movieData.pixelSize_);

            outputFile{3+i,1} = [p.OutputDirectory filesep 'BEMParams ' num2str(i) ' frame.mat'];

            disp(['saving forward map and custom force mesh at ' outputFile{3+i,1} '...'])
            save(outputFile{3+i,1},'forceMesh','M','-v7.3');
            display(['Done: solution for frame: ',num2str(i)]);
            % Fill in the values to be stored:
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            save(outputFile{1},'forceField');
        end
    elseif i>1
        display('Loading BEM parameters... ')
        load(outputFile{3},'forceMesh','M','sol_mats');
        for i=frameSequence
            if p.usePaxImg && length(movieData.channels_)>1
                paxImage=movieData.channels_(2).loadImage(i);
                [pos_f,force,~]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,p.regParam,[],[], 'paxImg', paxImage, 'useLcurve', p.useLcurve);
            else
                [pos_f,force,~]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,sol_mats.L,[],[]);
            end
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            % Save each iteration (for recovery of unfinished processes)
            save(outputFile{1},'forceField');
            display(['Done: solution for frame: ',num2str(i)]);
            %Update the waitbar
            if mod(i,5)==1 && ishandle(wtBar)
                ti=toc;
                waitbar(i/nFrames,wtBar,sprintf([logMsg timeMsg(ti*nFrames/i-ti)]));
            end
        end
    else
                %         if ishandle(wtBar)
%             waitbar(0,wtBar,sprintf([logMsg ' for first frame']));
%         end
        if p.useLcurve
            [pos_f, force, forceMesh, M, pos_u, u, sol_coef,  sol_mats]=...
                reg_FastBEM_TFM(grid_mat, displField, i, ...
                p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                'useLcurve',p.useLcurve>0, 'LcurveFactor',p.LcurveFactor,'thickness',p.thickness/movieData.pixelSize_,...
                'LcurveDataPath',outputFile{5,1},'LcurveFigPath',outputFile{4,1},'fwdMap',M,...
                'lcornerOptimal',p.lcornerOptimal);
            params = parseProcessParams(forceFieldProc,paramsIn);
            params.regParam = sol_mats.L;
            p.regParam = sol_mats.L;
            forceFieldProc.setPara(params);
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            save(outputFile{1},'forceField');
        else
            [pos_f, force, forceMesh, M, pos_u, u, sol_coef,  sol_mats]=...
                reg_FastBEM_TFM(grid_mat, displField, i, ...
                p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                'useLcurve',p.useLcurve>0, 'thickness',p.thickness/movieData.pixelSize_,'fwdMap',M);
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            save(outputFile{1},'forceField');
        end
        %             display('The total time for calculating the FastBEM solution: ')

        % The following values should/could be stored for the BEM-method.
        % In most cases, except the sol_coef this has to be stored only
        % once for all frames!
        if p.saveBEMparams && strcmpi(reuseFwdMap,'No')
            disp(['saving forward map and force mesh at ' outputFile{3} '...'])
            save(outputFile{3},'forceMesh','M','sol_mats','pos_u','u','-v7.3');
        end
        for i=frameSequence(2:end)
            % since the displ field has been prepared such
            % that the measurements in different frames are ordered in the
            % same way, we don't need the position information any
            % more. The displ. measurements are enough.
            display('5.) Re-evaluate the solution:... ')
            if p.usePaxImg && length(movieData.channels_)>1
                paxImage=movieData.channels_(2).loadImage(i);
                [pos_f,force,~]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,p.regParam,[],[], 'paxImg', paxImage, 'useLcurve', p.useLcurve);
            else
                [pos_f,force,~]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,sol_mats.L,[],[]);
            end
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            % Save each iteration (for recovery of unfinished processes)
            save(outputFile{1},'forceField');
            display(['Done: solution for frame: ',num2str(i)]);
            %Update the waitbar
            if mod(i,5)==1 && ishandle(wtBar)
                ti=toc;
                waitbar(i/nFrames,wtBar,sprintf([logMsg timeMsg(ti*nFrames/i-ti)]));
            end
        end
    end
else % FTTC
    for i=frameSequence
        [grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
        if p.useLcurve
            [rho,eta,reg_corner,alphas] = calculateLcurveFTTC(grid_mat, iu_mat, p.YoungModulus,...
                p.PoissonRatio, gridSpacing, i_max, j_max, p.regParam,p.LcurveFactor);
            [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, p.YoungModulus,...
                p.PoissonRatio, movieData.pixelSize_/1000, gridSpacing, i_max, j_max, reg_corner);
            [reg_corner,ireg_corner,~,hLcurve]=regParamSelecetionLcurve(rho,eta,alphas,reg_corner,'manualSelection',true);
            save(outputFile{5,1},'rho','eta','reg_corner','ireg_corner');
            saveas(hLcurve,outputFile{4,1});
            close(hLcurve)
        else
            [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, p.YoungModulus,...
                p.PoissonRatio, movieData.pixelSize_/1000, gridSpacing, i_max, j_max, p.regParam);
        end
        forceField(i).pos=pos_f;
        forceField(i).vec=force;
    end
end
% For calculation of traction map
% The drift-corrected frames should have independent channel
% ->StageDriftCorrectionProcess
[tMapIn, tmax, tmin, cropInfo] = generateHeatmapShifted(forceField,displField,0);
% Insert traction map in forceField.pos 
disp('Generating traction maps ...')
for ii=frameSequence
    % starts with original size of beads
    cur_tMap = zeros(size(firstMask));
    cur_tMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{ii};
    tMap{ii} = cur_tMap;
end

% Fill in the values to be stored:
clear grid_mat;
clear iu;
clear iu_mat;
                
disp('Saving ...')
save(outputFile{1},'forceField');
save(outputFile{2},'tMap'); % need to be updated for faster loading. SH 20141106
forceFieldProc.setTractionMapLimits([tmin tmax])

% Close waitbar
% if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished calculating force field!')
end
