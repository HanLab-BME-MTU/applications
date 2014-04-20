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
p.saveBEMparams = false;
p.useLcurve = true;
p.lastToFirst = false;
%% --------------- Initialization ---------------%%
% if feature('ShowFigureWindows'),
%     wtBar = waitbar(0,'Initializing...','Name',forceFieldProc.getName());
%     wtBarArgs={'wtBar',wtBar};
% else
    wtBar=-1;
    wtBarArgs={};
% end

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
outputFile{2,1} = [p.OutputDirectory filesep 'BEMParams.mat'];
% if  ~strcmp(p.solMethodBEM,'QR')
if p.useLcurve
    outputFile{3,1} = [p.OutputDirectory filesep 'Lcurve.fig'];
    outputFile{4,1} = [p.OutputDirectory filesep 'LcurveData.mat'];
end
mkClrDir(p.OutputDirectory);
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
end
   
% Prepare displacement field for BEM
if strcmpi(p.method,'fastBEM')
    displField(end).par=0; % for compatibility with Achim parameter saving
    displField=prepDisplForBEM(displField,'linear');
end


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
    [reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField,1.5,0); %denser force mesh for ROI
end

forceField(nFrames)=struct('pos','','vec','','par','');

disp('Calculating force field...')
logMsg = 'Please wait, calculating force field';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

if p.lastToFirst 
    frameSequence = nFrames:-1:1;
else
    frameSequence = 1:nFrames;
end
jj = 0;
for i=frameSequence
    jj=jj+1;
    [grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
    
%     if min(min(maskArray(:,:,1))) == 0 %when ROI was used, make a denser displacement field by interpolation
%         % vectorized position
%         posx=reshape(grid_mat(:,:,1),[],1);
%         posy=reshape(grid_mat(:,:,2),[],1);
% 
%         dispMat = [displField(i).pos(:,2:-1:1) displField(i).pos(:,2:-1:1)+displField(i).vec(:,2:-1:1)];
%         intDisp = vectorFieldSparseInterp(dispMat,[posy posx],20,20,[]);
%         displField(i).vec = intDisp(:,4:-1:3) - intDisp(:,2:-1:1);
%         displField(i).pos = [posx posy];
%         
% %         vecx=reshape(iu_mat(:,:,1),[],1);
% %         vecy=reshape(iu_mat(:,:,2),[],1);
% %         displField(i).vec =  [vecx vecy];
%     end

%     posx=reshape(grid_mat(:,:,1),[],1);
%     posy=reshape(grid_mat(:,:,2),[],1);
%     displField(i).pos = [posx posy];
%     vecx=reshape(iu_mat(:,:,1),[],1);
%     vecy=reshape(iu_mat(:,:,2),[],1);
%     displField(i).vec =  [vecx vecy];

    if strcmpi(p.method,'FastBEM')
        % If grid_mat=[], then an optimal hexagonal force mesh is created
        % given the bead locations defined in displField:
        if p.usePaxImg && length(movieData.channels_)>1
            paxImage=movieData.channels_(2).loadImage(i);
            [pos_f, force, forceMesh, M, pos_u, u, sol_coef, sol_mats]=...
                reg_FastBEM_TFM(grid_mat, displField, i, ...
                p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                'thickness',p.thickness/movieData.pixelSize_,'paxImg',paxImage,'pixelSize',movieData.pixelSize_);
            
            outputFile{3+i,1} = [p.OutputDirectory filesep 'BEMParams ' num2str(i) ' frame.mat'];

            disp(['saving forward map and custom force mesh at ' outputFile{3+i,1} '...'])
            save(outputFile{3+i,1},'forceMesh','M','-v7.3');
            display(['Done: solution for frame: ',num2str(i)]);
            
        else
            if jj==1
%                 if ishandle(wtBar)
%                     waitbar(0,wtBar,sprintf([logMsg ' for first frame']));
%                 end
                if p.useLcurve
                    [pos_f, force, forceMesh, M, pos_u, u, sol_coef, sol_mats]=...
                        reg_FastBEM_TFM(grid_mat, displField, i, ...
                        p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                        'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                        'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                        'useLcurve',p.useLcurve,'thickness',p.thickness/movieData.pixelSize_,...
                        'LcurveDataPath',outputFile{4,1},'LcurveFigPath',outputFile{3,1});
                    p.regParam = sol_mats.L;
                else
                    [pos_f, force, forceMesh, M, pos_u, u, sol_coef, sol_mats]=...
                        reg_FastBEM_TFM(grid_mat, displField, i, ...
                        p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                        'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                        'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                        'useLcurve',p.useLcurve,'thickness',p.thickness/movieData.pixelSize_);
                end
                %             display('The total time for calculating the FastBEM solution: ')
                
                % The following values should/could be stored for the BEM-method.
                % In most cases, except the sol_coef this has to be stored only
                % once for all frames!
                if p.saveBEMparams
                    disp(['saving forward map and force mesh at ' outputFile{2} '...'])
                    save(outputFile{2},'forceMesh','M','sol_mats','pos_u','u','-v7.3');
                end
                
                %             % Calculate L-curve
                %             if ~strcmp(p.solMethodBEM,'QR')
                %                 hLcurve=plotLcurve(M,sol_mats,u,forceMesh,p.LcurveFactor,...
                %                     'wtBar',wtBar);
                %                 saveas(hLcurve,outputFile{3});
                %                 close(hLcurve)
                %             end
            else
                % since the displ field has been prepared such
                % that the measurements in different frames are ordered in the
                % same way, we don't need the position information any
                % more. The displ. measurements are enough.
                display('5.) Re-evaluate the solution:... ')
                %             % pull the new u-vector:
                %             u=vertcat(displField(i).vec(:,1),displField(i).vec(:,2));
                % recalculate the solution for the new displacement vec:
                if p.usePaxImg && length(movieData.channels_)>1
                    paxImage=movieData.channels_(2).loadImage(i);
                    [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                        displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,p.regParam,[],[], paxImage);
                else
                    [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                        displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,p.regParam,[],[]);
                end
                display(['Done: solution for frame: ',num2str(i)]);
            end
        end
    else
        [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, p.YoungModulus,...
            p.PoissonRatio, movieData.pixelSize_/1000, gridSpacing, i_max, j_max, p.regParam);
    end   
    
    
    % Fill in the values to be stored:
    forceField(i).pos=pos_f;
    forceField(i).vec=force;
    clear grid_mat;
    clear iu;
    clear iu_mat;
    
    % Update the waitbar
%     if mod(i,5)==1 && ishandle(wtBar)
%         ti=toc;
%         waitbar(i/nFrames,wtBar,sprintf([logMsg timeMsg(ti*nFrames/i-ti)]));
%     end
end

save(outputFile{1},'forceField');

% Close waitbar
% if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished calculating force field!')
