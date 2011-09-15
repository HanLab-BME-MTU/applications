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

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows'),
    wtBar = waitbar(0,'Initializing...','Name',forceFieldProc.getName());
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check displacement field process
iDisplFieldProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,1);     
if isempty(iDisplFieldProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])   
end

displFieldProc=movieData.processes_{iDisplFieldProc};
if ~displFieldProc.checkChannelOutput
    error(['The channel must have a displacement field ! ' ...
        'Please calculate displacement field to all needed channels before '...
        'running force field calculation!'])
end

inFilePaths{1} = displFieldProc.outFilePaths_{1};
forceFieldProc.setInFilePaths(inFilePaths);

% Set up the output directories

outputFile = [p.OutputDirectory filesep 'forceField.mat'];
mkClrDir(p.OutputDirectory);
forceFieldProc.setOutFilePaths(outputFile);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting calculating force  field...')
displField=displFieldProc.loadChannelOutput;

% Prepare displ for BEM
if strcmpi(p.method,'fastBEM')
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
[reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField);

forceField(nFrames)=struct('pos','','vec','','par','');

disp('Calculating force field...')
logMsg = 'Please wait, calculating force field';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;


for i=1:nFrames
    [grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
    
    if strcmpi(p.method,'FastBEM')
        % If grid_mat=[], then an optimal hexagonal force mesh is created
        % given the bead locations defined in displField:

        if i==1 || displField(i).par.prep4fastBEM==0;
            [pos_f, force, forceMesh, M, pos_u, u, sol_coef, sol_mats]=...
                reg_FastBEM_TFM(grid_mat, displField, i, ...
                p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM);
            display('The total time for calculating the FastBEM solution: ')
        elseif i>1 && displField(i).par.prep4fastBEM==1
            % since the displ field has been prepared such
            % that the measurements in different frames are ordered in the
            % same way, we don't need the position information any
            % more. The displ. measurements are enough.
            display('5.) Re-evaluate the solution:... ')
            % pull the new u-vector:
            u=vertcat(displField(i).vec(:,1),displField(i).vec(:,2));
            % recalculate the solution for the new displacement vec:
            [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,u,forceMesh,p.regParam,[],[]);
            display(['Done: solution for frame: ',num2str(i)]);
        end
        
        % The following values should/could be stored for the BEM-method.
        % In most cases, except the sol_coef this has to be stored only
        % once for all frames!
        if saveAllBEMpar==1 && ~savedSolMatsOnce
            forceField(i).par.forceMesh     = forceMesh;
            forceField(i).par.sol_coef      = sol_coef;
            forceField(i).par.M             = M; % This should not be saved every time! Although necessary to calculate the L-curve!
            forceField(i).par.sol_mats      = sol_mats; % This should not be saved every time! Although necessary to calculate the L-curve!
            forceField(i).par.pos           = pos_u;
            forceField(i).par.u             = u;  
            forceField(i).par.meshPtsFwdSol = meshPtsFwdSol;
            
            % don't save it again for the next frames. sol_coef can be
            % easily obtained from calcSolFromSolMatsFastBEM!
            savedSolMatsOnce = 1;
        elseif saveAllBEMpar==1 && displField(i).par.prep4fastBEM==0
            display('!!! The solution matrices cannot be saved for all frames (~1GB per frame)!!!')
            display('!!!        Choose manually the ones you are intereseted in               !!!')
        end
    else
        [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, yModu_Pa, pRatio, pixSize_mu, gridSpacing, i_max, j_max, regParam);
    end   
    
    
    % Fill in the values to be stored:
    forceField(i).pos=pos_f;
    forceField(i).vec=force;
    clear grid_mat;
    clear iu;
    clear iu_mat;
    
    % Update the waitbar
    if mod(i,5)==1 && feature('ShowFigureWindows')
        ti=toc;
        waitbar(i/nFrames,wtBar,sprintf([logMsg timeMsg(ti*nFrames/i-ti)]));
    end
end

save(outputFile,'forceField');

% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished calculating force field!')
