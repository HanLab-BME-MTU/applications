function [fitResults,errFlag] = armaxFitKalman(trajOut,trajIn,modelParamOrOrder,...
    minOpt)
%ARMAXFITKALMAN fits a number of ARMAX models to a time series which could depend on another time series.
%
%SYNOPSIS [fitResults,errFlag] = armaxFitKalman(trajOut,trajIn,modelParamOrOrder,...
%    minOpt)
%
%INPUT
%   Mandatory
%       trajOut          : Observations of output time series to be fitted. An 
%                          array of structures with fields:
%               .observations: 2D array of measurements and their uncertainties.
%                              Missing points should be indicated with NaN.
%               .weight      : The weight with which each movie belongs to the
%                              group. Optional. If field doesn't exist,
%                              all weights will be taken as 1.
%                          If there is only one series, it can also be
%                          input directly as a 2D array. In this case, its
%                          weight is 1.
%   Optional
%       trajIn           : Observations of input time series to be fitted. Either an 
%                          array of structures trajIn(1:nTraj).observations, or a
%                          2D array representing one single trajectory. 
%               .observations: 2D array of measurements and their uncertainties.
%                          Must not have any missing points.
%                          Enter as [] if there is no input series.
%                          Default: [].
%       modelParamOrOrder: Either
%                          (1)3x2 matrix indicating range of orders of models:
%                             [lowest AR order   highest AR order
%                              lowest MA order   highest MA order
%                              lowest X  order   highest X  order]
%                             Use 0, 0 and -1 for arOrder, maOrder and xOrder,
%                             respectively, to indicate the lack of any depedence.
%                          Or
%                          (2)1, 2 or 3D structure array of models to be tested, 
%                             where each entry consists of 3 elements (row vectors):
%               .arParamP0   : Initial guess of parameters related to 
%                              partial AR coefficients.
%               .maParamP0   : Initial guess of parameters related to 
%                              partial MA coefficients.
%               .xParam0     : Initial guess of X coefficients.
%                          Default: order range with input [0 3; 0 3; -1 3],
%                          order range without input [0 3; 0 3; -1 -1], with
%                          initial guesses of parameters picked randomly.
%       minOpt           : Minimization option: 
%                          -'ml' for Matlab local minimizer "fmincon";
%                          -'tl' for Tomlab local minimizer "ucSolve";
%                          -'tg' for Tomlab global minimizer "glbFast"' followed
%                           by Tomlab local minimizer "ucSolve"; -- DON'T USE
%                          -'nag' for NAG's local minimizer E04JAF. -- DON'T USE
%                          Default: 'tl'.
%
%OUTPUT fitResults: Structure array containing fitting results. 
%                   Will have the dimensions indicated in modelParamOrOrder. 
%                   Contains the fields:
%           .arParamP0 : Initial guess of parameters related to partial AR coefficients.
%           .maParamP0 : Initial guess of parameters related to partial MA coefficients.
%           .xParam0   : Initial guess of X coefficients.
%           .arParamK  : Estimated AR coefficients (1st row) and parameters 
%                        related to partial AR coefficients (2nd row)
%                        using likelihood maximization.
%           .maParamK  : Estimated MA coefficients (1st row) and parameters 
%                        related to partial MA coefficients (2nd orw) 
%                        using likelihood maximization.
%           .xParamK   : Estimated X coefficients using likelihood maximization.
%           .arParamL  : Estimated AR coefficients using least squares fitting.
%           .maParamL  : Estimated MA coefficients using least squares fitting.
%           .xParamL   : Estimated X coefficients using least squares fitting.
%           .varCovMatL: Variance-covariance matrix of ARMAX coefficients, 
%                        estimated via least squares fitting.
%           .varCovMatF: Variance-covariance matrix of ARMAX coefficients,
%                        estimated via the Fisher Information Matrix.
%           .wnVariance: Estimated variance of white noise in process.
%           .wnVector  : Structure array containing the field:
%               .observations: Estimated white noise series in corresponding 
%                              trajectory.
%           .selectCrit: Model selection criterion. Contains 3 fields:
%               .aic      : Akaike's Information Criterion.
%               .aicc     : Bias-Corrected Akaike's Information Criterion.
%               .bic      : Bayesian Information Criterion.
%           .pVCompKL  : P-value for comparing maximum likelihood estimate
%                        with least squares estimate.
%           .pVPort    : P-value from portmanteau test of the uncorrelatedness
%                        of the white noise series.
%           .numObserve: Number of observations used in estimation.
%           .success   : 1 if fitting executed normally, 0 otherwise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARKS Data is fitted using armaxCoefKalman beginning with the different
%        input models.
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--armaxFitKalman: You should AT LEAST supply the output series!');
    errFlag = 1;
    fitResults = [];
    return
end

%check "trajOut" and turn into struct if necessary
weightZero = false; % remembers whether any output-trajectory was removed
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    numTraj = 1; %number of trajectories supplied
    trajOut.weight = 1; %weight indicating association with group
    clear tmp
else
    if ~isfield(trajOut,'observations')
        disp('--armaxFitKalman: Please input trajOut in fields ''observations''!')
        errFlag = 1;
    end
    numTraj = length(trajOut); %number of trajectories supplied
    %assign default weights if not supplied
    if ~isfield(trajOut,'weight')
        [trajOut.weight] = deal(1);
    else
        % loop and check whether there are empty weights
        for i=1:numTraj
            if isempty(trajOut(i).weight)
               trajOut(i).weight = 1;
            end
        end
        
        % check for weight==0 and remove the data
        weights = cat(1,trajOut.weight);
        noWeightIdx = weights == 0;
        if any(noWeightIdx)
            weightZero = true;
            trajOut(noWeightIdx) = [];
            % make sure that there are trajectories left
            if isempty(trajOut)
                disp('--armaxFitKalman: All trajectories have weight 0!')
                errFlag = 1;
            end
        end
        
    end
end

%check "trajIn" and turn into struct if necessary
if nargin < 2 || isempty(trajIn)
    trajIn = [];
else
    if ~isstruct(trajIn)
        tmp = trajIn;
        clear trajIn
        trajIn.observations = tmp;
        clear tmp
    elseif ~isfield(trajIn,'observations')
        disp('--armaxFitKalman: Please input trajIn in fields ''observations''!')
        errFlag = 1;
    end
    
    % if there were output trajectories with zero weight: remove
    if weightZero
        trajIn(noWeightIdx) = [];
    end
end

%assign default values of optional variables
orderValArma_def = 0:3; %AR and MA order values
if isempty(trajIn)
    orderValX_def = -1; %X order values
else
    orderValX_def = -1:3; %X order values
end
repeat_def = 10; %times to repeat local minimization if initial guess not supplied
minOpt_def = 'tl'; %minimization option

%check models to be tested
if nargin < 3 || isempty(modelParamOrOrder) %if no models to fit were input

    %indicate that initial guess was not supplied
    suppliedIG = 0;

    %randomly generate initial ARMAX parameter guesses
    for l=1:repeat_def
        for k=length(orderValX_def):-1:1
            for j=length(orderValArma_def):-1:1
                for i=length(orderValArma_def):-1:1
                    modelParam(i,j,k,l).arParamP0 = 4*rand(1,orderValArma_def(i)) - 2;
                    modelParam(i,j,k,l).maParamP0 = 4*rand(1,orderValArma_def(j)) - 2;
                    modelParam(i,j,k,l).xParam0 = 4*rand(1,orderValX_def(k)+1) - 2;
                end
            end
        end
    end
        
else %if models or model orders were supplied

    if isstruct(modelParamOrOrder) %if initial guesses were input

        %indicate that initial guess was supplied
        suppliedIG = 1;
        
        %assign model parameters
        modelParam = modelParamOrOrder;
        
        %place [] in fields not provided
        if ~isfield(modelParam,'arParamP0')
            modelParam(1).arParamP0 = [];
        end
        if ~isfield(modelParam,'maParamP0')
            modelParam(1).maParamP0 = [];
        end
        if ~isfield(modelParam,'xParam0')
            modelParam(1).xParam0 = [];
        end

    else %if range of AR and MA orders was input

        %indicate that initial guess was not supplied
        suppliedIG = 0;
        
        %check size of matrix
        [nRow,nCol] = size(modelParamOrOrder);
        if nRow ~= 3 || nCol ~= 2
            disp('--armaxFitKalman: Matrix indicating range of AR, MA and X orders should be a 3x2 matrix!');
            errFlag = 1;
        else
            
            %generate random initial parameter guesses for the order range specified
            for k=modelParamOrOrder(3,1):modelParamOrOrder(3,2)
                k1 = k - modelParamOrOrder(3,1) + 1;
                for j=modelParamOrOrder(2,1):modelParamOrOrder(2,2)
                    j1 = j - modelParamOrOrder(2,1) + 1;
                    for i=modelParamOrOrder(1,1):modelParamOrOrder(1,2)
                        i1 = i - modelParamOrOrder(1,1) + 1;

                        for l=1:repeat_def
                            modelParam(i1,j1,k1,l).arParamP0 = 4*rand(1,i) - 2;
                            modelParam(i1,j1,k1,l).maParamP0 = 4*rand(1,j) - 2;
                            modelParam(i1,j1,k1,l).xParam0 = 4*rand(1,k+1) - 2;
                        end

                    end
                end
            end
        end %(if nRow ~= 3 || nCol ~= 2 ... else ...)

    end %(if isstruct(modelParamOrOrder) ... else ...)

end %(if nargin < 3 || isempty(modelParamOrOrder) ... else ...)

%check minimization option
if nargin < 4 || isempty(minOpt) %if minimization option was not input
    minOpt = minOpt_def; %assign default
else %if minimization option was input, check its value
%     if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl') ...
%             && ~strcmp(minOpt,'tg') && ~strcmp(minOpt,'nag'))
%         disp('--armaxFitKalman: "minOpt" should be either "ml", "tl", "tg" or "nag"!');
%         errFlag = 1;
%     end
    if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl'))
        disp('--armaxFitKalman: "minOpt" should be either "ml" or "tl"!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--armaxFitKalman: Please fix input data!');
    fitResults = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% turn off warning:divide by zero
warningState = warning;
warning('off','MATLAB:divideByZero');

%initialize progress text
progressText(0,'fitting ARMA models');

%get total number of runs
numRunsTot = length(modelParam(:));

%initialize variable counting number of runs
iRun = 0;

%go over all suggested models
for k=size(modelParam,3):-1:1
    for j=size(modelParam,2):-1:1
        for i=size(modelParam,1):-1:1

            if strcmp(minOpt,'tg') || suppliedIG %if global minimization or if initial guess was supplied

                %assign model data
                arParamP0 = modelParam(i,j,k,1).arParamP0;
                maParamP0 = modelParam(i,j,k,1).maParamP0;
                xParam0   = modelParam(i,j,k,1).xParam0;

                %estimate ARMA coeffients and white noise variance
                [arParamK,maParamK,xParamK,arParamL,maParamL,xParamL,...
                    varCovMatL,varCovMatF,wnVariance,wnVector,selectCrit,...
                    pVCompKL,pVPort,errFlag] = armaxCoefKalman(trajOut,...
                    trajIn,arParamP0,maParamP0,xParam0,[],minOpt);

                %display progress
                iRun = iRun + 1;
                progressText(iRun/numRunsTot,'fitting ARMA models');

            else %if local minimization and initial guess was not supplied

                %anitialize output
                neg2LL = 10^20;
                arParamK   = [];
                maParamK   = [];
                xParamK    = [];
                arParamL   = [];
                maParamL   = [];
                xParamL    = [];
                varCovMatL = [];
                varCovMatF = [];
                wnVariance = [];
                wnVector   = [];
                selectCrit = [];
                pVCompKL   = [];
                pVPort     = [];
                errFlag    = [];

                for l=1:repeat_def

                    %assign model data
                    arParamP0 = modelParam(i,j,k,l).arParamP0;
                    maParamP0 = modelParam(i,j,k,l).maParamP0;
                    xParam0   = modelParam(i,j,k,l).xParam0;

                    %estimate ARMA coeffients and white noise variance
                    [arParamK1,maParamK1,xParamK1,arParamL1,maParamL1,...
                        xParamL1,varCovMatL1,varCovMatF1,wnVariance1,...
                        wnVector1,selectCrit1,pVCompKL1,pVPort1,errFlag1]...
                        = armaxCoefKalman(trajOut,trajIn,arParamP0,...
                        maParamP0,xParam0,[],minOpt);

                    if ~isempty(selectCrit1)

                        %get -2ln(likelihood) of estimated model
                        neg2LL1 = selectCrit1.aic - ...
                            2*(length([arParamP0 maParamP0 xParam0])+1);

                        %choose this as the best model if its
                        %-2ln(likelihood) is the smallest so far
                        if neg2LL1 < neg2LL
                            neg2LL     = neg2LL1;
                            arParamK   = arParamK1;
                            maParamK   = maParamK1;
                            xParamK    = xParamK1;
                            arParamL   = arParamL1;
                            maParamL   = maParamL1;
                            xParamL    = xParamL1;
                            varCovMatL = varCovMatL1;
                            varCovMatF = varCovMatF1;
                            wnVariance = wnVariance1;
                            wnVector   = wnVector1;
                            selectCrit = selectCrit1;
                            pVCompKL   = pVCompKL1;
                            pVPort     = pVPort1;
                            errFlag    = errFlag1;
                        end

                    end

                    %display progress
                    iRun = iRun + 1;
                    progressText(iRun/numRunsTot,'fitting ARMA models');
                    
                end %(for l=1:repeat_def)

            end %(if strcmp(minOpt,'tg') || suppliedIG ... else ...)

            %get number of observations used in parameter estimation
            tmp = vertcat(wnVector.observations);
            numObserve = length(find(~isnan(tmp(:,1))));

            %write output as fields in fitResults
            fitResults(i,j,k) = struct(...
                'arParamP0',arParamP0,...
                'maParamP0',maParamP0,...
                'xParam0',xParam0,...
                'arParamK',arParamK,...
                'maParamK',maParamK,...
                'xParamK',xParamK,...
                'arParamL',arParamL,...
                'maParamL',maParamL,...
                'xParamL',xParamL,...
                'varCovMatL',varCovMatL,...
                'varCovMatF',varCovMatF,...
                'wnVariance',wnVariance,...
                'wnVector',wnVector,...
                'selectCrit',selectCrit,...
                'pVCompKL',pVCompKL,...
                'pVPort',pVPort,...
                'numObserve',numObserve,...
                'success',1-errFlag);

        end %(for k=1:size(modelParam,3))
    end %(for j=1:size(modelParam,2))
end %(for i=1:size(modelParam,1))


% put warnings back to what they were before
warning(warningState);


%%%%% ~~ the end ~~ %%%%%
