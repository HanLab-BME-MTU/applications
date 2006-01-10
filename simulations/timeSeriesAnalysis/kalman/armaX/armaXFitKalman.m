function [fitResults,errFlag] = armaXFitKalman(trajOut,trajIn,modelParamOrOrder,...
    minOpt)
%ARMAXFITKALMAN fits a number of ARMAX models to a time series which could depend on another time series.
%
%SYNOPSIS [fitResults,errFlag] = armaXFitKalman(trajOut,trajIn,modelParamOrOrder,...
%    minOpt)
%
%INPUT  trajOut          : Observations of output time series to be fitted. Either an 
%                          array of structures trajOut(1:nTraj).observations, or a
%                          2D array representing one single trajectory. 
%               .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%       trajIn           : Observations of input time series to be fitted. Either an 
%                          array of structures trajIn(1:nTraj).observations, or a
%                          2D array representing one single trajectory. 
%               .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%       modelParamOrOrder: Either
%                          (1)3x2 matrix indicating range of orders of models to fit:
%                             [lowest AR order   highest AR order
%                              lowest MA order   highest MA order
%                              lowest X  order   highest X  order]
%                          Or
%                          (2)1D or 2D structure array of models to be tested, 
%                             where each entry consists of 3 elements (row vectors):
%               .arParamP0   : Initial guess of parameters related to partial AR coefficients.
%               .maParamP0   : Initial guess of parameters related to partial MA coefficients.
%               .xParam0     : Initial guess of X coefficients.
%                          Optional. Default: order range [0 3; 0 3; 0 3] with 
%                          initial guesses of parameters picked randomly.
%       minOpt           : Minimization option: 
%                          -'ml' for Matlab local minimizer "fmincon";
%                          -'tl' for Tomlab local minimizer "ucSolve";
%                          -'tg' for Tomlab global minimizer "glbFast"' followed
%                           by Tomlab local minimizer "ucSolve";
%                          -'nag' for NAG's local minimizerE04JAF.
%                          Optional. Default: 'ml'.
%
%OUTPUT fitResults: Structure array containing fitting results. 
%                   Will have same dimensions as modelParam. 
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
%           .varCovMat : Variance-covariance matrix of ARMAX coefficients, 
%                        estimated via least squares fitting.
%           .wnVariance: Estimated variance of white noise in process.
%           .wnVector  : Structure array containing the field:
%               .observations: Estimated white noise series in corresponding 
%                              trajectory.
%           .selectCrit: Model selection criterion. Contains 3 fields:
%               .aic         : Akaike's Information Criterion.
%               .aicc        : Bias-Corrected Akaike's Information Criterion.
%               .bic         : Bayesian Information Criterion.
%           .success   : 1 if fitting executed normally, 0 otherwise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARKS Data is fitted using "armaCoefKalman" beginning with the different
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
if nargin < 2
    disp('--armaXFitKalman: You should AT LEAST supply the input and output series!');
    errFlag = 1;
    fitResults = [];
    return
end

%assign default values of optional variables
minOrder_def = 0;  %minimum AR, MA and X orders
maxOrder_def = 3;  %maximum AR, MA and X orders
minOpt_def = 'ml'; %minimization option
repeat_def = 5;    %number of times to repeat local minim. if initial guess not supplied

%check "trajOut" and turn into struct if necessary
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    clear tmp
elseif ~isfield(trajOut,'observations')
    disp('--armaXFitKalman: Please input trajOut in fields ''observations''!')
    errFlag = 1;
end

%check "trajIn" and turn into struct if necessary
if ~isstruct(trajIn)
    tmp = trajIn;
    clear trajIn
    trajIn.observations = tmp;
    clear tmp
elseif ~isfield(trajIn,'observations')
    disp('--armaXFitKalman: Please input trajIn in fields ''observations''!')
    errFlag = 1;
end

%check models to be tested
if nargin < 3 || isempty(modelParamOrOrder) %if no models to fit were input

    %randomly generate initial parameter guesses for AR, MA and X orders from 0:5
    for i=maxOrder_def:-1:minOrder_def
        for j=maxOrder_def:-1:minOrder_def
            for k=maxOrder_def:-1:minOrder_def
                modelParam(i+1,j+1,k+1).arParamP0 = 4*rand(1,i) - 2;
                modelParam(i+1,j+1,k+1).maParamP0 = 4*rand(1,j) - 2;
                modelParam(i+1,j+1,k+1).xParam0   = 4*rand(1,k) - 2;
            end
        end
    end
        
else %if models or model order were supplied

    if isstruct(modelParamOrOrder) %if initial guesses were input

        %assign model parameters
        modelParam = modelParamOrOrder;
        
        %indicate that initial guess was supplied
        suppliedIG = 1;

    else %if range of AR and MA orders was input

        %indicate that initial guess was not supplied
        suppliedIG = 0;
        
        %check size of matrix
        [nRow,nCol] = size(modelParamOrOrder);
        if nRow ~= 3 || nCol ~= 2
            disp('--armaXFitKalman: Matrix indicating range of AR, MA and X orders should be a 3x2 matrix!');
            errFlag = 1;
        else
            %generate random initial parameter guesses for the order range specified
            for i=modelParamOrOrder(1,1):modelParamOrOrder(1,2)
                for j=modelParamOrOrder(2,1):modelParamOrOrder(2,2)
                    for k=modelParamOrOrder(3,1):modelParamOrOrder(3,2)
                        
                        i1 = i - modelParamOrOrder(1,1) + 1;
                        j1 = j - modelParamOrOrder(2,1) + 1;
                        k1 = k - modelParamOrOrder(3,1) + 1;

                        for l=1:repeat_def
                            modelParam(i1,j1,k1,l).arParamP0 = 4*rand(1,i) - 2;
                            modelParam(i1,j1,k1,l).maParamP0 = 4*rand(1,j) - 2;
                            modelParam(i1,j1,k1,l).xParam0   = 4*rand(1,k) - 2;
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
    if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl') ...
            && ~strcmp(minOpt,'tg') && ~strcmp(minOpt,'nag'))
        disp('--armaXFitKalman: "minOpt" should be either "ml", "tl", "tg" or "nag"!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--armaXFitKalman: Please fix input data!');
    fitResults = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%go over all suggested models
for i=1:size(modelParam,1)
    for j=1:size(modelParam,2)
        for k=1:size(modelParam,3)

            if strcmp(minOpt,'tg') || suppliedIG %if global minimization or if initial guess was supplied

                %assign model data
                arParamP0 = modelParam(i,j,k,1).arParamP0;
                maParamP0 = modelParam(i,j,k,1).maParamP0;
                xParam0   = modelParam(i,j,k,1).xParam0;

                %estimate ARMA coeffients and white noise variance
                [arParamK,maParamK,xParamK,arParamL,maParamL,xParamL,...
                    varCovMat,wnVariance,wnVector,selectCrit,pVCompKL,...
                    pVPort,errFlag] = armaXCoefKalman(trajOut,trajIn,...
                    arParamP0,maParamP0,xParam0,[],minOpt);

            else %if local minimization and initial guess was not supplied

                neg2LL = 10^20;
                for l=1:repeat_def

                    %assign model data
                    arParamP0 = modelParam(i,j,k,l).arParamP0;
                    maParamP0 = modelParam(i,j,k,l).maParamP0;
                    xParam0   = modelParam(i,j,k,l).xParam0;

                    %estimate ARMA coeffients and white noise variance
                    [arParamK1,maParamK1,xParamK1,arParamL1,maParamL1,...
                        xParamL1,varCovMat1,wnVariance1,wnVector1,selectCrit1,...
                        pVCompKL1,pVPort1,errFlag1] = armaXCoefKalman(trajOut,...
                        trajIn,arParamP0,maParamP0,xParam0,[],minOpt);

                    if ~isempty(selectCrit1)
                        neg2LL1 = selectCrit1.aic - ...
                            2*(length([arParamP0 maParamP0 xParam0])+1);
                    end

                    if neg2LL1 < neg2LL
                        neg2LL     = neg2LL1;
                        arParamK   = arParamK1;
                        maParamK   = maParamK1;
                        xParamK    = xParamK1;
                        arParamL   = arParamL1;
                        maParamL   = maParamL1;
                        xParamL    = xParamL1;
                        varCovMat  = varCovMat1;
                        wnVariance = wnVariance1;
                        wnVector   = wnVector1;
                        selectCrit = selectCrit1;
                        pVCompKL   = pVCompKL1;
                        pVPort     = pVPort1;
                        errFlag    = errFlag1;
                    end

                end

            end

            %write output as fields in fitResults
            fitResults(i,j,k) = struct('arParamP0',arParamP0,'maParamP0',maParamP0,...
                'xParam0',xParam0,'arParamK',arParamK,'maParamK',maParamK,...
                'xParamK',xParamK,'arParamL',arParamL,'maParamL',maParamL,...
                'xParamL',xParamL,'varCovMat',varCovMat,'wnVariance',wnVariance,...
                'wnVector',wnVector,'selectCrit',selectCrit,'pVCompKL',pVCompKL,...
                'pVPort',pVPort,'success',1-errFlag);

        end %(for k=1:size(modelParam,3))
    end %(for j=1:size(modelParam,2))
end %(for i=1:size(modelParam,1))


%%%%% ~~ the end ~~ %%%%%
