function [fitResults,errFlag] = armaFitKalman(trajectories,modelParamOrder,...
    minOpt)
%ARMAFITKALMAN fits a number of ARMA models (i.e. determines their coefficients and white noise variance) to a time series which could have missing data points.
%
%SYNOPSIS [fitResults,errFlag] = armaFitKalman(trajectories,modelParamOrder,...
%    minOpt)
%
%INPUT  trajectories   : Observations of time series to be fitted. Either an 
%                        array of structures traj(1:nTraj).observations, or a
%                        2D array representing one single trajectory. 
%               .observations: 2D array of measurements and their uncertainties.
%                        Missing points should be indicated with NaN.
%       modelParamOrder: Either
%                        (1)2x2 matrix indicating range of orders of models to fit:
%                           [lowest AR order   highest AR order
%                            lowest MA order   highest MA order]
%                        Or
%                        (2)1D or 2D structure array of models to be tested, 
%                           where each entry consists of 2 elements:
%               .arParamP0   : Initial guess of parameters related to AR coef. (row vector).
%               .maParamP0   : Initial guess of parameters related to MA coef. (row vector).
%                        Optional. Default: order range [0 5; 0 5] with 
%                        initial guesses of parameters are picked randomly.
%       minOpt         : Minimization option: 
%                        -'ml' for Matlab local minimizer "fmincon";
%                        -'tl' for Tomlab local minimizer "ucSolve";
%                        -'tg' for Tomlab global minimizer "glbFast"' followed
%                         by Tomlab local minimizer "ucSolve";
%                        -'nag' for NAG's local minimizerE04JAF.
%                        Optional. Default: 'tl'.
%
%OUTPUT fitResults: Structure array containing fitting results. 
%                   Will have same dimensions as modelParam. 
%                   If orderRange is input insteadHas the fields:
%           .arParamP0   : Initial guess of AR coefficients (row vector).
%           .maParamP0   : Initial guess of MA coefficients (row vector).
%           .arParamK    : Estimated AR parameters using likelihood maximization.
%           .maParamK    : Estimated MA parameters using likelihood maximization.
%           .arParamL    : Estimated AR parameters using least squares fitting.
%           .maParamL    : Estimated MA parameters using least squares fitting.
%           .varCovMat   : Variance-covariance matrix of ARMA coefficients, 
%                          estimated via least squares fitting.
%           .wnVariance  : Estimated variance of white noise in process.
%           .wnVector    : Structure array containing the field:
%               .observations: Estimated white noise series in corresponding 
%                              trajectory.
%           .selectCrit  : Model selection criterion. Contains 3 fields:
%               .aic         : Akaike's Information Criterion.
%               .aicc        : Bias-Corrected Akaike's Information Criterion.
%               .bic         : Bayesian Information Criterion.
%           .success     : 1 if fitting executed normally, 0 otherwise.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARKS Data is fitted using "armaCoefKalman" beginning with the different
%        input models.
%
%Khuloud Jaqaman, January 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--armaFitKalman: You should AT LEAST input the trajectory to be fitted!');
    errFlag = 1;
    fitResults = [];
    return
end

%check "trajectories" and turn into struct if necessary
if ~isstruct(trajectories)
    tmp = trajectories;
    clear trajectories
    trajectories.observations = tmp;
    clear tmp
elseif ~isfield(trajectories,'observations')
    disp('--armaFitKalman: Please input the trajectories in fields ''observations''!')
    errFlag = 1;
end

%check models to be tested
if nargin < 2 || isempty(modelParamOrder) %if no models to fit were input

    %randomly generate initial parameter guesses for AR and MA orders from 0:5
    for i=0:5
        for j=0:5
            modelParam(i+1,j+1).arParamP0 = 4*rand(1,i) - 2;
            modelParam(i+1,j+1).maParamP0 = 4*rand(1,j) - 2;
        end
    end

else %if models or model order were supplied

    if isstruct(modelParamOrder) %if initial guesses were input

        modelParam = modelParamOrder;

    else %if range of AR and MA orders was input

        %check size of matrix
        [nRow,nCol] = size(modelParamOrder);
        if nRow ~= 2 || nCol ~= 2
            disp('--armaFitKalman: Matrix indicating range of AR and MA order should be a 2x2 matrix!');
            errFlag = 1;
        else
            %generate random initial parameter guesses for the order range specified
            for i=1:modelParamOrder(1,2)-modelParamOrder(1,1)+1
                for j=1:modelParamOrder(2,2)-modelParamOrder(2,1)+1
                    modelParam(i,j).arParamP0 = 4*rand(1,i) - 2;
                    modelParam(i,j).maParamP0 = 4*rand(1,j) - 2;
                end
            end
        end

    end

end %(if nargin < 2 || isempty(modelParamOrder))

%check minimization option
if nargin < 3 || isempty(minOpt) %if minimization option was not input
    minOpt = 'tl'; %assign default
else %if minimization option was input, check its value
    if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl') ...
            && ~strcmp(minOpt,'tg') && ~strcmp(minOpt,'nag'))
        disp('--armaFitKalman: "minOpt" should be either "ml", "tl", "tg" or "nag"!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--armaFitKalman: Please fix input data!');
    fitResults = [];
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%go over all suggested models
for i=1:size(modelParam,1)
    for j=1:size(modelParam,2)

        %assign model data
        arParamP0 = modelParam(i,j).arParamP0;
        maParamP0 = modelParam(i,j).maParamP0;

        %estimate ARMA coeffients and white noise variance
        [arParamK,maParamK,arParamL,maParamL,varCovMat,wnVariance,...
            wnVector,selectCrit,pVCompKL,pVPort,errFlag] = ...
            armaCoefKalman(trajectories,arParamP0,maParamP0,[],minOpt);

        %write output as fields in fitResults
        fitResults(i,j) = struct('arParamP0',arParamP0,'maParamP0',...
            maParamP0,'arParamK',arParamK,'maParamK',maParamK,...
            'arParamL',arParamL,'maParamL',maParamL,'varCovMat',...
            varCovMat,'wnVariance',wnVariance,'wnVector',wnVector,...
            'selectCrit',selectCrit,'pVCompKL',pVCompKL,'pVPort',pVPort,...
            'success',1-errFlag);

    end %(for j=1:size(modelParam,2))
end %(for i=1:size(modelParam,1))


%%%%% ~~ the end ~~ %%%%%
