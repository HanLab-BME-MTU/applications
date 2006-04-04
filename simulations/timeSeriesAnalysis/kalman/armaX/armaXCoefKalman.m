function [arParamK,maParamK,xParamK,arParamL,maParamL,xParamL,varCovMatL,...
    varCovMatF,wnVariance,wnVector,selectCrit,pVCompKL,pVPort,errFlag] ...
    = armaXCoefKalman(trajOut,trajIn,arParamP0,maParamP0,xParam0,...
    constParam,minOpt)
%ARMAXCOEFKALMAN fits an ARMAX model to a time series which could depend on another series.
%
%SYNOPSIS [arParamK,maParamK,xParamK,arParamL,maParamL,xParamL,varCovMatLS,...
%    varCovMatFI,wnVariance,wnVector,selectCrit,pVCompKL,pVPort,errFlag] ...
%    = armaXCoefKalman(trajOut,trajIn,arParamP0,maParamP0,xParam0,...
%    constParam,minOpt)
%
%INPUT
%   Mandatory
%       trajOut   : Observations of output time series to be fitted. Either an
%                   array of structures trajOut(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Missing points should be indicated with NaN.
%       trajIn    : Observations of input time series. Either an
%                   array of structures trajIn(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Must not have any missing points.
%                   Enter as [] if there is no input series.
%       arParamP0 : Initial guess of parameters determining the AR coef. (row vector).
%                   They are related to the partial AR coef by the
%                   equation: partial AR coef. =
%                   (1-exp(arParamP0))/(1+exp(arParamP0)).
%                   Enter as [] if there is no AR part in model.
%   Optional
%       maParamP0 : Initial guess of parameters determining the MA coef. (row vector).
%                   They are related to the partial MA coef by the
%                   equation: partial MA coef. =
%                   (1-exp(maParamP0))/(1+exp(maParamP0)).
%                   Enter as [] if there is no MA part in model.
%                   Default: [].
%       xParam0   : Initial guess of coefficients of dependence on input.
%                   Enter as [] if there is no dependence on input.
%                   Default: [].
%       constParam: Set of constrained parameters. Constains 2 fields:
%           .ar          : 2D array. 1st column is AR parameter number and
%                          2nd column is parameter value. No need to input if
%                          there are no constraints on AR parameters.
%           .ma          : 2D array. 1st column is MA parameter number and
%                          2nd column is parameter value.No need to input if
%                          there are no constraints on MA parameters.
%                   Default: []
%       minOpt    : Minimization option:
%                   -'ml' for Matlab local minimizer "fmincon";
%                   -'tl' for Tomlab local minimizer "ucSolve";
%                   -'tg' for Tomlab global minimizer "glbFast"' followed
%                     by Tomlab local minimizer "ucSolve";
%                   -'nag' for NAG's local minimizerE04JAF.
%                   Default: 'ml'
%
%OUTPUT arParamK  : Estimated AR coefficients (1st row) and parameters related
%                   to partial AR coefficients (2nd row) using likelihood maximization.
%       maParamK  : Estimated MA coefficients (1st row) and parameters related
%                   to partial MA coefficients (2nd row) using likelihood maximization.
%       xParamK   : Estimated X coefficients using likelihood maximization,
%                   indicating dependence on input.
%       arParamL  : Estimated AR coefficients using least squares fitting.
%       maParamL  : Estimated MA coefficients using least squares fitting.
%       xParamL   : Estimated X coefficients using least squares fitting.
%       varCovMatL: Variance-covariance matrix of ARMAX coefficients,
%                   estimated via least squares fitting:
%           .cofactorMat  : Cofactor matrix.
%           .posterioriVar: A posteriori estimate of residuals' variance.
%       varCovMatF: Variance-covariance matrix of ARMAX coefficients,
%                   estimated via the Fisher Information Matrix.
%       wnVariance: Estimated variance of white noise in process.
%       wnVector  : Structure array containing the field:
%           .observations : Estimated white noise series in corresponding
%                           trajectory.
%       selectCrit: Model selection criterion. Contains 3 fields:
%            .aic         : Akaike's Information Criterion.
%            .aicc        : Bias-Corrected Akaike's Information Criterion.
%            .bic         : Bayesian Information Criterion.
%       pVCompKL  : P-value for comparing ARMA coefficients obtained by
%                   Kalman filter and maximum likelihood estimation with
%                   those obtained using least squares fitting.
%       pVPort    : P-value of the portmanteau test of the residuals.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%REMARKS The Kalman filter & likelihood maximization algorithm implemented
%        here is an ARMAX generalization of the algorithm presented in
%        R. H. Jones, "Maximum Likelihood Fitting of ARMA Models to Time
%        Series with Missing Observations," Technometrics 22: 389-395
%        (1980). The ARMAX generalization is based on material in C. K.
%        Chui and G. Chen, "Kalman Filtering with Real-Time Applications,"
%        3rd ed. (1999), Ch. 2.
%
%        Equation numbers used are those in the Jones paper.
%
%        One main difference between this algorithm and that in the paper
%        is that I do not estimate the observational error variance,
%        but use that obtained via image analysis in the case of
%        experimental data or simulated in the case of simulated data.
%
%        Another main difference is that I use more than one trajectory to
%        calculate one set of coefficients for the condition they
%        represent. Thus Eqs. 3.14 and 3.15 in the paper are modified
%        accordingly.
%
%        Then the problem is re-formulated as a least squares problem in
%        order to obtain the variance-covariance matrix of the estimated
%        coefficients. The output is regressed onto its past values, the
%        white noise series estimated above, and the input series.
%
%        CONSTRAINED MINIMIZATION MUST BE UPDATED. DON'T USE!
%        CURRENTLY I SUBTRACT THE MEAN. THINK ABOUT IT!
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
errFlag    =  0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether all mandatory input variables were supplied
if nargin < 3
    disp('--armaXCoefKalman: The function requires at least 3 input variables!');
    errFlag  = 1;
    return
end

%assign default values of optional variables
maParamP0_def  = [];
xParam0_def    = [];
constParam_def = [];
minOpt_def     = 'ml';

%check "trajOut" and turn it into struct if necessary
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    clear tmp
elseif ~isfield(trajOut,'observations')
    disp('--armaXCoefKalman: Please input trajOut in fields ''observations''!')
    errFlag = 1;
end

%get number of trajectories supplied
numTraj = length(trajOut);

%shift trajectories so that each trajectory's mean = 0
for i=1:numTraj
    traj = trajOut(i).observations(:,1);
    traj = traj - nanmean(traj);
    trajOut(i).observations(:,1) = traj;
end

%add column for observational error of output if not provided
trajOriginal = trajOut; %needed for least squares
for i=1:numTraj
    traj = trajOut(i).observations;
    [trajLength,nCol] = size(traj);
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            traj = [traj zeros(trajLength,1)];
        else
            disp('--armaXCoefKalman: "trajOut.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
    trajOut(i).observations = traj;
end

%check "trajIn", turn into struct if necessary and add column for
%observational error if not provided
if isempty(trajIn) %if there is no input

    for i=1:numTraj
        trajIn(i).observations = [];
    end

else %if there is an input series

    if ~isstruct(trajIn) %turn into struct
        tmp = trajIn;
        clear trajIn
        trajIn.observations = tmp;
        clear tmp
    elseif ~isfield(trajIn,'observations')
        disp('--armaXCoefKalman: Please input trajIn in fields ''observations''!')
        errFlag = 1;
    end

    for i=1:numTraj %check for observational error
        traj = trajIn(i).observations;
        [trajLength,nCol] = size(traj);
        if nCol ~= 2
            if nCol == 1 %if no error is supplied, assign it the value 0
                traj = [traj zeros(trajLength,1)];
            else
                disp('--armaXCoefKalman: "trajIn.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
                errFlag = 1;
            end
        end
        trajIn(i).observations = traj;
    end

end

%get arOrder
[nRow,arOrder] = size(arParamP0);
if ~isempty(arParamP0)
    if nRow ~= 1
        disp('--armaXCoefKalman: "arParamP0" should be a row vector!');
        errFlag = 1;
    end
end

%get maOrder
if nargin < 4 || isempty(maParamP0)
    maOrder = 0;
    maParamP0 = maParamP0_def;
else
    [nRow,maOrder] = size(maParamP0);
    if nRow ~= 1
        disp('--armaXCoefKalman: "maParamP0" should be a row vector!');
        errFlag = 1;
    end
end

%get xOrder
if nargin < 5 || isempty(xParam0)
    xOrder = -1;
    xParam0 = xParam0_def;
else
    if isempty(trajIn)
        disp('--armaXCoefKalman: "xParam0" is supplied but there is no "trajIn"!');
        errFlag = 1;
    end
    [nRow,xOrder] = size(xParam0);
    if nRow ~= 1
        disp('--armaXCoefKalman: "xParam0" should be a row vector!');
        errFlag = 1;
    end
    xOrder = xOrder - 1;
end

%check parameter constraints
if nargin < 6 || isempty(constParam) %if no constraints were entered
    constParam = constParam_def;
else
    if isfield(constParam,'ar')
        nCol = size(constParam.ar,2);
        if nCol ~= 2
            disp('--armaXCoefKalman: constParam.ar should have 2 columns!');
            errFlag = 1;
        else
            if min(constParam.ar(:,1)) < 1 || max(constParam.ar(:,1)) > arOrder
                disp('--armaXCoefKalman: Wrong AR parameter numbers in constraint!');
                errFlag = 1;
            end
        end
    else
        constParam.ar = zeros(0,2);
    end
    if isfield(constParam,'ma')
        nCol = size(constParam.ma,2);
        if nCol ~= 2
            disp('--armaXCoefKalman: constParam.ma should have 2 columns!');
            errFlag = 1;
        else
            if min(constParam.ma(:,1)) < 1 || max(constParam.ma(:,1)) > maOrder
                disp('--armaXCoefKalman: Wrong MA parameter numbers in constraint!');
                errFlag = 1;
            end
        end
    else
        constParam.ma = zeros(0,2);
    end
end %(nargin < 6 || isempty(constParam) ... else ...)

%check whether minOpt has one of the required values
if nargin < 7 || isempty(minOpt)
    minOpt = minOpt_def;
else
    if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl') ...
            && ~strcmp(minOpt,'tg') && ~strcmp(minOpt,'nag'))
        disp('--armaXCoefKalman: "minOpt" should be either "ml", "tl", "tg" or ''nag''!');
        errFlag = 1;
    end
end

%exit if there are problems in input data
if errFlag
    disp('--armaXCoefKalman: Please fix input data!');
    return
end

%obtain number of available observations, per trajectory and total
numAvail = zeros(1,numTraj);
for i=1:numTraj
    traj = trajOut(i).observations(:,1);
    numAvail(i) = length(find(~isnan(traj)));
end
totAvail = sum(numAvail);

%get an initial estimate of the white noise variance
tmp = vertcat(trajOut.observations);
wnVariance0 = nanvar(tmp(:,1)); %variance from "previous" iteration
wnVariance = 0.8*wnVariance0; %variance from "current" iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maximum likelihood estimation of model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%since in Jones' algorithm the variance of the observational error must be
%divided by the unknown variance of the white noise, the minimization is
%done iteratively until the white noise variance does not change
%significantly between iterations.

%assign initial guess of parameters to paramT
paramT = [arParamP0 maParamP0 xParam0];

%while the variance changes by more than 5% from one iteration to the next
while abs(wnVariance-wnVariance0)/wnVariance0 > 0.05

    %update wnVariance0
    wnVariance0 = wnVariance;

    %divide the observational error by the standard deviation of the white
    %noise in all trajectories
    trajOut2 = trajOut;
    wnStd = sqrt(wnVariance);
    for i=1:numTraj
        trajOut2(i).observations(:,2) = trajOut2(i).observations(:,2)/wnStd;
    end

    %get initial guess of parameters
    param0 = paramT;
    
    %if there are coefficients to determine
    if arOrder + maOrder ~= 0 || ~isempty(xParam0)

        %decide which minimization algorithm to use
        switch minOpt

            case 'ml' %local minimization using Matlab's fmincon

                %define optimization options.
                options = optimset('Display','final','MaxFunEvals',10000,...
                    'MaxIter',1000);

                %define structure containing additional parameters
                %note that it is written in Tomlab notation for convenience
                prob.user.arOrder = arOrder;
                prob.user.maOrder = maOrder;
                prob.user.trajOut = trajOut2;
                prob.user.trajIn  = trajIn;
                prob.user.numAvail = totAvail;

                try

                    %assign lower and upper bounds of variables
                    boundLow  = [-10*ones(1,arOrder+maOrder) -2*ones(1,xOrder+1)];
                    boundHigh = [10*ones(1,arOrder+maOrder) 2*ones(1,xOrder+1)];

                    %                     %minimize -2ln(likelihood) using fmincon
                    %                     [params,fval,exitFlag] = fmincon(@neg2LnLikelihoodX,param0,...
                    %                         [],[],[],[],boundLow,boundHigh,[],options,prob);

                    %minimize -2ln(likelihood) using fminunc
                    [params,fval,exitFlag] = fminunc(@neg2LnLikelihoodX,param0,...
                        options,prob);

                catch

                    exitFlag = 0;

                end

                %proceed if minimization was successful
                if exitFlag > 0 %successful
                    proceed = 1;
                else %not successful
                    proceed = 0;
                end

            case 'tl' %local minimization using Tomlab's ucSolve

                %                 if isempty(constParam) %if there are no constraints

                prob = conAssign('neg2LnLikelihoodX',[],[],[],...
                    [-10*ones(1,arOrder+maOrder) -2*ones(1,xOrder+1)],...
                    [10*ones(1,arOrder+maOrder) 2*ones(1,xOrder+1)],...
                    'locMinNegLik',param0);
                prob.PriLevOpt = 0;
                prob.user.arOrder = arOrder;
                prob.user.maOrder = maOrder;
                prob.user.trajOut = trajOut2;
                prob.user.trajIn  = trajIn;
                prob.user.numAvail = totAvail;

                %minimize -2ln(likelihood) using Tomlab's ucSolve
                result = tomRun('ucSolve',prob,[],2);

                %                 else %if there are constraints
                %
                %                     %define local minimizaton problem with constraints
                %                     prob = conAssign('neg2LnLikelihood',[],[],[],...
                %                         -9*ones(1,arOrder+maOrder),9*ones(1,arOrder+maOrder),...
                %                         'locMinNegLik',param0,[],[],[],[],[],'minKalmanConstraint',...
                %                         [],[],[],[constParam.ar(:,2);constParam.ma(:,2)],...
                %                         [constParam.ar(:,2);constParam.ma(:,2)]);
                %                     prob.PriLevOpt = 0;
                %                     prob.user.arOrder = arOrder;
                %                     prob.user.trajectories = trajectories;
                %                     prob.user.numAvail = totAvail;
                %                     prob.user.constParam.ar = constParam.ar(:,1);
                %                     prob.user.constParam.ma = constParam.ma(:,1);
                %
                %                     %minimize -2ln(likelihood) using Tomlab's conSolve
                %                     result = tomRun('conSolve',prob,[],2);
                %
                %                 end

                %proceed if minimization was successful
                if result.ExitFlag == 0
                    params = result.x_k';
                    proceed = 1;
                else
                    proceed = 0;
                end

        end %(switch minOpt)

    else %i.e. if there are no AR, MA or X coefficients to determine

        params = [];
        proceed = 1;

    end %(if arOrder + maOrder + xOrder ~= 0 ... else ...)

    %if minimization was successful
    if proceed

        %assign parameter values to paramT to use them as starting guess in
        %next iteration
        paramT = params;
        
        %assign parameters obtained through minimization
        arParamP = params(1:arOrder);
        maParamP = params(arOrder+1:arOrder+maOrder);
        xParamK  = params(arOrder+maOrder+1:end);

        %get AR and MA coefficients from the partial AR and MA coefficients, respectively
        if ~isempty(arParamP)
            [arParamK,errFlag] = levinsonDurbinAR(arParamP);
        else
            arParamK = [];
        end
        if ~isempty(maParamP)
            [maParamK,errFlag] = levinsonDurbinMA(maParamP);
        else
            maParamK = [];
        end

        %obtain likelihood, white noise sequence and white noise variance
        sum1 = 0;
        sum2 = 0;
        wnVarianceSamp = zeros(numTraj,1);
        for i = 1:numTraj %for all trajectories

            %get the innovations, their variances and the estimated white noise series
            %using Kalman prediction and filtering
            [innovation,innovationVar,wnVector(i).observations,dummy1,...
                dummy2,errFlag] = armaXKalmanInnov(trajOut2(i).observations,...
                trajIn(i).observations,arParamK,maParamK,xParamK);
            if errFlag
                disp('--armaXCoefKalman: "armaXKalmanInnov" did not function properly!');
                return
            end

            %calculate variance of white noise in current trajectory
            wnVarianceSamp(i) = nanmean(innovation.^2./innovationVar);

            %1st sum in Eq. 3.15
            sum1 = sum1 + nansum(log(innovationVar));
            %2nd sum in Eq. 3.15
            sum2 = sum2 + nansum(innovation.^2./innovationVar);

        end %(for i = 1:numTraj)

        %calculate -2ln(likelihood) (Eq. 3.15)
        neg2LnLikelihoodV = sum1 + totAvail*log(sum2);

        %calculate mean white noise variance of all trajectories (Eq. 3.14)
        wnVariance = (numAvail*wnVarianceSamp)/totAvail;

        %get number of parameters estimated: arOrder AR coefficients, maOrder MA
        %coefficients, xOrder+1``` X coefficients, and white noise variance
        numParam = arOrder + maOrder + xOrder + 2;

        %evaluate Akaike's Information Criterion
        selectCrit.aic = neg2LnLikelihoodV + 2*numParam;

        %evaluate the bias-corrected Akaike's Information Criterion
        selectCrit.aicc = neg2LnLikelihoodV + 2*numParam*totAvail/(totAvail-numParam-1);

        %evaluate the Bayesian Information Criterion
        selectCrit.bic = neg2LnLikelihoodV + log(totAvail)*numParam;

        %put partial coefficients in 2nd row of coefficients matrix
        arParamK(2,:) = arParamP;
        maParamK(2,:) = maParamP;

    else %if minimization was not successful

        errFlag = 1;
        return

    end %(if proceed)

end %(while abs(wnVariance-wnVariance0)/wnVariance0 > 0.05)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Least squares fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reformulate the problem as a least squares fitting and obtain the
%variance-covariance matrix of the estimated ARMA coefficients
[varCovMatL,arParamL,maParamL,xParamL,errFlag] = armaXLeastSquares(...
    trajOriginal,trajIn,wnVector,arOrder,maOrder,xOrder,constParam,...
    wnVariance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimation of variance-covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty([arParamK(1,:) maParamK(1,:) xParamK(1,:)])

    %calculate the model's Fisher information matrix
    [fishInfoMat,errFlag] = armaXFisherInfoMatrix(trajOut,trajIn,...
        arParamK(1,:),maParamK(1,:),xParamK,wnVariance);

    %get the variance-covariance matrix of the ARMA coefficients
    varCovMatF = inv(fishInfoMat/totAvail)/totAvail;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if parameters found through least squares fitting are statistically
%equivalent to those found through maximum likelihood estimation with
%Kalman filtering. If they are not equivalent, then results cannot be
%trusted and model is discarded

%prepare input
armaXCoef1.arParam = arParamK(1,:);
armaXCoef1.maParam = maParamK(1,:);
armaXCoef1.xParam = xParamK;
armaXCoef2.arParam = arParamL;
armaXCoef2.maParam = maParamL;
armaXCoef2.xParam = xParamL;

%compare parameters
[H,pVCompKL,errFlag] = armaXCoefComp(armaXCoef1,armaXCoef2,varCovMatL,varCovMatL,...
    'global');
if errFlag
    pVCompKL = 0;
end

%report failure of fit and do not consider results if coefficients are significantly different
if H == 1
    %     disp('--armaXCoefKalman: Discrepency between least squares and maximum likelihood!')
    errFlag = 1;
end

%use the portmanteau test to check whether residuals are white noise.
[H,pVPort,errFlag] = portmanteau(wnVector,10,0.01);

%report failure of fit and do not consider results if residuals are not white noise
if H == 1
    %     disp('--armaXCoefKalman: Residuals did not pass portmanteau test!')
    errFlag = 1;
end


%%%%% ~~ the end ~~ %%%%%


%THESE ARE THE OTHER ALGORITHMS THAT IN PRINCIPLE CAN BE USED FOR
%OPTIMIZATION!
%
%         case 'tg' %global minimization using Tomlab's glbFast and ucSolve
%
%             %initial parameter values
%             param0 = [arParamP0 maParamP0];
%
%             if ~isempty(param0)
%
%                 %define global minimization problem
%                 prob = glcAssign('neg2LnLikelihood',-9*ones(1,arOrder+maOrder),...
%                     9*ones(1,arOrder+maOrder),'globMinNegLik');
%                 prob.PriLevOpt = 2;
%                 %                 prob.optParam.MaxFunc = 15000;
%                 prob.user.arOrder = arOrder;
%                 prob.user.trajectories = trajectories;
%                 prob.user.numAvail = totAvail;
%
%                 %find global minimum of -2ln(likelihood) using Tomlab's glbFast
%                 result = tomRun('glbFast',prob,[],2);
%
%             else
%
%                 %avoid global minimization if ARMA(0,0) is being "optimized"!
%                 result.ExitFlag = 1;
%
%             end
%
%             if result.ExitFlag == 0 %if global minimization was successful
%                 paramI = result.x_k(:,1)'; %use its results as initial guess for local minimization
%             else %if global minimization failed
%                 paramI = param0; %use user's initial guess as initial guess for local minimization
%             end
%
%             %define local minimizaton problem
%             prob = conAssign('neg2LnLikelihood',[],[],[],...
%                 -10*ones(1,arOrder+maOrder),10*ones(1,arOrder+maOrder),...
%                 'locMinNegLik',paramI);
%             prob.PriLevOpt = 0;
%             prob.user.arOrder = arOrder;
%             prob.user.trajectories = trajectories;
%             prob.user.numAvail = totAvail;
%
%             %refine minimum using ucSolve
%             result = tomRun('ucSolve',prob,[],2);
%
%             %proceed if minimization was successful
%             if result.ExitFlag == 0
%                 params = result.x_k';
%                 proceed = 1;
%             else
%                 proceed = 0;
%             end
%
%         case 'nag' %local minimization using NAG's E04JAF
%
%             %define structure containing parameters required for function
%             %evaluation; they are written in Tomlab notation for convenience
%             prob.user.arOrder = arOrder;
%             prob.user.trajectories = trajectories;
%             prob.user.numAvail = totAvail;
%
%             %save "prob" in file "funct1Input" so that funct1 loads the
%             %variables when called.
%             save('funct1Input','prob');
%
%             %initial parameter values
%             params = [arParamP0 maParamP0];
%
%             [params,fval,lowerB,upperB,exitFlag] = ...
%                 e04jaf(params,[-0.99*ones(1,arOrder+maOrder)],...
%                 [0.99*ones(1,arOrder+maOrder)],0);
%
%             %proceed if minimization was successful
%             if (exitFlag == 0 || exitFlag == 5 || exitFlag == 6)
%                 params = params';
%                 proceed = 1;
%             else
%                 proceed = 0;
%             end
