function [arParamK,maParamK,arParamL,maParamL,varCovMat,wnVariance,...
    wnVector,selectCrit,pVCompKL,pVPort,errFlag] = armaCoefKalman(...
    trajectories,arParamP0,maParamP0,constParam,minOpt)
%ARMACOEFKALMAN fits an ARMA(p,q) model to a time series which could have missing data points.
%
%SYNOPSIS [arParamK,maParamK,arParamL,maParamL,varCovMat,wnVariance,...
%    wnVector,selectCrit,pVCompKL,pVPort,errFlag] = armaCoefKalman(...
%    trajectories,arParamP0,maParamP0,constParam,minOpt)
%
%INPUT  trajectories: Observations of time series to be fitted. Either an 
%                     array of structures traj(1:nTraj).observations, or a
%                     2D array representing one single trajectory. 
%           .observations: 2D array of measurements and their uncertainties.
%                     Missing points should be indicated with NaN.
%       arParamP0   : Initial guess of parameters determining the AR coef. (row vector).
%                     They are related to the partial AR coef by the
%                     equation: partial AR coef. =
%                     (1-exp(arParamP0))/(1+exp(arParamP0))
%       maParamP0   : Initial guess of parameters derermining the MA coef. (row vector).
%                     They are related to the partial MA coef by the
%                     equation: partial MA coef. =
%                     (1-exp(maParamP0))/(1+exp(maParamP0))
%       constParam  : Set of constrained parameters. Constains 2 fields:
%           .ar     : 2D array. 1st column is AR parameter number and
%                     2nd column is parameter value. No need to input if
%                     there are no constraints on AR parameters.
%           .ma     : 2D array. 1st column is MA parameter number and
%                     2nd column is parameter value.No need to input if
%                     there are no constraints on MA parameters.
%                     Optional. Default: 0
%       minOpt      : Optional. Minimization option: 
%                     -'ml' for Matlab local minimizer "fmincon";
%                     -'tl' for Tomlab local minimizer "ucSolve";
%                     -'tg' for Tomlab global minimizer "glbFast"' followed
%                       by Tomlab local minimizer "ucSolve";
%                     -'nag' for NAG's local minimizerE04JAF.
%                     Default: 'tl'
%
%OUTPUT arParamK    : Estimated AR parameters using likelihood maximization.
%       maParamK    : Estimated MA parameters using likelihood maximization.
%       arParamL    : Estimated AR parameters using least squares fitting.
%       maParamL    : Estimated MA parameters using least squares fitting.
%       varCovMat   : Variance-covariance matrix of ARMA coefficients, 
%                     estimated via least squares fitting:
%           .cofactorMat  : Cofactor matrix.
%           .posterioriVar: A posteriori estimate of residuals' variance.
%       wnVariance  : Estimated variance of white noise in process.
%       wnVector    : Structure array containing the field:
%           .observations : Estimated white noise series in corresponding 
%                           trajectory.
%       selectCrit  : Model selection criterion. Contains 3 fields:
%            .aic         : Akaike's Information Criterion.
%            .aicc        : Bias-Corrected Akaike's Information Criterion.
%            .bic         : Bayesian Information Criterion.
%       pVCompKL    : P-value for comparing ARMA coefficients obtained by
%                     Kalman filter and maximum likelihood estimation with
%                     those obtained using least squares fitting.
%       pVPort      : P-value of the portmanteau test of the residuals.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS The Kalman filter & likelihood maximization algorithm implemented 
%        here is that presented in R. H. Jones, "Maximum Likelihood Fitting
%        of ARMA Models to Time Series with Missing Observations",
%        Technometrics 22: 389-395 (1980). All equation numbers used here 
%        are those in that paper. The main difference is that I do not 
%        estimate the observational error variance, but use that obtained 
%        from experimental data or simulated trajectories (and is thus 
%        time-dependent). 
%
%        After the ARMA coefficients and white noise in process are
%        estimated using the above algorithm, the latter is used to do a
%        least squares fitting of an ARMA model to the data. From this, one
%        gets another estimate of the ARMA cofficients, as well as the
%        variance-covariance matrix of the estimated coefficients. Note
%        that the values of the ARMA coeffients obtained via the two
%        methods should agree with each other.
%
%        Finally, trajectories are shifted to get zero mean before analysis
%        is done.
%
%        THE 'ml', 'tg' AND 'nag' MINIMIZATION OPTIONS SHOULD BE MODIFIED
%        TO ACCOUNT FOR CONSTRAINTS AND FOR NEW MINIMIZATION PARAMETERS! 
%        ONLY 'tl' DOES RIGHT NOW!
%
%Khuloud Jaqaman, July 2004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arParamK = [];
maParamK = [];
arParamL = [];
maParamL = [];
varCovMat = [];
wnVariance = [];
wnVector = [];
selectCrit = [];
pVCompKL = [];
pVPort = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 3
    disp('--armaCoefKalman: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check trajectory and turn it into struct if necessary
if ~isstruct(trajectories)
    tmp = trajectories;
    clear trajectories
    trajectories.observations = tmp;
    clear tmp
elseif ~isfield(trajectories,'observations')
    disp('--armaCoefKalman: Please input the trajectories in fields ''observations''!')
    errFlag = 1;
    return
end

%get number of trajectories and add column for observational error if not input
trajOriginal = trajectories;
for i=1:length(trajectories);
    traj = trajectories(i).observations;
    [trajLength,nCol] = size(traj);
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            traj = [traj zeros(trajLength,1)];
        else
            disp('--armaCoefKalman: "trajectories.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
            errFlag = 1;
        end
    end
    trajectories(i).observations = traj;
end

%get arOrder
[nRow,arOrder] = size(arParamP0);
if ~isempty(arParamP0)
    if nRow ~= 1
        disp('--armaCoefKalman: "arParamP0" should be a row vector!');
        errFlag = 1;
    end
end

%get maOrder
[nRow,maOrder] = size(maParamP0);
if ~isempty(maParamP0)
    if nRow ~= 1
        disp('--armaCoefKalman: "maParamP0" should be a row vector!');
        errFlag = 1;
    end
end

%check parameter constraints
if nargin < 4 || isempty(constParam) %if no constraints were entered
    constParam = [];
else
    if isfield(constParam,'ar')
        [nRow,nCol] = size(constParam.ar);
        if nCol ~= 2
            disp('--armaCoefKalman: constParam.ar should have 2 columns!');
        else
            if min(constParam.ar(:,1)) < 1 || max(constParam.ar(:,1)) > arOrder
                disp('--armaCoefKalman: Wrong AR parameter numbers in constraint!');
            end
        end
    else
        constParam.ar = zeros(0,2);
    end
    if isfield(constParam,'ma')
        [nRow,nCol] = size(constParam.ma);
        if nCol ~= 2
            disp('--armaCoefKalman: constParam.ma should have 2 columns!');
        else
            if min(constParam.ma(:,1)) < 1 || max(constParam.ma(:,1)) > maOrder
                disp('--armaCoefKalman: Wrong MA parameter numbers in constraint!');
            end
        end
    else
        constParam.ma = zeros(0,2);
    end
end

%check whether minOpt has one of the required values
if nargin < 5 || isempty(minOpt)
    minOpt = 'tl';
else
    if (~strcmp(minOpt,'ml') && ~strcmp(minOpt,'tl') ...
            && ~strcmp(minOpt,'tg') && ~strcmp(minOpt,'nag'))
        disp('--armaCoefKalman: "minOpt" should be either "ml", "tl" or "tg"!');
        errFlag = 1;
    end
end %(if nargin < 5 || isempty(minOpt))

%exit if there are problems in input data
if errFlag
    disp('--armaCoefKalman: Please fix input data!');
    return
end

%obtain number of available observations and shift trajectories to get zero mean
for i=1:length(trajectories)
    
    traj = trajectories(i).observations;
    numAvail(i) = length(find(~isnan(traj(:,1)))); %get # of available points
    traj(:,1) = traj(:,1) - nanmean(traj(:,1)); %shift trajectory
    trajectories(i).observations = traj;

end
totAvail = sum(numAvail); %calculate total number of available points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maximum likelihood estimation of model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if arOrder + maOrder ~= 0

    switch minOpt

        case 'ml' %local minimization using Matlab's fmincon

            %define optimization options.
            options = optimset('Display','final','DiffMaxChange',1e-3,...
                'DiffMinChange',1e-8,'TolFun',1e-4,'TolX',1e-4,...
                'maxFunEvals',4000,'maxIter',3000);

            %define structure containing additional parameters
            %note that it i0s written in Tomlab notation for convenience
            prob.user.arOrder = arOrder;
            prob.user.trajectories = trajectories;
            prob.user.numAvail = totAvail;

            %initial parameter values
            param0 = [arParamP0 maParamP0];

            %minimize -2ln(likelihood) using fmincon
            [params,fval,exitFlag] = fmincon(@neg2LnLikelihood,param0,[],...
                [],[],[],-0.99*ones(1,arOrder+maOrder),0.99*ones(1,...
                arOrder+maOrder),[],options,prob);

            %proceed if minimization was successful
            if exitFlag > 0
                proceed = 1;
            else
                proceed = 0;
            end

        case 'tl' %local minimization using Tomlab's ucSolve

            %initial parameter values
            param0 = [arParamP0 maParamP0];

            if isempty(constParam) %if there are no constraints

                prob = conAssign('neg2LnLikelihood',[],[],[],[],[],...
                    'locMinNegLik',param0);
                prob.PriLevOpt = 0;
                prob.user.arOrder = arOrder;
                prob.user.trajectories = trajectories;
                prob.user.numAvail = totAvail;

                %minimize -2ln(likelihood) using Tomlab's ucSolve
                result = tomRun('ucSolve',prob,[],2);

            else %if there are constraints

                %define local minimizaton problem with constraints
                prob = conAssign('neg2LnLikelihood',[],[],[],[],[],...
                    'locMinNegLik',param0,[],[],[],[],[],'minKalmanConstraint',...
                    [],[],[],[constParam.ar(:,2);constParam.ma(:,2)],...
                    [constParam.ar(:,2);constParam.ma(:,2)]);
                prob.PriLevOpt = 0;
                prob.user.arOrder = arOrder;
                prob.user.trajectories = trajectories;
                prob.user.numAvail = totAvail;
                prob.user.constParam.ar = constParam.ar(:,1);
                prob.user.constParam.ma = constParam.ma(:,1);

                %minimize -2ln(likelihood) using Tomlab's conSolve
                result = tomRun('conSolve',prob,[],2);

            end
                
            %proceed if minimization was successful
            if result.ExitFlag == 0
                params = result.x_k';
                proceed = 1;
            else
                proceed = 0;
            end

        case 'tg' %global minimization using Tomlab's glbFast and ucSolve

            %initial parameter values
            param0 = [arParamP0 maParamP0];

            if ~isempty(param0)

                %define global minimization problem
                prob = glcAssign('neg2LnLikelihood',-10*ones(1,arOrder+maOrder),...
                    10*ones(1,arOrder+maOrder),'globMinNegLik');
                prob.PriLevOpt = 2;
                prob.optParam.MaxFunc = 15000;
                prob.user.arOrder = arOrder;
                prob.user.trajectories = trajectories;
                prob.user.numAvail = totAvail;

                %find global minimum of -2ln(likelihood) using Tomlab's glbFast
                result = tomRun('glbFast',prob,[],2);

            else

                %avoid global minimization if ARMA(0,0) is being "optimized"!
                result.ExitFlag = 1;

            end

            if result.ExitFlag == 0 %if global minimization was successful
                paramI = result.x_k(:,1)'; %use its results as initial guess for local minimization
            else %if global minimization failed
                paramI = param0; %use user's initial guess as initial guess for local minimization
            end

            %define local minimizaton problem
            prob = conAssign('neg2LnLikelihood',[],[],[],[],[],...
                'locMinNegLik',paramI);
            prob.PriLevOpt = 0;
            prob.user.arOrder = arOrder;
            prob.user.trajectories = trajectories;
            prob.user.numAvail = totAvail;

            %refine minimum using ucSolve
            result = tomRun('ucSolve',prob,[],2);

            %proceed if minimization was successful
            if result.ExitFlag == 0
                params = result.x_k';
                proceed = 1;
            else
                proceed = 0;
            end

        case 'nag' %local minimization using NAG's E04JAF

            %define structure containing parameters required for function
            %evaluation; they are written in Tomlab notation for convenience
            prob.user.arOrder = arOrder;
            prob.user.trajectories = trajectories;
            prob.user.numAvail = totAvail;

            %save "prob" in file "funct1Input" so that funct1 loads the
            %variables when called.
            save('funct1Input','prob');

            %initial parameter values
            params = [arParamP0 maParamP0];

            [params,fval,lowerB,upperB,exitFlag] = ...
                e04jaf(params,[-0.99*ones(1,arOrder+maOrder)],...
                [0.99*ones(1,arOrder+maOrder)],0);

            %proceed if minimization was successful
            if (exitFlag == 0 || exitFlag == 5 || exitFlag == 6)
                params = params';
                proceed = 1;
            else
                proceed = 0;
            end

    end %(switch minOpt)

else %if arOrder+maOrder=0
    
    params = [];
    proceed = 1;
    
end %(if arOrder + maOrder ~= 0)

%if minimization was successful
if proceed

    %assign parameters obtained through minimization
    arParamP = params(1:arOrder);
    maParamP = params(arOrder+1:end);

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

    %check for causality and invertibility of estimated model
    %these two criteria should be taken care of by using the partial AR and
    %MA coefficients, but I do this check just in case something goes wrong
    r = abs(roots([-arParamK(end:-1:1) 1]));
    if ~isempty(find(r<=1))
        disp('--armaCoefKalman: Predicted model not causal!');
        errFlag = 1;
    end
    r = abs(roots([maParamK(end:-1:1) 1]));
    if ~isempty(find(r<=1))
        disp('--armaCoefKalman: Predicted model not invertible!');
        errFlag = 1;
    end

    %obtain likelihood, white noise sequence and white noise variance
    sum1 = 0;
    sum2 = 0;
    for i = 1:length(trajectories)

        %get the innovations, their variances and the estimated white noise series
        %using Kalman prediction and filtering
        [innovation,innovationVar,wnVector(i).observations,errFlag] = ...
            armaKalmanInnov(trajectories(i).observations,arParamK,maParamK);
        if errFlag
            disp('--armaCoefKalman: "armaKalmanInnov" did not function properly!');
            return
        end

        %calculate variance of white noise in current trajectory
        wnVarianceSamp(i) = nanmean(innovation.^2./innovationVar);

        %1st sum in Eq. 3.15
        sum1 = sum1 + nansum(log(innovationVar));
        %2nd sum in Eq. 3.15
        sum2 = sum2 + nansum(innovation.^2./innovationVar);

    end %(for i = 1:length(trajectories))
    
    %calculate -2ln(likelihood)
    neg2LnLikelihoodV = sum1 + totAvail*log(sum2);

    %calculate mean white noise variance of all trajectories
    wnVariance = sum(numAvail.*wnVarianceSamp)/totAvail;

    %get number of parameters estimated: arOrder AR coefficients, maOrder MA
    %coefficients and white noise variance
    numParam = arOrder + maOrder + 1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Least squares fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reformulate the problem as a least squares fitting and obtain the
%variance-covariance matrix of the estimated ARMA coefficients
[varCovMat,arParamL,maParamL,errFlag] = armaLeastSquares(trajOriginal,...
    wnVector,length(arParamK(1,:)),length(maParamK(1,:)),constParam,...
    wnVariance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Post-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if parameters found through least squares fitting are statistically
%equivalent to those found through maximum likelihood estimation with
%Kalman filtering. If they are not equivalent, then results cannot be
%trusted and model is skipped

%prepare input
armaCoef1.arParam = arParamK(1,:);
armaCoef1.maParam = maParamK(1,:);
armaCoef2.arParam = arParamL;
armaCoef2.maParam = maParamL;

%compare parameters
[H,pVCompKL,errFlag] = armaCoefComp(armaCoef1,armaCoef2,varCovMat,varCovMat,...
    'global');
if errFlag
    pVCompKL = 0;
end

%report failure of fit and do not consider results if coefficients are significantly different
if H == 1
%     disp('--armaCoefKalman: Discrepency between least squares and maximum likelihood!')
    errFlag = 1;
end

%use the portmanteau test to check whether residuals are white noise. 
[H,pVPort,errFlag] = portmanteau(wnVector,10,0.01);

%report failure of fit and do not consider results if residuals are not white noise
if H == 1
%     disp('--armaCoefKalman: Residuals did not pass portmanteau test!')
    errFlag = 1;
end


%%%%% ~~ the end ~~ %%%%%
