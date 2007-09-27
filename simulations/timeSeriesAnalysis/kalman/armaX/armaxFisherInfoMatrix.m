function [fishInfoMat,errFlag] = armaxFisherInfoMatrix(trajOut,trajIn,...
    arParam,maParam,xParam,wnVariance)
%ARMAXFISHERINFOMATRIX calculates the Fisher information of an ARMAX model estimated by likelihood maximization
%
%SYNOPSIS [fishInfoMat,errFlag] = armaxFisherInfoMatrix(trajOut,trajIn,...
%    arParam,maParam,xParam,wnVariance)
%
%INPUT  trajOut   : Observations of output time series to be fitted. An 
%                   array of structures with fields:
%           .observations: 2D array of measurements and their uncertainties.
%                          Missing points should be indicated with NaN.
%           .weight      : The weight with which each movie belongs to the
%                          group. Optional. If field doesn't exist,
%                          all weights will be taken as 1.
%                   If there is only one series, it can also be input 
%                   directly as a 2D array. In this case, its weight is 1.
%       trajIn    : Observations of input time series. Either an
%                   array of structures trajIn(1:nTraj).observations, or a
%                   2D array representing one single trajectory.
%           .observations: 2D array of measurements and their uncertainties.
%                   Must not have any missing points.
%                   Enter as [] if there is no input series.
%       arParam   : Autoregressive coefficients (row vector). Enter as []
%                   if there is no AR part.
%       maParam   : Moving average coefficients (row vector). Enter as []
%                   if there is no MA part.
%       xParam    : Coefficients of dependence on input (row vector). Enter
%                   as [] if there is no X part.
%       wnVariance: White noise variance.
%
%OUTPUT fishInfoMat: Fisher information matrix.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%REMARKS Let J be the Fisher Information Matrix, x(t) and V(t) the 
%        innovation and its variance at time t and a(i) the ith parameter 
%        in the model. Then the i,jth element of J is given by
%        J(i,j) = sum(over t=1...N) {(dV(t)/da(i))*(dV(t)/da(j))/(2V(t)^2)
%                                    +(dx(t)/da(i))*(dx(t)/da(j))/V(t)}
%
%MATLAB VERSION (originally written on): 7.0.4.352 (R14) Service Pack 2 
%
%Khuloud Jaqaman, February 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fishInfoMat = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether all input variables were supplied
if nargin < 6
    disp('--armaxFisherInfoMatrix: The function requires 6 input arguments!');
    errFlag  = 1;
    return
end

%check "trajOut" and turn it into struct if necessary
if ~isstruct(trajOut)
    tmp = trajOut;
    clear trajOut
    trajOut.observations = tmp;
    numTraj = 1; %number of trajectories supplied
    trajOut.weight = 1; %weight indicating association with group
    clear tmp
else
    if ~isfield(trajOut,'observations')
        disp('--armaxFisherInfoMatrix: Please input trajOut in fields ''observations''!')
        errFlag = 1;
    end
    numTraj = length(trajOut); %number of trajectories supplied
    if ~isfield(trajOut,'weight') %assign default weights if not supplied
        for i=1:numTraj
            trajOut(i).weight = 1;
        end
    end
end

%add column for observational error of output if not provided
for i=1:numTraj
    traj = trajOut(i).observations;
    [trajLength,nCol] = size(traj);
    if nCol ~= 2
        if nCol == 1 %if no error is supplied, it is assumed that there is no observational error
            traj = [traj zeros(trajLength,1)];
        else
            disp('--armaxFisherInfoMatrix: "trajOut.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
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
        disp('--armaxFisherInfoMatrix: Please input trajIn in fields ''observations''!')
        errFlag = 1;
    end
    
    if ~isempty(trajIn(1).observations)
        for i=1:numTraj %check for observational error
            traj = trajIn(i).observations;
            [trajLength,nCol] = size(traj);
            if nCol ~= 2
                if nCol == 1 %if no error is supplied, assign it the value 0
                    traj = [traj zeros(trajLength,1)];
                else
                    disp('--armaxFisherInfoMatrix: "trajIn.observations" should have either 1 column for measurements, or 2 columns: 1 for measurements and 1 for measurement uncertainties!');
                    errFlag = 1;
                end
            end
            trajIn(i).observations = traj;
        end
    end

end

%exit if there are problems in input data
if errFlag
    disp('--armaxFisherInfoMatrix: Please fix input data!');
    return
end

%obtain number of available observations, per trajectory and total
numAvail = zeros(1,numTraj);
for i=1:numTraj
    traj = trajOut(i).observations(:,1);
    numAvail(i) = length(find(~isnan(traj)));
end
totAvail = sum(numAvail);

%calculate the normalization constant "weight0" for the weights, so that
%sum_i=1^numTraj(weight(i)*numAvail(i)/weight0) = sum_i=1^numTraj(numAvail(i)) = totAvail
weight0 = sum([trajOut.weight].*numAvail)/totAvail;

%normalize the weight with weight0
for i=1:numTraj
    trajOut(i).weight = trajOut(i).weight/weight0;
end

%get orders and number of parameters
arOrder = length(arParam);
maOrder = length(maParam);
xOrder  = length(xParam) - 1;
numParam = arOrder + maOrder + xOrder + 1;

%shift trajectories so that each trajectory's mean = 0
%shift only if pure ARMA (no X)
if xOrder == -1
    for i=1:numTraj
        traj = trajOut(i).observations(:,1);
        traj = traj - nanmean(traj);
        trajOut(i).observations(:,1) = traj;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computation of Fisher information matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get maxOrder = size of matrices and vectors
maxOrder = max(arOrder,maOrder+1);

%add zeros to ends of arParam and maParam to get vectors of length maxOrder
arParamMod = [arParam zeros(1,maxOrder-arOrder)];
maParamMod = [maParam zeros(1,maxOrder-maOrder)];

%construct the transition matrix, F
F = diag(ones(maxOrder-1,1),1);
F(end,:) = arParamMod(end:-1:1); 

%initialize matrix where transition matrix derivatives are stored
FDeriv = zeros(numParam*maxOrder,maxOrder);

%get the transition matrix partial derivatives with respect to the parameters
%go only over the AR coefficients, since the derivative of F wrt MA and X
%coefficients is zero
for i=1:arOrder
    FDeriv(i*maxOrder,maxOrder-i+1) = 1;
end

%get the vector G, which is the covariance of the process at time t 
%with the white noise at time t+i-1, normalized by the white noise variance
G = ones(maxOrder,1);
for i = 2:maxOrder
    dummy = maParamMod(i-1) + arParamMod(1:i-1)*G(i-1:-1:1);
    G(i) = dummy;
end

%initialize the vector where the partial derivatives of G are stored
GDeriv = zeros(numParam*maxOrder,1);

%get the partial derivatives of G wrt the AR coefficients
for i=1:arOrder

    %initialize the temporary vector of derivatives
    derivTemp = zeros(maxOrder,1);
    
    %differentiate AR coefficients wrt arParam(i)
    paramDeriv = zeros(1,maxOrder);
    paramDeriv(i) = 1;

    %differentiate all components of G
    for j=2:maxOrder
        dummy = paramDeriv(1:j-1)*G(j-1:-1:1) ...
            + arParamMod(1:j-1)*derivTemp(j-1:-1:1);
        derivTemp(j) = dummy;
    end

    %store the derivatives in the big vector GDeriv
    GDeriv((i-1)*maxOrder+1:i*maxOrder) = derivTemp;
    
end

%get the partial derivatives of G wrt the MA coefficients
for i=1:maOrder

    %initialize the temporary vector of derivatives
    derivTemp = zeros(maxOrder,1);

    %differentiate MA coefficients wrt maParam(i)
    paramDeriv = zeros(1,maxOrder);
    paramDeriv(i) = 1;

    %differentiate all components of G
    for j=2:maxOrder
        dummy = paramDeriv(j-1) + arParamMod(1:j-1)*derivTemp(j-1:-1:1);
        derivTemp(j) = dummy;
    end

    %store the derivatives in the big vector GDeriv
    GDeriv((i+arOrder-1)*maxOrder+1:(i+arOrder)*maxOrder) = derivTemp;

end

%construct the matrix of dependence on input
B = zeros(maxOrder,xOrder+1);
B(end,:) = xParam(end:-1:1);

%initialize matrix where derivatives of B are stored
BDeriv = zeros(numParam*maxOrder,xOrder+1);
for i=1:xOrder+1
    BDeriv((i+arOrder+maOrder)*maxOrder,xOrder-i+2) = 1;
end

%construct observation vector H
H = zeros(1,maxOrder);
H(1) = 1;

%put numParam copies of F along the diagonal of one big matrix
FBigMat = [];
for i=1:numParam
    FBigMat = blkdiag(FBigMat,F);
end

%put all derivatives of F along the diagonal of one big matrix
FDerivBigMat = [];
for i=1:numParam
    FDerivBigMat = blkdiag(FDerivBigMat,FDeriv((i-1)*maxOrder+1:i*maxOrder,:));
end

%put numParam copies of wnVariance*G*G' along the diagonal of one big matrix
dummy = wnVariance*G*G';
GGprimeBigMat = [];
for i=1:numParam
    GGprimeBigMat = blkdiag(GGprimeBigMat,dummy);
end

%get the contribution of G and its derivatives to the derivative of the
%covariance matrix
GtimesGDeriv = [];
for i=1:numParam
    derivTemp = GDeriv((i-1)*maxOrder+1:i*maxOrder);
    GtimesGDeriv = blkdiag(GtimesGDeriv,wnVariance*(derivTemp*G'+G*derivTemp'));
end

%put numParam copies of B along the diagonal of one big matrix
BBigMat = [];
for i=1:numParam
    BBigMat = blkdiag(BBigMat,B);
end

%put all derivatives of B along the diagonal of one big matrix
BDerivBigMat = [];
for i=1:numParam
    BDerivBigMat = blkdiag(BDerivBigMat,BDeriv((i-1)*maxOrder+1:i*maxOrder,:));
end

%put numParam copies of H'*H along the diagonal of one big matrix
dummy = H'*H;
HtimesHBigMat = [];
for i=1:numParam
    HtimesHBigMat = blkdiag(HtimesHBigMat,dummy);
end

%get the initial state covariance matrix
[stateCovMat00,errFlag] = covKalmanInit(arParam,maParam,G,arOrder,...
    maOrder,maxOrder);
stateCovMat00 = stateCovMat00*wnVariance;

%get the partial derivatives of the initial state covariance matrix
[stateCovMatDeriv00,errFlag] = covKalmanInitDeriv(F,FDeriv,...
    G,GDeriv,wnVariance,xOrder,stateCovMat00,1e-5);

%put numParam copies of covariance matrix along the diagonal of a big matrix
dummy = [];
for j=1:numParam
    dummy = blkdiag(dummy,stateCovMat00);
end
stateCovMat00 = dummy;

%put derivatives of covariance matrix along diagonal of a big matrix
dummy = [];
for j=1:numParam
    dummy = blkdiag(dummy,stateCovMatDeriv00((j-1)*maxOrder+1:j*maxOrder,:));
end
stateCovMatDeriv00 = dummy;

%initialize the Fisher informaiton matrix
fishInfoMat = zeros(numParam);

%go over all trajectories
for i=1:numTraj
    
    %get trajectories
    trajOutI = trajOut(i).observations;
    trajInI  = trajIn(i).observations;
    
    %get trajectory length
    trajLength = size(trajOutI,1);

    %add maxOrder zeros to end of input series
    trajInI = [trajInI; zeros(maxOrder,2)];

    %assign initial state vector and it derivatives
    stateVecT_T = zeros(maxOrder*numParam,1);
    stateVecDerivT_T = zeros(maxOrder*numParam,1);

    %assign initial state covariance matrix and its derivatives
    stateCovMatT_T= stateCovMat00;
    stateCovMatDerivT_T = stateCovMatDeriv00;

    %initialize vector of contributions from input
    inputContr = zeros(numParam*maxOrder,1);
    inputContrDeriv = zeros(numParam*maxOrder,1);

    %go over all time points in the trajectory
    for j=1:trajLength

        %calculate contribution of input series
        if ~isempty(xParam)
            
            %get modified xOrder to deal with first few points
            modXOrder = j - 1 + maxOrder - max(1,j-1+maxOrder-xOrder);

            %get input vector
            inputVector = trajInI(j-1+maxOrder-modXOrder:j-1+maxOrder,1);
            inputVector = [zeros(xOrder+1-length(inputVector),1); inputVector];
            inputVector = repmat(inputVector,numParam,1);

            %get contribution of input to the state vector and its derivative
            inputContr = BBigMat*inputVector;
            inputContrDeriv = BDerivBigMat*inputVector;
            
        end

        %predict the state vector
        stateVecT1_T = FBigMat*stateVecT_T + inputContr;

        %get the derivative of the predicted state vector
        stateVecDerivT1_T = FDerivBigMat*stateVecT_T + ...
            FBigMat*stateVecDerivT_T + inputContrDeriv;
        
        %obtain the predicted state's covariance matrix
        stateCovMatT1_T = FBigMat*stateCovMatT_T*FBigMat' + GGprimeBigMat;

        %get the partial derivatives of the covariance matrix
        stateCovMatDerivT1_T = FDerivBigMat*stateCovMatT_T*FBigMat' + ...
            FBigMat*stateCovMatDerivT_T*FBigMat' + ...
            FBigMat*stateCovMatT_T*FDerivBigMat' + GtimesGDeriv;

        %update covariance matrix given observation at current time point
        if isnan(trajOutI(j,1)) %if observation at this time point is missing

            %"update" state vector and its derivative
            stateVecT_T = stateVecT1_T;
            stateVecDerivT_T = stateVecDerivT1_T;
            
            %"update" state covariance matrix and its derivatives
            stateCovMatT_T = stateCovMatT1_T;
            stateCovMatDerivT_T = stateCovMatDerivT1_T;

        else %if there is an observation
            
            %predict observable
            observableP = stateVecT1_T(1);

            %calculate innovation
            innovation = trajOutI(j,1) - observableP;

            %get the innovation variance
            innovationVar = stateCovMatT1_T(1,1) + trajOutI(j,2)^2;
            
            %get the partial derivatives of the innovation
            innovDeriv = -stateVecDerivT1_T(1:maxOrder:end);
            
            %get the partial derivatives of the innovation variance
            innovVarDeriv = diag(stateCovMatDerivT1_T);
            innovVarDeriv = innovVarDeriv(1:maxOrder:end);
            
            %calculate delta*H and delta
            deltaH = stateCovMatT1_T*HtimesHBigMat/innovationVar;
            delta = deltaH*repmat(H',numParam,1);
            
            %calculate (derivative of delta) * H and (derivative of delta)
            deltaDerivH = stateCovMatDerivT1_T*HtimesHBigMat/innovationVar ...
                - stateCovMatT1_T*HtimesHBigMat*stateCovMatDerivT1_T*...
                HtimesHBigMat/innovationVar^2;
            deltaDeriv = deltaDerivH*repmat(H',numParam,1);

            %update the state vector
            stateVecT_T = stateVecT1_T + delta*innovation;

            %update the derivative of the state vector
            innovDerivTmp = repmat(innovDeriv',maxOrder,1);
            stateVecDerivT_T = stateVecDerivT1_T + deltaDeriv*innovation ...
                + delta.*innovDerivTmp(:);

            %update state covariance matrix and make sure it's symmetric
            stateCovMatT_T = stateCovMatT1_T - deltaH*stateCovMatT1_T;
            stateCovMatT_T = (stateCovMatT_T+stateCovMatT_T')/2;

            %update derivatives of state covariance matrix
            stateCovMatDerivT_T = stateCovMatDerivT1_T - deltaDerivH*...
                stateCovMatT1_T - deltaH*stateCovMatDerivT1_T;
            stateCovMatDerivT_T = (stateCovMatDerivT_T+stateCovMatDerivT_T')/2;
            
            %update the Fisher information matrix
            fishInfoMat = fishInfoMat + (innovVarDeriv*innovVarDeriv'/...
                innovationVar^2/2 + innovDeriv*innovDeriv'/innovationVar)*...
                trajOut(i).weight;

        end %(if isnan(trajOutI(j,1)) ... else ...)

    end %(for j=1:trajLength)

end %(for i=1:numTraj)


%%%%% ~~ the end ~~ %%%%%




% %
% %Calculate, for all i and j, the
% %sum over t=1...N of { (1/V(t)) * < dyTilde(t)/dtheta(i) * dyTilde(t)/dtheta(j) > } 
% %
% %This problem is solved by contructing a new state vector 
% %Z* = [Z; dZ/dtheta1; dZ/dtheta2; ... dZ/dthetan] and considering its
% %equation of state and observational equation
% %
% 
% %construct new transition matrix
% FStar = blkdiag(F,FBigMat);
% FStar(maxOrder+1:end,1:maxOrder) = FDeriv;
% 
% %construct new process-error covariance vector
% GStar = [G; GDeriv];
% 
% %construct new observation vector
% HStar = zeros(1,maxOrder*(numParam+1));
% HStar(1) = 1;
% 
% %devise matrix for selecting terms from state covariance matrix 
% %needed for constructing Fisher information matrix
% selectionMat = zeros(maxOrder*numParam,numParam);
% for i=1:numParam
%     selectionMat((i-1)*maxOrder+1,i) = 1;
% end
% 
% %get the initial state covariance matrix
% [stateCovMat00,errFlag] = covKalmanInitGen(FStar,GStar,wnVariance,1e-5);
% 
% %go over all trajectories
% for i=1:numTraj
%     
%     %get trajectories
%     trajOutI = trajOut(i).observations;
%     trajInI  = trajIn(i).observations;
%     
%     %get trajectory length
%     trajLength = size(trajOutI,1);
% 
%     %initialize state vector and its covariance matrix
%     ZStarT_T = zeros(maxOrder*(numParam+1),1);
%     PStarT_T = stateCovMat00;
% 
%     %go over all time points in the trajectory
%     for j=1:trajLength
% 
%         %predict state at time j given state at time j-1
%         ZStarT1_T = FStar*ZStarT_T;
% 
%         %obtain the predicted state's covariance matrix
%         PStarT1_T = FStar*PStarT_T*FStar' ...
%             + wnVariance*GStar*GStar';
%         PStarT1_T = (PStarT1_T+PStarT1_T')/2;
% 
%         %predict observable at time j
%         observableP = ZStarT1_T(1);
% 
%         if isnan(trajOutI(j,1)) %if observation at this time point is missing
% 
%             %cannot modify state vector and its covariance matrix predicted
%             %from previous timepoint since there is no observation
%             ZStarT_T = ZStarT1_T;
%             PStarT_T = PStarT1_T;
% 
%         else %if there is an observation
% 
%             %get innovation
%             innovation = trajOutI(j,1) - observableP;
% 
%             %and its variance
%             innovationVar = PStarT1_T(1,1) + trajOutI(j,2)^2;
% 
%             %calculate delta
%             delta = PStarT1_T*HStar'/innovationVar;
% 
%             %modify state vector prediction using observation
%             ZStarT_T = ZStarT1_T + delta*innovation;
% 
%             %update state covariance matrix and make sure it's symmetric
%             PStarT_T = PStarT1_T - delta*HStar*PStarT1_T;
%             PStarT_T = (PStarT_T+PStarT_T')/2;
%             
%             %update the Fisher information matrix
%             fishInfoMat = fishInfoMat + selectionMat'*PStarT1_T(...
%                 maxOrder+1:end,maxOrder+1:end)*selectionMat/innovationVar;
% 
%         end %(if isnan(trajOutI(j,1)) ... else ...)
% 
%     end %(for j=1:trajLength)
% 
% end %(for i=1:numTraj)
