function [traj,errFlag] = simSetarma(vThresholds,delay,arParam,maParam,...
    noiseSigma,trajLength,trajInit)
%SIMSETARMA simulates a Self Exciting Threshold Autoregressive Moving Average trajectory
%
%SYNOPSIS [traj,errFlag] = simSetarma(vThresholds,delay,arParam,maParam,...
%    noiseSigma,trajLength,trajInit)
%
%INPUT  vThresholds : Column vector of values of thresholds, sorted in increasing order.
%       delay       : Time lag of value compared to vThresholds, set to [] if there is 
%                     only 1 regime (i.e. 0 thresholds).
%       arParam     : Matrix of autoregression parameters. Entries for coefficients
%                     beyond AR order should be filled with NaN.
%       maPAram     : Matrix of moving average parameters. Entries for coefficients
%                     beyond MA order should be filled with NaN.
%       noiseSigma  : Column vector of standard deviations of normally 
%                     distributed white noise used in simulation (noise mean = 0).
%       trajLength  : Length of trajectory to be simulated.
%       trajInit    : Trajectory in first max(max(AR order),delay) time points.
%
%OUTPUT traj        : Simulated trajectory.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('simSetarma')
    disp('--simSetarma: Incorrect number of input arguments!');
    errFlag  = 1;
    traj = [];
    return
end

%check input data
if ~isempty(vThresholds)
    [nThresholds,dummy] = size(vThresholds);
    if dummy ~= 1
        disp('--simSetarma: Variable "vThresholds" should be a column vector!');
        errFlag = 1;
    else
        if min(vThresholds(2:end)-vThresholds(1:end-1)) <= 0
            disp('--simSetarma: Entries in "vThresholds" should be sorted in increasing order, with no two elements equal!');
            errFlag = 1;
        end
    end
    if delay <= 0
        disp('--simSetarma: "delay" should be a positive integer!');
        errFlag = 1;
    end
else
    nThresholds = 0;
    delay = 1;
end
if ~isempty(arParam)
    dummy = size(arParam,1);
    if dummy ~= nThresholds + 1
        disp('--simSetarma: Wrong number of rows in "arParam"!');
        errFlag = 1;
    else
        for i = 1:nThresholds+1
            arOrder(i) = length(find(~isnan(arParam(i,:))));
            r = abs(roots([-arParam(i,arOrder(i):-1:1) 1]));
            if ~isempty(find(r<=1.00001))
                disp('--simSetarma: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
                errFlag = 1;
            end
        end
    end
else
    arOrder = zeros(1,nThresholds+1);
    arParam = zeros(nThresholds+1,1);
end
if ~isempty(maParam)
    dummy = size(maParam,1); 
    if dummy ~= nThresholds + 1
        disp('--simSetarma: Wrong number of rows in "maParam"!');
        errFlag = 1;
    else
        for i = 1:nThresholds+1
            maOrder(i) = length(find(~isnan(maParam(i,:))));
            r = abs(roots([maParam(i,maOrder(i):-1:1) 1]));
            if ~isempty(find(r<=1.00001))
                disp('--simSetarma: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
                errFlag = 1;
            end
        end
    end
else
    maOrder = zeros(1,nThresholds+1);
    maParam = zeros(nThresholds+1,1);
end
[dummy,dummy2] = size(noiseSigma);
if dummy2 ~= 1
    disp('--simSetarma: "noiseSigma" should be a column vector!');
    errFlag = 1;
end
if dummy ~= nThresholds + 1
    disp('--simSetarma: Wrong number of rows in "noiseSigma"!');
    errFlag = 1;
end
if min(noiseSigma) < 0
    disp('--simSetarma: "noiseSigma" should be nonnegative!');
    errFlag = 1;
end
if trajLength <= 0
    disp('--simSetarma: "trajLength" should be a nonnegative integer!');
    errFlag = 1;
end
if length(trajInit) ~= max([arOrder delay])
    disp('--simSetarma: Number of time points initialized should equal max(arOrder,delay)!');
    errFlag = 1;
end
if errFlag
    disp('--simSetarma: Please fix input data!');
    traj = [];
    return
end

%put +/- infinity at ends of thresholds vector
vThresholds = [-Inf; vThresholds; Inf];

%number of previous time points needed due to dependence on the past.
shift = max([arOrder maOrder delay]);

%actual length of trajectory simulated. The added 10*shift time
%points in the beginning are deleted in the end of the calculation and are 
%used only to remove any artificiality in the initial conditions.
tempL = trajLength+10*shift;

%reserve memory for trajectory and noise
traj = zeros(tempL,1);
noise = zeros(tempL,1);

%enter values of first few time steps
traj(shift:-1:shift-max([arOrder delay])+1) = trajInit(end:-1:1);

%obtain trajectory
for i = shift+1:tempL
    level = find(((traj(i-delay)>vThresholds(1:end-1)) + ... %determine level
        (traj(i-delay)<=vThresholds(2:end))) == 2);
    noise(i) = noiseSigma(level)*randn(1);  
    traj(i) = arParam(level,1:arOrder(level))*traj(i-1:-1:i-arOrder(level))... %AR
        + maParam(level,1:maOrder(level))*noise(i-1:-1:i-maOrder(level))...    %MA
        + noise(i);                                                            %noise
end

%get rid of initial 10*shift time points
traj = traj(10*shift+1:end);

%shift trajectory so that it has zero mean
%traj = traj - mean(traj);
