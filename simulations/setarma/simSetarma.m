    function [traj,errFlag] = simSetarma(nThresholds,vThresholds,delay,arOrder,...
    maOrder,arParam,maParam,noiseSigma,trajLength,trajInit)
%SIMSETARMA simulates a Self Exciting Threshold Autoregressive Moving Average trajectory
%
%SYNOPSIS [traj,errFlag] = simSetarma(nThresholds,vThresholds,delay,arOrder,...
%    maOrder,arParam,maParam,noiseSigma,trajLength,trajInit)
%
%INPUT  nThresholds : Number of thresholds.
%       vThresholds : Column vector (of size nThresholds) of values of thresholds, 
%                     sorted in increasing order.
%       delay       : Time lag of value compared to vThresholds.
%       arOrder     : Order of autoregressive part.
%       maOrder     : Order of moving average part.
%       arParam     : (nThresholds+1) by arOrder matrix of autoregression parameters.
%       maPAram     : (nThresholds+1) by maOrder matrix of moving average parameters.
%       noiseSigma  : Standard deviation of noise used in simulation (noise mean = 0).
%       trajLength  : Length of trajectory to be simulated.
%       trajInit    : Trajectory in first max(arOrder,delay) time points.
%
%OUTPUT traj        : Simulated trajectory.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, February 2004

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= nargin('simSetarma')
    disp('--simSetarma: Incorrect number of input arguments!');
    errFlag  = 1;
    traj = [];
    return
end

%check input data
if nThresholds < 0
    disp('--simSetarma: "nThresholds" should be a nonnegative integer!');
    errFlag = 1;
end
if length(vThresholds) ~= nThresholds
    disp('--simSetarma: Number of elements in "vThresholds" should equal "nThresholds"!');
    errFlag = 1;
else
    if min(vThresholds(2:end)-vThresholds(1:end-1)) <= 0
        disp('--simSetarma: Entries in "vThresholds" should be sorted in increasing order, with no two elements alike!');
        errFlag = 1;
    end
end
if nThresholds ~= 0
    if delay <= 0
        disp('--simSetarma: "delay" should be a positive integer!');
        errFlag = 1;
    end
else
    if delay < 0
        disp('--simSetarma: "delay" should be a nonnegative integer!');
        errFlag = 1;
    end
end
    
if arOrder < 0
    disp('--simSetarma: "arOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if maOrder < 0
    disp('--simSetarma: "maOrder" should be a nonnegative integer!');
    errFlag = 1;
end
if arOrder == 0 && maOrder == 0
    disp('--simSetarma: Either "arOrder" or "maOrder" should be > 0!');
    errFlag = 1;
end
if errFlag
    disp('--simSetarma: Please fix input data!');
    traj = [];
    return
end
if arOrder ~= 0
    [nRow,nCol] = size(arParam);
    if nRow ~= nThresholds+1
        disp('--simSetarma: Wrong number of rows in "arParam"!');
        errFlag = 1;
    end
    if nCol ~= arOrder
        disp('--simSetarma: Wrong number of columns in "arParam"!');
        errFlag = 1;
    end
    for i = 1:nThresholds+1
        r = abs(roots([-arParam(i,end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--simSetarma: Causality requires the polynomial defining the autoregressive part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if maOrder ~= 0
    [nRow,nCol] = size(maParam);
    if nRow ~= nThresholds+1
        disp('--simSetarma: Wrong number of rows in "maParam"!');
        errFlag = 1;
    end
    if nCol ~= maOrder
        disp('--simSetarma: Wrong number of columns in "maParam"!');
        errFlag = 1;
    end
    for i = 1:nThresholds+1
        r = abs(roots([maParam(i,end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--simSetarma: Invertibility requires the polynomial defining the moving average part of the model not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
if noiseSigma < 0
    disp('--simSetarma: "noiseSigma" should be nonnegative!');
    errFlag = 1;
end
if trajLength <= 0
    disp('--simSetarma: "trajLength" should be a nonnegative integer!');
    errFlag = 1;
end
if length(trajInit) ~= max(arOrder,delay)
    disp('--simSetarma: Number of time points initialized should equal max(arOrder,delay)!');
    errFlag = 1;
end
if errFlag
    disp('--simSetarma: Please fix input data!');
    traj = [];
    return
end

%put +/- infinity on ends of thresholds for comparison later on
vThresholds = [-Inf; vThresholds; Inf];

%number of previous time points needed due to dependence on the past.
shift = max(max(arOrder,maOrder),delay);

%actual length of trajectory simulated. The added 10*shift time
%points in the beginning are deleted in the end of the calculation and are 
%used only to remove any artificiality in the initial conditions.
tempL = trajLength+10*shift;

%noise in simulation
noise = idinput(tempL,'rgs',[],[-noiseSigma noiseSigma]);

%reserve memory for trajectory
traj = zeros(tempL,1);

%enter values of first few time steps
traj(shift:-1:shift-max(arOrder,delay)+1) = trajInit(end:-1:1);

%obtain trajectory
for t = shift+1:tempL
    level = find(((vThresholds(1:end-1)<=traj(t-delay)) + ... %determine level
        (vThresholds(2:end)>traj(t-delay))) == 2);
    traj(t) = noise(t);                                     %noise at time point t
    for i = 1:arOrder                                       %AR part
        traj(t) = traj(t) + arParam(level,i)*traj(t-i);
    end
    for i = 1:maOrder                                       %MA part
        traj(t) = traj(t) + maParam(level,i)*noise(t-i);
    end
end

%get rid of initial 10*shift time points
traj = traj(10*shift+1:end);

%shift trajectory so that it has zero mean
traj = traj - mean(traj);
