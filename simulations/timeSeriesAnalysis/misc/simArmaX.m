function [trajOut,errFlag] = simArmaX(arParam,maParam,xParam,noiseSigma,...
    trajIn)
%SIMARMAX generates an Autoregressive Moving Average series with depedence on an input series
%
%SYNOPSIS [trajOut,errFlag] = simArmaX(arParam,maParam,xParam,noiseSigma,...
%    trajIn)
%
%INPUT  arParam   : Row vector of autoregressive coefficients.
%       maParam   : Row vector of moving average coefficients.
%       xParam    : Row vector of coefficients of dependence on input.
%       noiseSigma: Standard deviation of white noise. Assumed to be normally
%                   distributed with mean zero and given std.
%       trajIn    : Input series.
%
%OUTPUT trajOut   : Generated trajectory.
%       errFlag   : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, January 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trajOut = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if correct number of arguments were used when function was called
if nargin ~= nargin('simArmaX')
    disp('--simArmaX: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check AR coefficients
if ~isempty(arParam)
    [nRow,arOrder] = size(arParam);
    if nRow ~= 1
        disp('--simArmaX: "arParam" should be a row vector!');
        errFlag = 1;
    else
        r = abs(roots([-arParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--simArmaX: Causality requires the AR polynomial not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
else
    arOrder = 0;
    arParam = zeros(1,0);
end

%check MA coefficients
if ~isempty(maParam)
    [nRow,maOrder] = size(maParam);
    if nRow ~= 1
        disp('--simArmaX: "maParam" should be a row vector!');
        errFlag = 1;
    else
        r = abs(roots([maParam(end:-1:1) 1]));
        if ~isempty(find(r<=1))
            disp('--simArmaX: Invertibility requires the MA polynomial not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
else
    maOrder = 0;
    maParam = zeros(1,0);
end

%check X coefficients (indicating dependence on input)
if ~isempty(xParam)
    [nRow,xOrder] = size(xParam);
    if nRow ~= 1
        disp('--simArmaX: "xParam" should be a row vector!');
        errFlag = 1;
    end
    xOrder = xOrder - 1;
else
    xOrder = -1;
    xParam = zeros(1,0);
end

%make sure that white noise standard deviation is not negative
if noiseSigma < 0
    disp('--simArmaX: "noiseSigma" should be nonnegative!');
    errFlag = 1;
end

%get length of input series = length of output series
trajLength = length(trajIn);

%exit if there are problems in input data
if errFlag
    disp('--simArmaX: Please fix input data!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trajectory generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reserve memory for trajectory and noise
trajOut = zeros(trajLength,1);

%initialize normal random number generator
randn('state',sum(100*clock));

%get white noise vector
noise = noiseSigma*randn(trajLength,1);

%construct trajectory
for i = max([arOrder maOrder xOrder])+1:trajLength
    trajOut(i) = arParam*trajOut(i-1:-1:i-arOrder) ...% AR
        + maParam*noise(i-1:-1:i-maOrder) ...         % MA
        + noise(i) ...                                % noise
        + xParam*trajIn(i:-1:i-xOrder);               % X
end


%%%%% ~~ the end ~~ %%%%%
