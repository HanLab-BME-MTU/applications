function [aicV,errFlag] = tarAic(traj,vThresholds,delay,tarParam)
%TARAIC calculates Akaike's Information Criterion to evaluate a model's fit to a given trajectory (which could have missin gdata points).
%
%SYNOPSIS [aicV,errFlag] = tarAic(traj,vThresholds,delay,arParam)
%
%INPUT  traj         : Trajectory to be compared to model.
%                      Missing points should be indicated with NaN.
%       vThresholds  : Column vector of thresholds, sorted in increasing order.
%       delayTest    : Row vector of values of delay parameter.
%       tarParam     : AR parameters in each regime.
%
%OUTPUT aicV         : Value of Akaike's Information Criterion.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, April 2004

%initialize output
errFlag = 0;
aicV = [];

%check if correct number of arguments was used when function was called
if nargin ~= nargin('tarAic')
    disp('--tarAic: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check input data
if ~isempty(vThresholds)
    [nThresholds,dummy] = size(vThresholds);
    if dummy ~= 1
        disp('--tarAic: Variable "vThresholds" should be a column vector!');
        errFlag = 1;
    else
        if min(vThresholds(2:end)-vThresholds(1:end-1)) <= 0
            disp('--tarAic: Entries in "vThresholds" should be sorted in increasing order, with no two elements equal!');
            errFlag = 1;
        end
    end
    if delay <= 0
        disp('--tarAic: "delay" should be a positive integer!');
        errFlag = 1;
    end
else
    nThresholds = 0;
    delay = 1;
end
dummy = size(tarParam,1);
if dummy ~= nThresholds + 1
    disp('--tarAic: Wrong number of rows in "tarParam"!');
    errFlag = 1;
else
    for i = 1:nThresholds+1
        tarOrder(i) = length(find(~isnan(tarParam(i,:))));
        r = abs(roots([-tarParam(i,tarOrder(i):-1:1) 1]));
        if ~isempty(find(r<=1.00001))
            disp('--tarAic: Causality requires the autoregressive polynomial not to have any zeros for z <= 1!');
            errFlag = 1;
        end
    end
end
[trajLength,dummy] = size(traj);
if dummy ~= 1
    disp('--tarAic: "traj" should be a column vector!');
    errFlag = 1;
end
if trajLength < max(tarOrder)
    disp('--tarAic: Length of trajectory should be larger than model order!');
    errFlag = 1;
end

%exit if there are problems in input data
if errFlag
    disp('--tarAic: Please fix input data!');
    return
end

%put -/+ infinity at ends of thresholds vector
vThresholds = [-Inf; vThresholds; Inf];

%find data points to be used in aic calculation and classify them in the different
%regimes (steps 1-5)

%1. get indices of available points
indx = find(~isnan(traj));

%2. discard first max(tarOrder) trajectory points
indx2 = find(indx>max(tarOrder)); %note that indx2 is an array of indices of elements in indx
indx2 = indx2(find(indx2>max(tarOrder)));

%3. find all points whose "delay"-time-steps predecessors are available
indx2 = indx2(find(indx(indx2)>delay));
indx2 = indx2(find(~isnan(traj(indx(indx2)-delay))));

%4. classify points into the different regimes
indxClass = zeros(trajLength,nThresholds+1); %initialize vector
for level = 1:nThresholds+1 %go over all levels (regimes)
    temp = indx2(find((traj(indx(indx2)-delay)>vThresholds(level)) + ...
        (traj(indx(indx2)-delay)<=vThresholds(level+1)) == 2)); %find all points in this level
    aicLength(level) = length(temp);
    indxClass(1:aicLength(level),level) = temp;
end
%remove empty entries in indxClass
indxClass = indxClass(1:max(aicLength),:);
%delete variable indx2
clear indx2;

%5. find all points whose tarOrder(level) previous points are available
aicSet = zeros(size(indxClass));
for level = 1:nThresholds+1
    temp = indx(indxClass(1:aicLength(level),level))-indx(indxClass(1:aicLength(level),level)-tarOrder(level));
    temp = indx(indxClass(find(temp==tarOrder(level)),level));
    aicLength(level) = length(temp);
    aicSet(1:aicLength(level),level) = temp;
end
%remove empty entries in aicSet
aicSet = aicSet(1:max(aicLength),:);
%delete variable temp
clear temp indxClass;

%add to the beginning of each column the number of data points in that regime
aicSet = [aicLength; aicSet];

%initialize residuals vector
residuals = NaN*ones(trajLength,1);

for level = 1:nThresholds+1
    
    %construct matrix of previous points multiplying AR coefficients
    %[size: aicLength by tarOrder(level)]
    prevPoints = zeros(aicLength(level),tarOrder(level));
    for i=1:tarOrder(level)
        prevPoints(:,i) = traj(aicSet(2:aicLength(level)+1,level)-i);
    end
    
    %construct vector of points in this regime
    %[size: aicLength by 1]
    currentPoint = traj(aicSet(2:aicLength(level)+1,level));
    
    %get vector of residuals
    residuals(aicSet(1:aicLength(level),level)) = prevPoints*tarParam(level,...
        1:tarOrder(level))' - currentPoint;
    
    %get standard deviation of residuals
    noiseVar(level) = var(residuals(aicSet(1:aicLength(level),level)));
    
end %(for levels=1:nThresholds+1)

%calculate Akaike's Information Criterion
aicV = aicSet(1,:)*log(noiseVar)' + sum(2*tarOrder+2);
