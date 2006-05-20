function [trajSamp,errFlag] = sampleTraj(traj0,sampInt,avgInt)
%SAMPLETRAJ samples a trajectory at a specified rate, with possible averaging
%
%SYNOPSIS [trajSamp,errFlag] = sampleTraj(traj0,sampInt,avgInt)
%
%INPUT  traj0   : 2-column vector of times [s] and lengths [microns].
%       sampInt : Interval between sampling points.
%       avgInt  : Interval over which to average. Optional. Default: 0.
%
%OUTPUT trajSamp: Sampled trajectory. 2- or 3-column vector of times [s],
%                 lengths [microns] and length std [microns] in case of
%                 averaging.
%       errFlag : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, May 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trajSamp = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check input data
if nargin < 2
    disp('--sampleTraj: Wrong number of input variables!');
    errFlag = 1;
    return
end

%assign default if avgInt is not given
if nargin < 3
    avgInt = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%store original sampling interval
sampIntOrig = sampInt;

%shift time to start at 0
timeInit = traj0(1,1);
traj0(:,1) = traj0(:,1) - timeInit;

%in case of averaging, if the time step is not uniform or not small
%enough for averaging, determine an appropriate "pre-sampling" interval
tmp = diff(traj0(:,1));
if avgInt ~= 0 && (std(tmp) ~= 0 || mean(tmp) > avgInt/2)
    %indicate that presampling at a lower rate is needed before
    %sampling with averaging is done
    presample = 1;
    %get new sampling interval to allow averaging
    sampInt = avgInt/5;
else
    presample = 0;
end

%if there no averaging or presampling is needed, sample with specified 
%interval without averaging
if avgInt == 0 || presample == 1

    %get total time of original trajectory
    totalTime = traj0(end,1);

    %calculate number of time points to sample at
    numTP = floor(totalTime/sampInt);

    %obtain time at which to sample
    samplingTime = [0:numTP]'*sampInt;

    %get index in traj0 that corresponds to the transition points just before
    %each sampling time point.
    indxTran = 0;
    for i=1:length(traj0)
        indxTran = indxTran + (samplingTime-traj0(i,1)>=0);
    end

    %reserve memory for trajTemp
    trajTemp = zeros(numTP+1,1);

    %calculate trajectory value at each sampling point
    for i=1:numTP+1
        if indxTran(i) ~= length(traj0)
            tmp = traj0(indxTran(i)+1,1) - traj0(indxTran(i),1);
            trajTemp(i) = (traj0(indxTran(i)+1,1)-samplingTime(i))*traj0(indxTran(i),2)/tmp ...
                + (samplingTime(i)-traj0(indxTran(i),1))*traj0(indxTran(i)+1,2)/tmp;
        else
            trajTemp(i) = traj0(end,2);
        end
    end

    %add time column to trajTemp, and store in traj0
    traj0 = [samplingTime trajTemp];

    %assign back original sampling interval
    sampInt = sampIntOrig;
    
end

%if there is averaging, sample at user-specified rate with averaging
if avgInt ~= 0

    %get current sampling interval
    dt = traj0(2,1) - traj0(1,1);
    
    %remove time column from traj0
    traj0 = traj0(:,2);
    
    %calculate ratio of required sampling time step to input time step
    ratio = round(sampInt/dt);

    %get number of entries in original trajectory
    vectorLength = length(traj0);

    %determine number of entries to be removed from end of traj0, to make
    %its length an integer multiple of ratio
    trash = rem(vectorLength,ratio);

    %remove extra entries from the end
    traj0 = traj0(1:vectorLength-trash,:);

    %get new number of entries in traj0
    vectorLength = length(traj0);

    %get number of entries in new, sampled trajectory
    trajLength = floor(vectorLength/ratio);

    %reshape mtLength to get "ratio" trajectories of "trajLength" entries each
    trajectories = reshape(traj0,ratio,trajLength);

    %obtain ratio of averaging interval to simulation time step, rounded to
    %nearest integer
    ratio = round(avgInt/dt);

    %get interval that is averaged over
    trajectories = trajectories(1:ratio,:);

    %get average length and its std at eachtime point
    trajSamp = [mean(trajectories)' max(1e-15,abs(std(trajectories)'))];

    %add time column
    trajSamp = [[0:size(trajSamp,1)-1]'*sampInt trajSamp];
    
else %if there is no averaging
    
    trajSamp = traj0;

end

%shift starting time back to what was given in input
trajSamp(:,1) = trajSamp(:,1) + timeInit;

%%%%% ~~ the end ~~ %%%%%

