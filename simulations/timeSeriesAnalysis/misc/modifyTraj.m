
function [trajMod,errFlag] = modifyTraj(traj,interval0)

errFlag = 0;

%check if correct number of arguments were used when function was called
if nargin ~= 2
    disp('--modifyTraj: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check trajectory and turn into struct if necessary
if ~isstruct(traj)
    tmp = traj;
    clear traj
    traj.observations = tmp;
    clear tmp
elseif ~isfield(traj,'observations')
    disp('--modifyTraj: Please input the trajectories in fields ''observations''')
    errFlag = 1;
    return
end

%initialize trajectory to be modified
trajMod = traj;

for i=1:length(traj)

    %get indices of missing observations
    indxMiss = [0; find(isnan(traj(i).observations(:,1))); length(traj(i).observations(:,1))+1];

    %evaluate lengths of intervals between missing observations
    interval = diff(indxMiss) - 1;

    %find intervals that are > 0 and <= "interval0"
    indxShort = find((interval<=interval0) & (interval>0));

    %substitute NaN for the existing measurement values in these intervals
    %(<=interval0)
    for j=1:length(indxShort)
        trajMod(i).observations(indxMiss(indxShort(j))+1:indxMiss(indxShort(j)+1)-1,:) = NaN;
    end

    %remove missing observation at the beginning of a trajectory
    while ~isempty(trajMod(i).observations) & isnan(trajMod(i).observations(1,1))
        trajMod(i).observations = trajMod(i).observations(2:end,:);
    end
    
    %remove missing observation at the end of a trajectory
    while ~isempty(trajMod(i).observations) & isnan(trajMod(i).observations(end,1))
        trajMod(i).observations = trajMod(i).observations(1:end-1,:);
    end

   trajMod(i).length = length(trajMod(i).observations(:,1));
   trajMod(i).numAvail = length(find(~isnan(trajMod(i).observations(:,1))));
   trajMod(i).percentMiss = 100*(1-trajMod(i).numAvail/trajMod(i).length);
    
end
