function modeFracAboveBelowDivider = diffCoefDistr(diffModeCoef,numTraj,trajLength,diffModeDivider)

%get number of diffusion modes
numMode = length(diffModeCoef);

%simulate trajectories in each mode
traj = NaN(trajLength,2,numTraj,numMode);
for iMode = 1 : numMode
    for iTraj = 1 : numTraj
        trajTmp = brownianMotion(2,diffModeCoef(iMode),trajLength,0.01);
        trajTmp = trajTmp(1:100:end,:);
        trajTmp = trajTmp(2:end,:);
        traj(:,:,iTraj,iMode) = trajTmp;
    end
end
msd = squeeze( mean( sum( diff(traj).^2,2 ) ) );
diffCoef = msd / 4;

%display a histogram of the diffuion coefficients and indicate mode
%dividers
figure
hist(log10(diffCoef),100)
hold on
for iMode = 1 : numMode-1
    plot(log10(diffModeDivider(iMode))*[1 1],[0 250])
end

%calculate fraction of each mode included within dividers
modeFracAboveBelowDivider = NaN(numMode,2);
modeFracAboveBelowDivider(1,2) = length(find(diffCoef(:,1)<diffModeDivider(1)))/size(diffCoef,1);
for i = 2 : numMode-1
    modeFracAboveBelowDivider(i,1) = length(find(diffCoef(:,i)>diffModeDivider(i-1)))/size(diffCoef,1);
    modeFracAboveBelowDivider(i,2) = length(find(diffCoef(:,i)<diffModeDivider(i)))/size(diffCoef,1);
end
i = numMode;
modeFracAboveBelowDivider(i,1) = length(find(diffCoef(:,i)>diffModeDivider(i-1)))/size(diffCoef,1);
