%script test3

trajLength = 10000; %trajectory length
missing = 2000;
arOrder = 4;

maxErr = 0.025; %maximum measurement error standard deviation

[traj1,errFlag] = simSetarma(0,[],0,arOrder,0,[-0.5 -0.1 0.7 0.2],[],0.1,trajLength,[0 0 0 0]);
if errFlag
    return
end

for i = 1:1
    
    err = maxErr*(0.1+0.9*rand(trajLength,1)); %measurement standard deviation
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    
%    [arParam(i,:),varCovMat(:,:,i),residuals(:,i),noiseSigma(i),errFlag] = arlsestim0(traj,arOrder);
    
    traj2 = traj; %copy traj1 into traj2
    
    indx = arOrder+1+floor((trajLength-arOrder-1)*rand(missing,1)); %randomly choose observations to be deleted
    
    traj2(indx,:) = NaN; %delete chosen points
    
    actMiss = length(find(isnan(traj2(:,1)))); %find number of points actually deleted 
    %(this check is needed because sometimes the same number is randomly chosen twice)
    while actMiss < missing %repeat loop until required number of points is deleted
        indx = arOrder+1+floor((trajLength-arOrder-1)*rand(missing-actMiss,1));
        traj2(indx,:) = NaN;
        actMiss = length(find(isnan(traj2(:,1))));
    end
    
    %estimate parameters
    [arParam00(i,:),varCovMat00(:,:,i),residuals00(:,i),noiseSigma00(i),errFlag] = arlsestim0(traj2,arOrder);
    if errFlag
        return
    end
    
    [arParam0(i,:),varCovMat0(:,:,i),residuals0(:,i),noiseSigma0(i),trajP0(:,:,i),errFlag] = arlsIterEstim(traj2,arOrder,0.001,0);
    if errFlag
        return
    end
    [arParam1(i,:),varCovMat1(:,:,i),residuals1(:,i),noiseSigma1(i),trajP1(:,:,i),errFlag] = arlsIterEstim(traj2,arOrder,0.001,1);
    if errFlag
        return
    end
    [arParam2(i,:),varCovMat2(:,:,i),residuals2(:,i),noiseSigma2(i),trajP2(:,:,i),errFlag] = arlsIterEstim(traj2,arOrder,0.001,2);
    if errFlag
        return
    end
    [arParam3(i,:),varCovMat3(:,:,i),residuals3(:,i),noiseSigma3(i),trajP3(:,:,i),errFlag] = arlsIterEstim(traj2,arOrder,0.001,3);
    if errFlag
        return
    end
    [arParam4(i,:),varCovMat4(:,:,i),residuals4(:,i),noiseSigma4(i),trajP4(:,:,i),errFlag] = arlsIterEstim(traj2,arOrder,0.001,4);
    if errFlag
        return
    end
    
    [trajP00(:,:,i),errFlag] = missPointARPred(traj2,arParam00(i,:),0);
    if errFlag
        return
    end
    [trajP01(:,:,i),errFlag] = missPointARPred(traj2,arParam00(i,:),1);
    if errFlag
        return
    end
    [trajP02(:,:,i),errFlag] = missPointARPred(traj2,arParam00(i,:),2);
    if errFlag
        return
    end
    [trajP03(:,:,i),errFlag] = missPointARPred(traj2,arParam00(i,:),3);
    if errFlag
        return
    end
    [trajP04(:,:,i),errFlag] = missPointARPred(traj2,arParam00(i,:),4);
    if errFlag
        return
    end
    
end

%     diff0 = trajP0(:,1)-traj1;
%     diff1 = trajP1(:,1)-traj1;
%     diff2 = trajP2(:,1)-traj1;
%     diff3 = trajP3(:,1)-traj1;
%     diff4 = trajP4(:,1)-traj1;
% %     diff5 = trajP5(:,1)-traj1;
%     diff00 = trajP00(:,1)-traj1;
%     diff01 = trajP01(:,1)-traj1;
%     diff02 = trajP02(:,1)-traj1;
%     diff03 = trajP03(:,1)-traj1;
%     diff04 = trajP04(:,1)-traj1;
% %     diff05 = trajP05(:,1)-traj1;
%     
%     [gamma,lags] = xcov(traj1,500);
%     [gamma0,lags] = xcov(trajP0(:,1),500);
%     [gamma1,lags] = xcov(trajP1(:,1),500);
%     [gamma2,lags] = xcov(trajP2(:,1),500);
%     [gamma3,lags] = xcov(trajP3(:,1),500);
%     [gamma4,lags] = xcov(trajP4(:,1),500);
% %     [gamma5,lags] = xcov(trajP5(:,1),500);
%     [gamma00,lags] = xcov(trajP00(:,1),500);
%     [gamma01,lags] = xcov(trajP01(:,1),500);
%     [gamma02,lags] = xcov(trajP02(:,1),500);
%     [gamma03,lags] = xcov(trajP03(:,1),500);
%     [gamma04,lags] = xcov(trajP04(:,1),500);
% %     [gamma05,lags] = xcov(trajP05(:,1),500);
%     
%     gamma = gamma(201:end);
%     gamma0 = gamma0(201:end);
%     gamma1 = gamma1(201:end);
%     gamma2 = gamma2(201:end);
%     gamma3 = gamma3(201:end);
%     gamma4 = gamma4(201:end);
%     gamma5 = gamma5(201:end);
%     gamma00 = gamma00(201:end);
%     gamma01 = gamma01(201:end);
%     gamma02 = gamma02(201:end);
%     gamma03 = gamma03(201:end);
%     gamma04 = gamma04(201:end);
%     gamma05 = gamma05(201:end);
%     
%     lags = lags(201:end);
%     
%     
%     diffg0 = (abs(gamma0)-abs(gamma))./abs(gamma);
%     diffg1 = (abs(gamma1)-abs(gamma))./abs(gamma);
%     diffg2 = (abs(gamma2)-abs(gamma))./abs(gamma);
%     diffg3 = (abs(gamma3)-abs(gamma))./abs(gamma);
% 
%     [pacfV,errFlag] = pacf(gamma,10);
%     [pacfV0,errFlag] = pacf(gamma0,10);
%     [pacfV1,errFlag] = pacf(gamma1,10);
%     [pacfV2,errFlag] = pacf(gamma2,10);
%     [pacfV3,errFlag] = pacf(gamma3,10);
    