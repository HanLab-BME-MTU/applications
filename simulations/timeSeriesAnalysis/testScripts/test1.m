%script test1

trajLength = 1000; %trajectory length
maxErr = 0.0002; %maximum measurement error standard deviation

[traj1,errFlag] = simSetarma(0,[],0,4,0,[0.9 -0.8 0.7 -0.6],[],0.001,trajLength,[1 2 0 -1]);

%no noise, equal weights
[arParam0,noiseSigma0,varCovMat0,errFlag] = arlsestim0(traj1,4);

%noise, different weights
for i=1:100
    err = maxErr*(0.9+0.1*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam01(i,:) = arParam;
    noiseSigma01(i) = noiseSigma;
    varCovMat01(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam11(i,:) = arParam;
    noiseSigma11(i) = noiseSigma;
    varCovMat11(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.8+0.2*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam02(i,:) = arParam;
    noiseSigma02(i) = noiseSigma;
    varCovMat02(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam12(i,:) = arParam;
    noiseSigma12(i) = noiseSigma;
    varCovMat12(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.7+0.3*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam03(i,:) = arParam;
    noiseSigma03(i) = noiseSigma;
    varCovMat03(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam13(i,:) = arParam;
    noiseSigma13(i) = noiseSigma;
    varCovMat13(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.6+0.4*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam04(i,:) = arParam;
    noiseSigma04(i) = noiseSigma;
    varCovMat04(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam14(i,:) = arParam;
    noiseSigma14(i) = noiseSigma;
    varCovMat14(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.5+0.5*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam05(i,:) = arParam;
    noiseSigma05(i) = noiseSigma;
    varCovMat05(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam15(i,:) = arParam;
    noiseSigma15(i) = noiseSigma;
    varCovMat15(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.4+0.6*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam06(i,:) = arParam;
    noiseSigma06(i) = noiseSigma;
    varCovMat06(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam16(i,:) = arParam;
    noiseSigma16(i) = noiseSigma;
    varCovMat16(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.3+0.7*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam07(i,:) = arParam;
    noiseSigma07(i) = noiseSigma;
    varCovMat07(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam17(i,:) = arParam;
    noiseSigma17(i) = noiseSigma;
    varCovMat17(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.2+0.8*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam08(i,:) = arParam;
    noiseSigma08(i) = noiseSigma;
    varCovMat08(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam18(i,:) = arParam;
    noiseSigma18(i) = noiseSigma;
    varCovMat18(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.1+0.9*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam09(i,:) = arParam;
    noiseSigma09(i) = noiseSigma;
    varCovMat09(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam19(i,:) = arParam;
    noiseSigma19(i) = noiseSigma;
    varCovMat19(:,:,i) = varCovMat;
end

for i=1:100
    err = maxErr*(0.001+0.999*rand(trajLength,1)); %measurement standard deviation
    traj = [traj1 err]; %trajectory without perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam10(i,:) = arParam;
    noiseSigma10(i) = noiseSigma;
    varCovMat10(:,:,i) = varCovMat;
    
    perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
    traj = [traj1+perturb err]; %trajectory with perturbation
    [arParam,noiseSigma,varCovMat,errFlag] = arlsestim0(traj,4);
    arParam20(i,:) = arParam;
    noiseSigma20(i) = noiseSigma;
    varCovMat20(:,:,i) = varCovMat;
end

