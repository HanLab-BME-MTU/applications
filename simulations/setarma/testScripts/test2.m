%script test2

trajLength = 1000; %trajectory length
maxErr = 2.5; %maximum measurement error standard deviation

count0 = 0;
count1 = 0;
count2 = 0;

for j = 1:100
    
    [traj1,errFlag] = simSetarma(0,[],0,4,0,[0.9 -0.8 0.7 -0.6],[],10,trajLength,[1 2 0 -1]);
    
    %no noise, equal weights
    [arParam0(j,:),noiseSigma0(j),varCovMat0(:,:,j),errFlag,counter] = arlsestim0(traj1,4);
    if counter
        count0 = count0 + 1;
    end
   
    %noise, different weights
    for i=1:100
        err = maxErr*(0.05+0.95*rand(trajLength,1)); %measurement standard deviation
        traj = [traj1 err]; %trajectory without perturbation
        [arParam,noiseSigma,varCovMat,errFlag,counter] = arlsestim0(traj,4);
        arParam10(i,:,j) = arParam;
        noiseSigma10(i,j) = noiseSigma;
        varCovMat10(:,:,i,j) = varCovMat;
        if counter
            count1 = count1 + 1;
        end
        
        perturb = (2*rand(trajLength,1)-1)*2.*err; %perturbation due to measurement error
        traj = [traj1+perturb err]; %trajectory with perturbation
        [arParam,noiseSigma,varCovMat,errFlag,counter] = arlsestim0(traj,4);
        arParam20(i,:,j) = arParam;
        noiseSigma20(i,j) = noiseSigma;
        varCovMat20(:,:,i,j) = varCovMat;
        if counter
            count2 = count2 + 1;
        end
    end
    
end

