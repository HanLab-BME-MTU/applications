initialState = struct('mtLength0',0.44,'capSize0',3,'unitConc',20);

modelParam = struct('minLength',0.2,'maxLength',1.5,'kTOnTFree',0.4,...
    'addAmpT',0.2,'addWidT',10,'addLenT',0.65,'kTOff',0.0,'kTOnD',0.2,...
    'kDOffFree',14,'addAmpD',7,'addWidD',10,'addLenD',0.85,'kHydrolysis',10);

kTOnDEff = modelParam.kTOnD*initialState.unitConc;
kDOffFree = modelParam.kDOffFree;
kHydrolysis = modelParam.kHydrolysis;

runInfo = struct('maxNumSim',1,'totalTime',10000,'simTimeStep',0.005,...
    'timeEps',0.5,'expTimeStep',1,'aveInterval',0.6,'analyzeOpt',2);

saveTraj = struct('saveOrNot',0,'fileName',[]);
saveStats = struct('saveOrNot',0,'fileName',[]);

expData(1,:) = [1.86 0.03];  %growth speed (micrometers/min)
expData(2,:) = [1.69 0.03];  %shrinkage speed (micrometers/min)
expData(3,:) = [0.27 0.1];   %catastrophe frequency (s^-1)
expData(4,:) = [0.27 0.1];   %rescue frequency (s^-1)
expData(5,:) = [0.48 0.01];  %minumum length (micrometers)
expData(6,:) = [1.12 0.02];  %maximum length (micrometers)

objective = zeros(9,13)+NaN;

i = 0;
for kTOnTFree = [0.25:0.05:0.65]
    i = i + 1;
    j = 0;
    for addAmpT = [0.05:0.05:kTOnTFree-0.05]
        j = j + 1;
        
        try
            
            modelParam.kTOnTFree = kTOnTFree;
            modelParam.addAmpT = addAmpT;
            
            kTOnTFreeEff = modelParam.kTOnTFree*initialState.unitConc;
            
            runInfo.simTimeStep = min(0.99*runInfo.timeEps/max(...
                [kTOnTFreeEff kTOnDEff kDOffFree kHydrolysis]),runInfo.aveInterval/2);   
            
            [errFlag,dataStats] = analyzeMtTrajectory(2,modelParam,initialState,runInfo,...
                saveTraj,saveStats);
            if errFlag
                error('--objectiveGTPCapLDepK: analyzeMtTrajectory did not finish successfully!');
            end
            
            gSpeed = dataStats.growthSpeed;
            if isempty(gSpeed)
                gSpeed(1) = Inf;
            end
            sSpeed = abs(dataStats.shrinkageSpeed);
            if isempty(sSpeed)
                sSpeed(1) = Inf;
            end
            catFreq = dataStats.catFreq;
            if isempty(catFreq)
                catFreq = Inf;
            end
            resFreq = dataStats.resFreq;
            if isempty(resFreq)
                resFreq = Inf;
            end
            minLength = dataStats.minDistanceM5;
            maxLength = dataStats.maxDistanceM5;
            
            objective(i,j) = ((gSpeed(1)-expData(1,1))/expData(1,1))^2 + ((sSpeed(1)-expData(2,1))/expData(2,1))^2 ...
                + ((catFreq(1)-expData(3,1))/expData(3,1))^2 + ((resFreq(1)-expData(4,1))/expData(4,1))^2 ...
                + ((minLength(1)-expData(5,1))/expData(5,1))^2 + ((maxLength(1)-expData(6,1))/expData(6,1))^2;
            
        catch
        end
        
    end        
end