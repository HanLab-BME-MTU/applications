function batchRun1

tic;

cd M:\Kdata2;

fid = fopen(['logFile-' nowString '.txt'],'w');

initialState = struct('mtLength0',0.44,'capSize0',5,'unitConc',20);
unitConc = initialState.unitConc;

maxBigI = 3^8;
paramVec = zeros(maxBigI,10);
simTimeStep = zeros(maxBigI,1);

bigI = 0;
for kTOnTFree = [0.04 0.4 4] %vary kTonTFree
    addAmpT = kTOnTFree/2;
    for addWidT = [10 50 100]; %vary addWidT
        for addLenT = [0.3 0.45 0.6] %vary addLenT
            for kTOnD = [0.02 0.2 2] %vary kTOnD
                for kDOffFree = [1.4 14 140] %vary kDOffFree
                    addAmpD = kDOffFree/2;
                    for addWidD = [10 50 100] %vary addWidD
                        for addLenD = [0.6 0.75 0.9] %vary addLenD
                            for kHydrolysis = [1.2 12 120] %vary kHydrolysis
                                
                                bigI = bigI + 1;
                                
                                paramVec(bigI,:) = [kTOnTFree addAmpT addWidT addLenT kTOnD kDOffFree addAmpD addWidD addLenD kHydrolysis]; 
                                    
                                kTOnTFreeEff = kTOnTFree*unitConc;%effective rate constants of "unit" addition
                                kTOnDEff = kTOnD*unitConc;
                                
                                simTimeStep(bigI) = min(0.99*0.2/max([kTOnTFreeEff kTOnDEff kDOffFree kHydrolysis]),0.1);
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

save('runLogFile', 'paramVec', 'simTimeStep');

for bigI = maxBigI:-1:1
    
    fprintf(fid,'--Run %i:\n', bigI);
    
    kTOnTFree = paramVec(bigI,1);
    addAmpT = paramVec(bigI,2);
    addWidT = paramVec(bigI,3);
    addLenT = paramVec(bigI,4);
    kTOnD = paramVec(bigI,5);
    kDOffFree = paramVec(bigI,6);
    addAmpD = paramVec(bigI,7);
    addWidD = paramVec(bigI,8);
    addLenD = paramVec(bigI,9);
    kHydrolysis = paramVec(bigI,10);
    dt = simTimeStep(bigI);
    
    try
        
        modelParam = struct('minLength',0.2,'maxLength',1','kTOnTFree',kTOnTFree,...
            'addAmpT',addAmpT,'addWidT',addWidT,'addLenT',addLenT,'kTOff',0.0,'kTOnD',...
            kTOnD,'kDOffFree',kDOffFree,'addAmpD',addAmpD,'addWidD',addWidD,'addLenD',...
            addLenD,'kHydrolysis',kHydrolysis);
        
        runInfo = struct('maxNumSim',1,'totalTime',10000,'simTimeStep',dt,...
            'timeEps',0.2,'expTimeStep',1,'aveInterval',0.6);
        
        saveFile = sprintf('mtStat%i',bigI);
        
        errFlag = analyzeMtTrajectory(2,modelParam,initialState,runInfo,saveFile);
        
        if errFlag
            fprintf(fid,'    did not finish successfully.\n\n');
        else
            fprintf(fid,'    finished successfully.\n\n');
        end
        
    catch
        fprintf(fid,'    %s\n',lasterr);
        fprintf(fid,'    did not finish successfully.\n\n');
    end
    
end

fclose(fid);

toc;