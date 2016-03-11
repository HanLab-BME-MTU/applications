diffCoefRange = [0.025 0.05 0.075 0.1];
receptorDensityRange = [1 2 4 10];
numRepeat = 200;
[time2aggregFree,time2aggregConf3,fracAggregFree,fracAggregConf3] = deal(NaN(numRepeat,length(diffCoefRange),length(receptorDensityRange)));

for k = 1 : length(receptorDensityRange)
    for j = 1 : length(diffCoefRange)
        
        disp([num2str(k) ' ' num2str(j)])
        
        %         modelParamFree = struct('diffCoef',diffCoefRange(j),'confDim',NaN,'receptorDensity',receptorDensityRange(k),'aggregationDist',0.005);
        %         modelParamConf = struct('diffCoef',diffCoefRange(j),'confDim',0.35,'receptorDensity',receptorDensityRange(k),'aggregationDist',0.005);
        modelParamConf3 = struct('diffCoef',diffCoefRange(j),'confDim',0.1,'receptorDensity',receptorDensityRange(k),'aggregationDist',0.005);
        
        for i = 1 : numRepeat
            
            rndSeed = 100 + i;
            simParam = struct('probDim',2,'areaSideLen',10,'timeStep',2e-5,'maxTime',10,'numRep',1,'randNumGenSeeds',[rndSeed rndSeed]);
            
%             [time2aggregFree(i,j,k),fracAggregFree(i,j,k)] = receptorAggregationSuperSimple(modelParamFree,simParam);
            [time2aggregConf3(i,j,k),fracAggregConf3(i,j,k)] = receptorAggregationSuperSimple(modelParamConf3,simParam);
            
        end
        
    end
end
