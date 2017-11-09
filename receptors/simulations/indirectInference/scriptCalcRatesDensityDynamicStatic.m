%Script to calculate interactions rates and cluster densities using
%aggregation information. Results will be saved in same directory as
%compound tracks and aggregation state.
%
%Khuloud Jaqaman, May 2015

sourceRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/varyDissRate/target/targetIS_sT10_dT0p1';
% saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170403/superRes/target';


%Define strings for directory hierarchy as needed
rDDir ={'rD20','rD40','rD60','rD80','rD100','rD120','rD140','rD160'};%,'rD20','rD40','rD60','rD80','rD100','rD120','rD140''rD8',
aPDir = {'aP0p5'};%,'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}
outDirNum =1:30;
dRdir={'dR1'};
lRDir = {'lR1'};%{'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'};
%define space and time information
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: pay attention for the 'areaSideLen','timeStep',and 'sampleStep'
%values. For the 
%LOW DENSITY simulation 'areaSideLen=25' and 
%for HIGH DENSITY it is 12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% modification: Now infoSpaceTime has two new variables: systemState and
% 'intensityInfo'. systemState represents dynamic or static states. 
% The values should be 1 for dynamic and 0 for static. Luciana de Oliveira, February 2017.

infoSpaceTime = struct('probDim',2,'areaSideLen',12,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10],'systemState',0);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
        
        
        %iterate through the different runs
        for outDirIndx = 1 : length(outDirNum)
             for aPDirIndx = 1 : length(aPDir)
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),filesep,lRDir{lRDirIndx}];
          
            systemState=infoSpaceTime.systemState ;
            
            if systemState==1
                %load compTracksAggregState
                tmp = load([currDir,filesep,'compTracksAggregState.mat']);
                compTracksAggregState = tmp.compTracksAggregState;
            elseif systemState==0
                tmp = load([currDir,filesep,'detectionAggregState.mat']);
                compTracksAggregState = tmp.detectionAggregState;
            end
                clear tmp
                
                %get rates and densities
                [rateOnPerClust,rateOffPerClust,densityPerClust,...
    numClustForRateCalc,clustHistory,clustStats] = ...
    clusterOnOffRatesAndDensityDynamicStatic(compTracksAggregState,infoSpaceTime);
                
                
                
                %save results
                
                saveDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
                            
                save([saveDir,'/ratesAndDensity_dt0p1_T10'],'rateOnPerClust',...
                    'rateOffPerClust','densityPerClust','numClustForRateCalc',...
                    'clustHistory','clustStats','-v7.3');
                
                clear compTracksAggregState
                
            end %for each labelRatio
            
        end %for each outDir
        
    end %for each aP
    
 
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx});
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
