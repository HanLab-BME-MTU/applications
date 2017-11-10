% This script is for the calculation of Mahalonobis distance and p-value
% considering multiple entries of dynamic and static data.


%% Input
% TargetRoot for the multiple entries, it is a general script and it can be
% both dynamic and static data. Important note that the targetRoot now is
% write as a cell array, it allow the script to decide how to do the
% calculations for the different system states.

% target source root
targetRoot=cell(1,1);
targetRoot{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170526/superRes/bootstrapping/analysis/targetIS_sT25_dT0p1';
targetRoot{2}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170316/superResol/staticData/analysis/targetIS_sT25_dT0p1';
targetRoot{3}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170313/staticData/target/analysis';
targetRoot{4}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/target/analysis';
targetRoot{5}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/target/analysis';


% probe source root
sourceRootProbe=cell(1,1);
sourceRootProbe{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/varyDissRate/analysis/probeIS_sT25_dT0p1';
% sourceRootProbe{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/varyDissRate/analysis/probeIS_sT25_dT0p1';
%  sourceRootProbe{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170313/staticData/probe/analysis';
% sourceRootProbe{4} ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/probe/analysis';
% sourceRootProbe{5} ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/probe/analysis';



saveDirectory='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170526/results/staticData';
%
%
% % name of the target
rDtarget ={'rD4'};%,'rD20','rD40','rD60','rD80','rD100','rD120','rD140''rD8',
aPtarget = {'aP0p5'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
dRDirT={'dR0p5'};
dRDir={'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};%,'dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'
% the labeled fraction is write as a cell array, with each line is
% equivalent for the same set of dynamic and static data and each column is
% for different sets.
bootRepSize=100;
lRtarget = {'lR1'};
% this vector identifies the system state and it is equivalent to each
% different entry in the targetRoot. It also have the length of label
% ratio.

stateSystemVector= 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBE INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%probes

%Define strings for directory hierarchy as needed
rDDir = {'rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%,{'rD20','rD40','rD60','rD80','rD100','rD120','rD140'},'rD4','rD6','rD8','rD10','rD12','rD14','rD16'
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'


% labeled ratio for dynamic and static data
lRProbe = {'lR1'};%{'lR0p91';'lR0p93';'lR0p95';'lR0p97';'lR0p99';'lR1'
%     'lR0p1','lR1'
%     'lR0p2','lR1'
%     'lR0p3','lR1'
%     'lR0p4','lR1'
%     'lR0p5','lR1'
%     'lR0p6','lR1'
%            'lR0p04','lR0p08','lR0p09','lR0p16';

%

%% CALCULATIONS


%In this version I am calculating only one value rD and aP for the target.

sMatrixBoot=cell(bootRepSize,1);
pMatrixBoot=cell(bootRepSize,1);

for iboot=1: bootRepSize
    
    for apDirIndxT=1 : length (aPtarget)
        
        for dRindexT=1: length (dRDirT)
            for rDDirIndx = 1 : length(rDDir)
                
                for dRDirIndx=1: length(dRDir)
                    
                    %Iterate through association probability values per density
                    for aPDirIndx = 1 : length(aPDir)
                        
                        for lRIndex=1:size(lRProbe,1)
                            
                            for dataStateIndex=1: length(stateSystemVector)
                                
                                %define if it is dynamic or static data
                                systemState=stateSystemVector(dataStateIndex);
                                
                                
                                
                                if systemState==1
                                    
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DYNAMIC DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
                                    % target
                                    % load target paramMatrix
                                    
                                    % as each entry of the stateSystemVector coresponds to
                                    % one different value of label ratio, here the lRtarget
                                    % index is given by the dataStateIndex.
                                    
                                    currDir=[targetRoot{1},filesep,rDtarget{1},filesep,...
                                        aPtarget{apDirIndxT},filesep,lRtarget{dataStateIndex}];
                                    tmp =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                    paramMatrixTarget{1}=tmp.paramMatrix;
                                    
                                    %probe
                                    
                                    % load probe paramMatrix
                                    
                                    % Now the lRprobe is an array, and each row corresponds
                                    % to a different entry stateSystemVector, so to have the lRprobe values
                                    % corresponding to each probe we need to index for the lR.
                                    
                                    %                             currDir = [sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
                                    %                               dRDir{1},  aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                    currDir=[sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
                                        dRDir{dRDirIndx},filesep,aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                    
                                    %load ratesAndDensity that has the paramMatrix
                                    tmp = load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                    paramMatrixProbe{dataStateIndex}=tmp.paramMatrix;
                                    
                                    % calculate theta and V for the data set
                                    
                                    [thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVaramCovMatrixForMahalonobisDistanceLastFrame(paramMatrixTarget{dataStateIndex,1},paramMatrixProbe{dataStateIndex,1},systemState(dataStateIndex));
                                    
                                    sizeV(dataStateIndex)=length(vTarget);
                                    cellArrayThetaTarget{dataStateIndex,1}=thetaTarget;
                                    cellArrayThetaProbe{dataStateIndex,1}=thetaProbe;
                                    cellArrayVtarget{dataStateIndex,1}=vTarget;
                                    cellArrayVprobe{dataStateIndex,1}=vProbe;
                                    
                                elseif systemState==0
                                    
                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATIC DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    
                                    %target
                                    
                                    currDir=[targetRoot{dataStateIndex},filesep,rDtarget{1},filesep,...
                                        dRDirT{dRindexT},filesep,aPtarget{apDirIndxT},filesep,lRtarget{dataStateIndex}];
                                    tmp =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                    
                                    %% bootstrapping
                                    
                                    
                                    paramMatrixTarget{dataStateIndex,1}=tmp.paramMatrix{iboot};
                                    
                                    %probe
                                    
                                    %                             currDir=[sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
                                    %                                dRDir{1}, filesep,aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                    currDir=[sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
                                        dRDir{dRDirIndx},filesep,aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                    tmp =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                    
                                    paramMatrixProbe{dataStateIndex,1}=tmp.paramMatrix;
                                    
                                    % calculate theta and V for the data set
                                    
                                    [thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVaramCovMatrixForMahalonobisDistanceLastFrame(paramMatrixTarget{dataStateIndex,1},paramMatrixProbe{dataStateIndex,1},stateSystemVector(dataStateIndex));
                                    
                                    sizeV(dataStateIndex)=length(vTarget);
                                    cellArrayThetaTarget{dataStateIndex,1}=thetaTarget;
                                    cellArrayThetaProbe{dataStateIndex,1}=thetaProbe;
                                    cellArrayVtarget{dataStateIndex,1}=vTarget;
                                    cellArrayVprobe{dataStateIndex,1}=vProbe;
                                    
                                    
                                end
                                
                                
                                % calculate the combined theta and V for all the entries in stateSystemVector.
                                [ thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal ] = combineStaticDynamicForMahalonobisDistance(cellArrayThetaTarget,cellArrayThetaProbe,...
                                    cellArrayVtarget,cellArrayVprobe,sizeV);
                                
                                % calculate S and p-value
                                
                                [sMatrix(dRDirIndx,aPDirIndx,rDDirIndx,lRIndex),pMatrix(dRDirIndx,aPDirIndx,rDDirIndx,lRIndex)]=calculationMahalonobisDistanceandPvalueDynamicStatic(thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal);
                                %            [sMatrix(rDDirIndx,aPDirIndx,lRIndex)]=calculationMahalonobisDistanceandPvalueDynamicStatic(thetaTargetFinal,thetaProbeFinal,vTargetFinal,vProbeFinal);
                                
                                % save the results
                                
                            end
                        end
                    end
                end
            end
        end
        sMatrixBoot{iboot}=sMatrix;
        pMatrixBoot{iboot}=pMatrix;
        
        resultsDir=[saveDirectory,filesep,rDtarget{1},dRDirT{dRindexT},aPtarget{apDirIndxT},lRtarget{1}];
        mkdir(resultsDir)
        
        save([resultsDir,filesep,'sMatrix_',lRtarget{1,:}],'sMatrix','-v7.3');
        save([resultsDir,filesep,'pMatrix_',lRtarget{1,:}],'pMatrix','-v7.3');
        
        % save a file containing which lRProbe are used for each
        % target and also the stateSystemVector to identify the kind of
        % data.
        save([resultsDir,filesep,'lRandStatesProbe_',lRtarget{1,:}],'lRProbe','stateSystemVector','-v7.3');
        
        
    end
end