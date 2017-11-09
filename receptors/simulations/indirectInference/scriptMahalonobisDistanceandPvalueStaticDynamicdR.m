% This script is for the calculation of Mahalonobis distance and p-value
% considering multiple entries of dynamic and static data.


%% Input
% TargetRoot for the multiple entries, it is a general script and it can be
% both dynamic and static data. Important note that the targetRoot now is
% write as a cell array, it allow the script to decide how to do the
% calculations for the different system states.

% target source root
targetRoot{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170327/varyDissRate/analysis/targetIS_sT25_dT0p1';
targetRoot{2}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/varyDissRate/target/analysis/targetIS_sT25_dT0p1';
% targetRoot{3}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170313/staticData/target/analysis';
% targetRoot{4}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/target/analysis';
% targetRoot{5}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/target/analysis';
% 
% 
% % probe source root
sourceRootProbe=cell(2,1);
sourceRootProbe{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170501/analysis/probe';

sourceRootProbe{2}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/varyDissRate/analysis/probeIS_sT25_dT0p1';

%  sourceRootProbe{1}='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170313/staticData/probe/analysis';
% sourceRootProbe{4} ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/probe/analysis';
% sourceRootProbe{5} ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170220/probe/analysis';



saveDirectory='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170508/results/combinationDynamicStaticData';
%
%
% % name of the target
rDtarget ={'rD4'};%,'rD20','rD40','rD60','rD80','rD100','rD120','rD140''rD8',
aPtarget = {'aP0p5'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 dRDirT={'dR0p5'};
dRDir={'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};%,'dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};%{'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};
% the labeled fraction is write as a cell array, with each line is
% equivalent for the same set of dynamic and static data and each column is
% for different sets.

lRtarget = {'lR0p2','lR1'};
%     'lR0p4','lR1'};%,'lR0p09','lR0p12','lR0p18'

% this vector identifies the system state and it is equivalent to each
% different entry in the targetRoot. It also have the length of label
% ratio.

stateSystemVector=[1 0];%[1 0];%;


%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBE INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%probes

%Define strings for directo0ry hierarchy as needed
rDDir = {'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%{'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'},{'rD20','rD40','rD60','rD80','rD100','rD120','rD140'} 'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'


% labeled ratio for dynamic and static data
 lRProbe = {%'lR1'
         'lR0p1','lR1'
        'lR0p2','lR1'
          'lR0p3','lR1'
        'lR0p4','lR1'
          'lR0p5','lR1'
          'lR0p6','lR1'
    %            'lR0p04','lR0p08','lR0p09','lR0p16';
    };
%

%% CALCULATIONS


%%%%%%%%%%%%%%%%%%%%%%%%%%% target hierarchy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rDDirIndxT=1:length(rDtarget)
    
    for dRDirIndxT = 1 : length(dRDirT)
        
        for aPDirIndxT = 1 : length(aPtarget)
            
            for  lRIndexT=1: length(lRtarget)
              
               
                
                %%%%%%%%%%%%%%%%%%%% probe hierarchy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for dRDirIndx = 1 : length(dRDir)
                     fprintf('\nProcessing target rD = %s, dRT=%s,  aP = %s, dRP=%s ',rDtarget{rDDirIndxT},dRDirT{dRDirIndxT},aPtarget{aPDirIndxT },dRDir{dRDirIndx}); 
                    
                    for rDDirIndx = 1 : length(rDDir)
                        
                        %Iterate through association probability values per density
                        for aPDirIndx = 1 : length(aPDir)
                            
                            for lRIndex=1:size(lRProbe,1)
                                
                                for dataStateIndex=1: length(stateSystemVector)
                                    
                                    %define if it is dynamic or static data
                                    systemState=stateSystemVector(dataStateIndex);
                                    
                                     paramMatrixTarget=cell(length(stateSystemVector),1);
                                     paramMatrixProbe=cell(length(stateSystemVector),1);
                                    
                                    if systemState==1
                                        
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% DYNAMIC DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
                                        % target
                                        % load target paramMatrix
                                        
                                        % as each entry of the stateSystemVector coresponds to
                                        % one different value of label ratio, here the lRtarget
                                        % index is given by the dataStateIndex.
                                        
                                         currDir=[targetRoot{dataStateIndex},filesep,rDtarget{rDDirIndxT},filesep,...
                                             dRDirT{dRDirIndxT},filesep,aPtarget{aPDirIndxT },filesep,lRtarget{dataStateIndex}];
%                                                                  currDir=[targetRoot{1},filesep,rDtarget{1},filesep,...
%                                                                    aPtarget{1},filesep,lRtarget{dataStateIndex}];
                                        tmp =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                        paramMatrixTarget{dataStateIndex,1}=tmp.paramMatrix;
                                        
                                        %probe
                                        
                                        % load probe paramMatrix
                                        
                                        % Now the lRprobe is an array, and each row corresponds
                                        % to a different entry stateSystemVector, so to have the lRprobe values
                                        % corresponding to each probe we need to index for the lR.
                                        
%                                            currDir = [sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
%                                              aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                        
                                          currDir = [sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
                                           dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
%                                         
                                        %load ratesAndDensity that has the paramMatrix
                                        tmp = load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                        paramMatrixProbe{dataStateIndex,1}=tmp.paramMatrix;
                                        
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
                                        
                                      currDir=[targetRoot{dataStateIndex},filesep,rDtarget{rDDirIndxT},filesep,...
                                           dRDirT{dRDirIndxT},filesep,aPtarget{aPDirIndxT },filesep,lRtarget{dataStateIndex}];
%                                         currDir=[targetRoot{dataStateIndex},filesep,rDtarget{1},filesep,...
%                                                                    aPtarget{1},filesep,lRtarget{dataStateIndex}];
%                                         
                                        tmp =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                         paramMatrixTarget{dataStateIndex,1}=tmp.paramMatrix;
% paramMatrixTarget=tmp.paramMatrix;
                                        
                                        %probe
                                        
%                                          currDir=[sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
%                                            aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                        
                                        
                                      currDir=[sourceRootProbe{dataStateIndex},filesep,rDDir{rDDirIndx},filesep,...
                                           dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,lRProbe{lRIndex,dataStateIndex}];
                                        tmp =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
                                        
                                       paramMatrixProbe{dataStateIndex,1}=tmp.paramMatrix;
%                                          paramMatrixProbe=tmp.paramMatrix;
                                        % calculate theta and V for the data set
                                        
                                         [thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVaramCovMatrixForMahalonobisDistanceLastFrame(paramMatrixTarget{dataStateIndex,1},paramMatrixProbe{dataStateIndex,1},stateSystemVector(dataStateIndex));
%                                          [thetaTarget,thetaProbe,vTarget,vProbe]=calcThetaVaramCovMatrixForMahalonobisDistanceLastFrame(paramMatrixTarget,paramMatrixProbe,stateSystemVector(dataStateIndex));
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
                resultsDir=[saveDirectory,filesep,rDtarget{1},dRDirT{dRDirIndxT},aPtarget{aPDirIndxT},lRtarget{1},lRtarget{2}];
                      mkdir(resultsDir)
                    
                    save([resultsDir,filesep,'sMatrix_',lRtarget{1,:}],'sMatrix','-v7.3');
                    save([resultsDir,filesep,'pMatrix_',lRtarget{1,:}],'pMatrix','-v7.3');
                    
                    % save a file containing which lRProbe are used for each
                    % target and also the stateSystemVector to identify the kind of
                    % data.
                    
                    save([resultsDir,filesep,'lRandStatesProbe_',lRtarget{1,:}],'lRProbe','stateSystemVector','-v7.3');  
                
            end
        end
    end
end

