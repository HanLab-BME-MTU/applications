%This function calculates the intermediate statistics and the final output
%is in the format to calculate the Mahalanobis distance and 

%% INPUT

sourceRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2017/01/20170120/varyDissRate/probe';
% saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170510';
destinationRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170501';

%Define strings for directory hierarchy as needed
  rDDir = {'rD8'}; %{'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
  aPDir = {'aP0p7','aP0p8'}; %{'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 dRDir={'dR0p5','dR2'};%,'dR1p5','dR1p75'
outDirNum =1:10;
lRDir = {'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5','lR0p6'};%'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5','lR0p6'
intensityQuantum=[1 0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: pay attention for the 'areaSideLen','timeStep',and 'sampleStep'
%values. For the low density simulation 'areaSideLen=25' and for high densities
%it is 12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoSpaceTime = struct('probDim',2,'areaSideLen',25,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10]);


%% Calculations

%step 1: calculate rates and density

%   functionExtractMultipleCompTracks(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
% 
%   functionCalcMultipleAggregState(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,intensityQuantum)
% 
%           functionCalcMultipleRatesDensity(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)
%step 4: move rates and density
% % 
    functionMoveMultipleRatesDensity(sourceRoot,destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
% % % % % % % 
% % % % % % % %step 5: collect rates and density
           functionCollectMultipleRatesDensity(destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
%step 6: combine rates and density
% % % % % % % % 
    functionCombineMultipleRatesDensity(destinationRoot,rDDir,aPDir,dRDir,lRDir)
% % % % 
