%Script to calculate the intermediate statistics of superresolution multiple data 

%% INPUT
originRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170306/varyDissRate/probe';

sourceRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170517/superRes/varyLabelFraction/probe';

destinationRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170517/superRes/varyLabelFraction/analysis/probeIS_sT25_dT0p1';

%Define strings for directory hierarchy as needed
  rDDir = {'rD4'}; %{'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
  aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}; %{'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 dRDir={'dR1p25'};%,'dR1p5','dR1p75'
outDirNum =1:10;
lRDir = {'lR0p91';'lR0p93';'lR0p95';'lR0p97';'lR0p99'};
labelRatio=0.91:0.02:0.99;
intensityQuantum=[1 0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: pay attention for the 'areaSideLen','timeStep',and 'sampleStep'
%values. For the low density simulation 'areaSideLen=25' and for high densities
%it is 12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoSpaceTime = struct('probDim',2,'areaSideLen',25,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10],'systemState',0);

%% Calculations


%step 1: calculate detection
  getMultipleAggregStateLabeldReceptors(originRoot,sourceRoot,rDDir,dRDir,aPDir,outDirNum,lRDir,labelRatio);

%step 2: calculate rates and density

  functionCalcMultipleRatesDensityReceptorInfoAllStatic(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)

%step 3: move rates and density

functionMoveMultipleRatesDensity(sourceRoot,destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)

%step 4: collect rates and density
 
functionCollectMultipleRatesDensity(destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)

%step 5: combine rates and density

functionCombineMultopleRatesDensityStaticDynamicdR(destinationRoot,rDDir,aPDir,dRDir,lRDir)
