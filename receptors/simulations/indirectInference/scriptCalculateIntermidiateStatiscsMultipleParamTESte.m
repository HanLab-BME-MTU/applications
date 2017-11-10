%Script to calculate the intermediate statistics of multiple data

%% INPUT

currFile ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170501';

%load the cell array with all the information
temp=load([currFile,filesep,'cellInfoAllSimProbes.mat']);
cellInfoAllSimProbe=temp.cellInfoAllProbes;

for cellIndex=1: size(cellInfoAllSimProbe,1)

sourceRoot =cellInfoAllSimProbe{cellIndex,1};

%Define strings for directory hierarchy as needed
 rDDir =cellInfoAllSimProbe{cellIndex,2};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
 aPDir =  cellInfoAllSimProbe{cellIndex,3};%'aP0p6','aP0p7','aP0p8' ,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 dRDir=cellInfoAllSimProbe{cellIndex,4};%,'dR1p5','dR1p75'
outDirNum =1:10;
lRDir = cellInfoAllSimProbe{cellIndex,5};

intensityQuantum=[1 0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: pay attention for the 'areaSideLen','timeStep',and 'sampleStep'
%values. For the low density simulation 'areaSideLen=25' and for high densities
%it is 12.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infoSpaceTime = struct('probDim',2,'areaSideLen',25,'timeStep',0.01,'sampleStep',0.1,'firstLastTP',[0 10]);

%% Calculations


%step 1: extract compTracks:

 functionExtractMultipleCompTracks(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)

%step 2: calculate aggregation state 

 functionCalcMultipleAggregState(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,intensityQuantum)

%step 3: calculate rates and density

 functionCalcMultipleRatesDensity(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,infoSpaceTime)
 
end

