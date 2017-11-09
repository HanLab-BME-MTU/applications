function [maxPmatrix,minSmatrix,maxPvalueOriginal,minSvalueOriginal,pValueMean,sValueMean,pMean]=bootstrapAnalysis
% BOOTSTRAPANALYSIS is a function for the analysis of the bootstraping . It
% will take the resulting p-value and Mahalanobis distance and do the
% analysis of eficience of the model
%
% Luciana de Oliveira, december 2016
%% Paths information

%directory for the original target
targetDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160817/results/S_PMatrix_sT25_dT0p1';

% directory with pValue and Svalues matrices for the bootstrapping
currDir = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161117/results';
rDtarget = {'rD14'};%,'rD8','rD10','rD12','rD14','rD16'};
aPtarget = {'aP0p3'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
lRtarget ={'lR0p2'};%,'lR0p3','lR0p4','lR0p5'};
lRStr = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5'};
BootNumber=1:100;


% Load values and build matrices with the values for all boortstrap repetition

%% calculations for the original target

temp= load([targetDir,filesep,'target',rDtarget{1},aPtarget{1},lRtarget{1},filesep,'pMatrix.mat']);
pMatrix=temp.pMatrix;
temp= load([targetDir,filesep,'target',rDtarget{1},aPtarget{1},lRtarget{1},filesep,'sMatrix.mat']);
sMatrix=temp.sMatrix;

%first calculate the maximum value for each receptor label

maxPvalue=max(max(pMatrix));
minSvalue=max(min(sMatrix));

%to have the global maximum value we need to %put this values together as 
%a line vector and then calculates the maximum limite of p or s.
maxPvalueOriginal=reshape(maxPvalue,1,length(maxPvalue));
minSvalueOriginal=reshape(minSvalue,1,length(minSvalue));

clear temp maxPvalue minSvalue pMatrix sMatrix

%% calculations for the bootstrapping

% pre allocating variables
maxPmatrix=zeros(length(BootNumber),length(lRStr));
minSmatrix=zeros(length(BootNumber),length(lRStr));



for  bootIndx=1:length(BootNumber)
%load p matrix    
temp= load([currDir,filesep,'bootstrapping',int2str(BootNumber(bootIndx)),filesep,rDtarget{1},aPtarget{1},lRtarget{1},filesep,'pMatrix.mat']);
pMatrix=temp.pMatrix;
    
%load s matrix
temp= load([currDir,filesep,'bootstrapping',int2str(BootNumber(bootIndx)),filesep,rDtarget{1},aPtarget{1},lRtarget{1},filesep,'sMatrix.mat']);
sMatrix=temp.sMatrix;


%load dof
 temp= load([currDir,filesep,'bootstrapping',int2str(BootNumber(bootIndx)),filesep,rDtarget{1},aPtarget{1},lRtarget{1},filesep,'dof.mat']);
 dof=temp.dof;

%first calculate the maximum value for each receptor label

maxPvalue=max(max(pMatrix));
minSvalue=max(min(sMatrix));

%to have the global maximum value we need to %put this values together as 
%a line vector and then calculates the maximum limite of p or s.
maxPvector=reshape(maxPvalue,1,length(maxPvalue));
minSvector=reshape(minSvalue,1,length(minSvalue));

% build the matrix with rows corresponding to each
maxPmatrix(bootIndx,:)=maxPvector;

minSmatrix(bootIndx,:)=minSvector;

dofVector(bootIndx)=dof;
end % for each bootstrap number

%% mean values

pValueMean=mean(maxPmatrix,1);
sValueMean=mean(minSmatrix,1);
dofMean=round(mean(dofVector));


%% New pvalue calculated from mean S distance
for   i=1:length(lRStr)
pMean(i) = 1 - chi2cdf(sValueMean(i),dofMean);
end


%% Percentage of p-value greater than 0.95

percentPvalue=sum(maxPmatrix>=0.95);

%% Save results
resultsDirBoot=[currDir,filesep,'finalResults',filesep,rDtarget{1},aPtarget{1},lRtarget{1}];
mkdir(resultsDirBoot)
save([resultsDirBoot,filesep,'maxPvalueOriginal'],'maxPvalueOriginal','-v7.3');
save([resultsDirBoot,filesep,'minSvalueOriginal'],'minSvalueOriginal','-v7.3');
save([resultsDirBoot,filesep,'pValueMean'],'pValueMean','-v7.3');
save([resultsDirBoot,filesep,'sValueMean'],'sValueMean','-v7.3');
save([resultsDirBoot,filesep,'pMean'],'pMean','-v7.3');
save([resultsDirBoot,filesep,'percentPvalue'],'percentPvalue','-v7.3');
end

