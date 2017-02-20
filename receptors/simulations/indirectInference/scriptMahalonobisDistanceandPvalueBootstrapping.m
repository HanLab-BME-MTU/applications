%Script to compare target and probe intermediate statistics data througth
%the calculation of Mahalanobis distance and pvalue. s and pValue matrix 
% will be saved in same directory of the data in a file caled Results
%
% The outputs of this script are the sMatrix and pMatrix. They are 3D
% matrices where the rows are receptor densities, the colums are association probability
% and the third dimension is the label ratio.
%
%Luciana de Oliveira, July 2016.

resultsDirectory='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/bootstrapping/results';

%target 
sourceRootTarget ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/bootstrapping/analysis';

% name of the target
rDtarget = {'rD4'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p2'};%,'lR0p3','lR0p4','lR0p5'};
BootNumber=1:100;


%probes
sourceRootProbe ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161201/analysis/probeIS_sT25_dT0p1';
%Define strings for directory hierarchy as needed
%  fprintf('\nProcessing bootstrapping repetition %s' ,bootIndx);
rDDir = {'rD4','rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD60','rD80','rD120','rD140','rD160'}; 
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}; %'dR0p5','dR2p0','dR5p0','aP0p2','aP0p5','aP0p8','dC0p05','dC0p2'};
lRDir = {'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5','lR0p6'};


% call the original pMatrix and sMatrix to have all the information
% together as an unique matrix. The original data will be the first 5
% matrices in the final matrix.

sourceRootOriginalData ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161201/results/target_sT25_dT0p1';
currDir=[sourceRootOriginalData,filesep,rDtarget{1},...
                aPtarget{1},lRtarget{1}];
tmp = load(fullfile(currDir,filesep,'pMatrix.mat'));
pMatrixOrig=tmp.pMatrix;
clear tmp
tmp = load(fullfile(currDir,filesep,'sMatrix.mat'));
sMatrixOrig=tmp.sMatrix;
clear tmp

%define sMatrix and pMatrix
sMatrix=zeros(length(rDDir),length(aPDir),length(lRDir),length(BootNumber)+1);
pMatrix=zeros(length(rDDir),length(aPDir),length(lRDir),length(BootNumber)+1);

% replace the first element of the 4th dimention with the original data
sMatrix(:,:,:,1)=sMatrixOrig;
pMatrix(:,:,:,1)=pMatrixOrig;


for lRDirIndxTarget = 1 : length(lRtarget)
    
    
    %Iterate through association probability values per density
    for aPDirIndxTarget = 1 : length(aPtarget)
        
               
        %iterate through the different labeling ratios
        for rDDirIndxTarget = 1 : length(rDtarget)
            
 fprintf('\nProcessing rD=%s, lR = %s, aP = %s ',rDtarget{rDDirIndxTarget},lRtarget{lRDirIndxTarget},aPtarget{aPDirIndxTarget});

for  bootIndx=1:length(BootNumber)       
% load ratesAndDensity that has the paramMatrix 
currDir=[sourceRootTarget,filesep,rDtarget{1},filesep,...
                aPtarget{1},filesep,lRtarget{1},filesep,'bootstrapping',int2str(BootNumber(bootIndx))];
tmp = load(fullfile(currDir,filesep,'ratesAndDensityComb_dt0p1_T10.mat'));

% call paramMatrix and determine target cluster sizes

paramMatrixTarget=tmp.paramMatrix;
clusterSizeTargetOnRate=size(tmp.rateOnPerClust,1);
clusterSizeTargetOffRate=size(tmp.rateOffPerClust,1);
clusterSizeTargetDensity=size(tmp.densityPerClust,1);


%,'rD6','rD8','rD10','rD12','rD14','rD16'};%,'rD4','rD6','rD8','rD10','rD12','rD14','rD16'};

fprintf('\n===============================================================');


%The top level directory is that of receptor density
for lRDirIndx = 1 : length(lRDir)
    
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing lR = %s, aP = %s ',lRDir{lRDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different labeling ratios
        for rDDirIndx = 1 : length(rDDir)
            
            %name of current directory
            currDir = [sourceRootProbe,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
            
            %load ratesAndDensity that has the paramMatrix 
            tmp = load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
           
%           call paramMatrix and determine probe cluster sizes
            paramMatrixProbe=tmp.paramMatrix;
            clusterSizeProbeOnRate=size(tmp.rateOnPerClust,1);
            clusterSizeProbeOffRate=size(tmp.rateOffPerClust,1);
            clusterSizeProbeDensity=size(tmp.densityPerClust,1);    

            
%call function to compare target and probe data
[sMatrix(rDDirIndx,aPDirIndx,lRDirIndx,bootIndx+1),pMatrix(rDDirIndx,aPDirIndx,lRDirIndx,bootIndx+1)]=calculationMahalonobisDistanceandPvalue(paramMatrixTarget,paramMatrixProbe);      
resultsDirBoot=[resultsDirectory,filesep,rDtarget{rDDirIndxTarget},aPtarget{aPDirIndxTarget},lRtarget{lRDirIndxTarget}];
 
       end
    end
end



%save s and p matrix
%make a directory
save([resultsDirBoot,filesep,'sMatrix'],'sMatrix','-v7.3');
 save([resultsDirBoot,filesep,'pMatrix'],'pMatrix','-v7.3');
end
        end
    end
end