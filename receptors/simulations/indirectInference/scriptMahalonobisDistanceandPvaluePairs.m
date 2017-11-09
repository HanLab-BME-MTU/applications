%Script to compare target and probe intermediate statistics data througth
%the calculation of Mahalanobis distance and pvalue. s and pValue matrix 
% will be saved in same directory of the data in a file caled Results
%
% The outputs of this script are the sMatrix and pMatrix. They are 3D
% matrices where the rows are receptor densities, the colums are association probability
% and the third dimension is the label ratio.
%
%Luciana de Oliveira, July 2016.

target1Root='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161209/pairComparison/target/analysis';
target2Root='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161208/DiffLabelRatio/analysis/target';
saveDirectory='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/pairComparison/results';


% name of the target
rDtarget = {'rD100'};%,'rD8','rD10','rD12','rD14','rD16'};
aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
lRtarget ={'lR0p06','lR0p12'};%,'lR0p3','lR0p4','lR0p5'};



    
    tic
    %Iterate through association probability values per density
    for aPTIndx = 1 : length(aPtarget)
        
               
        %iterate through the different labeling ratios
        for rDTIndx = 1 : length(rDtarget)

% load ratesAndDensity that has the paramMatrix for the target 1
currDir=[target1Root,filesep,rDtarget{rDTIndx},filesep,...
                aPtarget{aPTIndx},filesep,lRtarget{1}];
tmp1 =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));

% call paramMatrix and determine target cluster sizes

paramMatrixTarget1=tmp1.paramMatrix;

currDir=[target2Root,filesep,rDtarget{rDTIndx},filesep,...
                aPtarget{aPTIndx},filesep,lRtarget{2}];
            
tmp2 =  load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
% call paramMatrix and determine target cluster sizes

paramMatrixTarget2=tmp2.paramMatrix;


% load ratesAndDensity that has the paramMatrix for the probe 1
sourceRootProbe1 ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161201/analysis/probeIS_sT10_dT0p1';
sourceRootProbe2 ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161208/DiffLabelRatio/analysis/probe';
%Define strings for directory hierarchy as needed
 fprintf('\nProcessing target %s',rDtarget{1});
rDDir = {'rD20','rD40','rD60','rD80','rD100','rD120','rD140',};%,{'rD20','rD40','rD60','rD80','rD120','rD140','rD160'}; 
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}; %'dR0p5','dR2p0','dR5p0','aP0p2','aP0p5','aP0p8','dC0p05','dC0p2'};
lRDir = {'lR0p06','lR0p12'};
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
            currDir1 = [sourceRootProbe1,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{1}];
            
            %load ratesAndDensity that has the paramMatrix 
            tmp1 = load(fullfile(currDir1,'ratesAndDensityComb_dt0p1_T10.mat'));
           
%           call paramMatrix and determine probe cluster sizes
            paramMatrixProbe1=tmp1.paramMatrix;
            
            %name of current directory
            currDir2 = [sourceRootProbe2,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{2}];
            
            %load ratesAndDensity that has the paramMatrix 
            tmp2 = load(fullfile(currDir2,'ratesAndDensityComb_dt0p1_T10.mat'));
           
%           call paramMatrix and determine probe cluster sizes
            paramMatrixProbe2=tmp2.paramMatrix;
           
%calculate theta and v for set 1:

[thetaTarget1,thetaProbe1,vTarget1,vProbe1]=calcThetaVaramCovMatrixForMahalonobisDistance(paramMatrixTarget1,paramMatrixProbe1);
  
%calculate theta and v for set 2:

[thetaTarget2,thetaProbe2,vTarget2,vProbe2]=calcThetaVaramCovMatrixForMahalonobisDistance(paramMatrixTarget2,paramMatrixProbe2);


%% new variance covariance matrix

%number of rows is equal number of rows v1 plus number of rows v2

%target
numRowsT=size(vTarget1,1)+size(vTarget2,1);
vTarget12=zeros(numRowsT);

%fill v1 and v2 into v12

vTarget12(1:size(vTarget1,1),1:size(vTarget1,1))=vTarget1;
vTarget12(size(vTarget1,1)+1:end,size(vTarget1,1)+1:end)=vTarget2;

% thetas
thetaTarget12=[thetaTarget1;thetaTarget2];

%probe
numRowsP=size(vProbe1,1)+size(vProbe2,1);
vProbe12=zeros(numRowsP);

%fill v1 and v2 into v12

vProbe12(1:size(vProbe1,1),1:size(vProbe1,1))=vProbe1;
vProbe12(size(vProbe1,1)+1:end,size(vProbe1,1)+1:end)=vProbe2;

%thetas

thetaProbe12=[thetaProbe1;thetaProbe2];

%calculate mahalonobis distance

[sMatrix(rDDirIndx,aPDirIndx),pMatrix(rDDirIndx,aPDirIndx)]=calculationMahalonobisDistanceandPvaluePairs(thetaTarget12,thetaProbe12,vTarget12,vProbe12);

savDir=[saveDirectory,filesep,rDtarget{1},aPtarget{1},lRtarget{1},lRtarget{2}];
mkdir(savDir)
   save([savDir,filesep,'sMatrix',lRDir{1},lRDir{2}],'sMatrix','-v7.3');
  save([savDir,filesep,'pMatrix',lRDir{1},lRDir{2}],'pMatrix','-v7.3');


        end
    end
end
        end
    end

%save s and p matrix
%make a directory
