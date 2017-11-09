%Script to compare target and probe intermediate statistics data througth
%the calculation of Mahalanobis distance and pvalue. s and pValue matrix 
% will be saved in same directory of the data in a file caled Results
%
% The outputs of this script are the sMatrix and pMatrix. They are 3D
% matrices where the rows are receptor densities, the colums are association probability
% and the third dimension is the label ratio.
%
%Luciana de Oliveira, July 2016.

resultsDirectory='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170207/results';


%target paramMatrix dynamic

sourceRootTarget ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/10/20161004/analysis';
% name of the target
rDtarget = {'rD100'};%,'rD60','rD80','rD120','rD140','rD160'};
aPtarget = {'aP0p5'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
lRtarget = {'lR0p03'};%,'lR0p3','lR0p4','lR0p5'};



for rDTarIndx = 1 : length(rDtarget)
    
    for aPTarIndx = 1 : length(aPtarget)
        
        for lRTarIndx = 1 : length(lRtarget)
% load ratesAndDensity that has the paramMatrix 
currDir=[sourceRootTarget,filesep,rDtarget{rDTarIndx},filesep,...
                aPtarget{aPTarIndx},filesep,lRtarget{lRTarIndx}];
tmp = load(fullfile(currDir,filesep,'ratesAndDensityComb_dt0p1_T10.mat'));

% call paramMatrix and determine target cluster sizes

paramMatrixTarget=tmp.paramMatrix;
clusterSizeTargetOnRate=size(tmp.rateOnPerClust,1);
clusterSizeTargetOffRate=size(tmp.rateOffPerClust,1);
clusterSizeTargetDensity=size(tmp.densityPerClust,1);

%probes
sourceRootProbe ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161201/analysis/probeIS_sT10_dT0p1';
%Define strings for directory hierarchy as needed
rDDir = {'rD20','rD40','rD60','rD80','rD100','rD120','rD140'};%,'rD8','rD10','rD12','rD14','rD16'};
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
lRDir = {'lR0p03'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
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

           %% reform param matrix dynamic
[paramMatrixTargetReform,paramMatrixProbeReform]=reformParamMatrix(paramMatrixTarget,paramMatrixProbe);
          
            %% static 

%target paramMatrix static (remember here we will vary only the labeled fraction)


lRtargetStatic = {'lR0p06'};%,'lR0p3','lR0p4','lR0p5'};
%number of simulations
    
%target paramMatrix static (remember here we will vary only the labeled fraction)

      %reserve memory for the new paramMatrix for the static data
%        paramMatrixStatic=cell(length(lRtargetStatic),1);
      
%       paramMatrixProbeStatic=zeros(1000,numSim);
            
%         for lRTarIndxStatic = 1 : length(lRtargetStatic)
% load ratesAndDensity that has the paramMatrix 
 currDir=[sourceRootTargetStatic,filesep,rDtarget{rDTarIndx},filesep,...
                 aPtarget{aPTarIndx},filesep,lRtargetStatic{lRTarIndxStatic}];
tmp = load(fullfile([currDir,filesep,'paramMatrixLastFrame_',lRtargetStatic{lRTarIndxStatic},'.mat']));


% call paramMatrix and determine target cluster sizes
 
  paramMatrixStaticTarget= tmp.matrixDensity;

%          end

% % Allocate space for the matrix with the densities for all 
%        
      
       % replace values for each labeled ratio
       % first label ratio
%        paramMatrixTargetStatic(1:size(paramMatrixStatic{1},1),:) = paramMatrixStatic{1};
%        
%         for iLr = 2 : length(lRtargetStatic)     
%          
%   paramMatrixTargetStatic(size(paramMatrixTargetStatic,1)+1:size(paramMatrixTargetStatic,1)+size(paramMatrixStatic{iLr},1),:) = paramMatrixStatic{iLr};
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROBE

 lRprobeStatic = {'lR0p03'};%,'lR0p3','lR0p4','lR0p5'};

fprintf('\n===============================================================');
        
            %reserve memory for the new paramMatrix for the static data
%        paramMatrixStatic=cell(length(lRprobeStatic),1);
      
%       paramMatrixTargetStatic=zeros(1000,numSim);
            
%         for lRprobeIndxStatic = 1 : length(lRprobeStatic)
            
% load ratesAndDensity that has the paramMatrix 

 currDir=[sourceRootProbeStatic,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRprobeStatic{lRprobeIndxStatic}];

            tmp = load(fullfile([currDir,filesep,'paramMatrixLastFrame_',lRprobeStatic{lRprobeIndxStatic},'.mat']));

          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% call paramMatrix and determine target cluster sizes
 
   paramMatrixStaticProbe = tmp.matrixDensity;

%          end

% % Allocate space for the matrix with the densities for all 
%        
%        matrixDensity=zeros(1000,length(paramMatrixTarget));
       
       % replace values for each labeled ratio
       % first label ratio
%        paramMatrixProbeStatic(1:size(paramMatrixStatic{1},1),:) = paramMatrixStatic{1};
%        
%         for iLr = 2 : length(lRprobeStatic)     
%          
%   paramMatrixProbeStatic(size(paramMatrixProbeStatic,1)+1:size(paramMatrixProbeStatic,1)+size(paramMatrixStatic{iLr},1),:) = paramMatrixStatic{iLr};
%         end
    

%% reform param matrix static

%update the value of cluster size for the matrix with 
 maxClusterSizeTarget=size(paramMatrixStaticTarget,1);
 maxClusterSizeProbe=size(paramMatrixStaticProbe,1);

 if maxClusterSizeTarget ~= maxClusterSizeProbe
     
      %calculates the difference between the sizes
     
     clusterDiffDensity=abs(maxClusterSizeTarget-maxClusterSizeProbe);
    
    if (maxClusterSizeTarget > maxClusterSizeProbe);
    paramMatrixStaticProbe = [paramMatrixStaticProbe;zeros( clusterDiffDensity,size(paramMatrixStaticProbe,2))];
    elseif (maxClusterSizeTarget < maxClusterSizeProbe)
    paramMatrixStaticTarget=[paramMatrixStaticTarget; zeros(clusterDiffDensity,size(paramMatrixStaticTarget,2))];
    end
 end
   
                     

%% put matrices together
 paramMatrixProbeNew=[paramMatrixProbeReform; paramMatrixStaticProbe];   
 paramMatrixTargetNew=[paramMatrixTargetReform; paramMatrixStaticTarget];   

  %% calculation           
%call function to compare target and probe data
[sMatrix(rDDirIndx,aPDirIndx,lRDirIndx),pMatrix(rDDirIndx,aPDirIndx,lRDirIndx)]=calculationMahalonobisDistanceandPvalueLastFrame(paramMatrixTargetNew,paramMatrixProbeNew);      


       end
    end
end
%save s and p matrix
%make a directory
resultsDirBoot=[resultsDirectory,filesep,rDtarget{rDTarIndx},aPtarget{aPTarIndx},lRtarget{lRTarIndx}];
 mkdir(resultsDirBoot)
 save([resultsDirBoot,filesep,'sMatrix'],'sMatrix','-v7.3');
 save([resultsDirBoot,filesep,'pMatrix'],'pMatrix','-v7.3');
        end 
    end
end
