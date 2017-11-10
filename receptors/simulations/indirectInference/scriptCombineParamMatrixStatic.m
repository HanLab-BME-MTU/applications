%Script to combine the paramMatrixLastFrame and have the format to add
%information in the intermidiate statistics for the calculation of
%mahalonobis distance
%
%Luciana de Oliveira, February 2017.


sourceRootTargetStatic='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170206';

% name of the target
rDtarget = {'rD100'};%,'rD60','rD80','rD120','rD140','rD160'};
aPtarget = {'aP0p5'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
lRtargetStatic = {'lR0p1';'lR0p08';'lR0p12'};%,'lR0p3','lR0p4','lR0p5'};
numSim=30; %number of simulations

for rDTarIndx = 1 : length(rDtarget)
    
    for aPTarIndx = 1 : length(aPtarget)
    
%target paramMatrix static (remember here we will vary only the labeled fraction)

      %reserve memory for the new paramMatrix for the static data
       paramMatrixStatic=cell(length(lRtargetStatic),1);
      
%       paramMatrixProbeStatic=zeros(1000,numSim);
            
        for lRTarIndxStatic = 1 : length(lRtargetStatic)
% load ratesAndDensity that has the paramMatrix 
currDir=[sourceRootTargetStatic,filesep,rDtarget{rDTarIndx},filesep,...
                aPtarget{aPTarIndx},filesep,lRtargetStatic{lRTarIndxStatic}];
tmp = load(fullfile([currDir,filesep,'paramMatrixLastFrame_',lRtargetStatic{lRTarIndxStatic},'.mat']));

% call paramMatrix and determine target cluster sizes
 
  paramMatrixStatic{lRTarIndxStatic} = tmp.matrixDensity;
         end

% % Allocate space for the matrix with the densities for all 
%        
%        matrixDensity=zeros(1000,length(paramMatrixTarget));
       
       % replace values for each labeled ratio
       % first label ratio
       paramMatrixTargetStatic(1:size(paramMatrixStatic{1},1),:) = paramMatrixStatic{1};
       
        for iLr = 2 : length(lRtargetStatic)     
         
  paramMatrixTargetStatic(size(paramMatrixTargetStatic,1)+1:size(paramMatrixTargetStatic,1)+size(paramMatrixStatic{iLr},1),:) = paramMatrixStatic{iLr};
        end
clear tmp paramMatrixStatic

%% probes load the last frame info
sourceRootProbe ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170206/probe';
%Define strings for directory hierarchy as needed
rDDir = {'rD100'};%,'rD8','rD10','rD12','rD14','rD16'};
aPDir = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
lRprobeStatic = {'lR0p1';'lR0p08';'lR0p12'};%,'lR0p3','lR0p4','lR0p5'};
numSim=30; %number of simulations

fprintf('\n===============================================================');

        
%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        
            %reserve memory for the new paramMatrix for the static data
       paramMatrixStatic=cell(length(lRprobeStatic),1);
      
%       paramMatrixTargetStatic=zeros(1000,numSim);
            
        for lRprobeIndxStatic = 1 : length(lRprobeStatic)
            
% load ratesAndDensity that has the paramMatrix 

currDir=[sourceRootProbe,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRprobeStatic{lRprobeIndxStatic}];
tmp = load(fullfile([currDir,filesep,'paramMatrixLastFrame_',lRprobeStatic{lRprobeIndxStatic},'.mat']));

% call paramMatrix and determine target cluster sizes
 
  paramMatrixStatic{lRprobeIndxStatic} = tmp.matrixDensity;
         end

% % Allocate space for the matrix with the densities for all 
%        
%        matrixDensity=zeros(1000,length(paramMatrixTarget));
       
       % replace values for each labeled ratio
       % first label ratio
       paramMatrixProbeStatic(1:size(paramMatrixStatic{1},1),:) = paramMatrixStatic{1};
       
        for iLr = 2 : length(lRprobeStatic)     
         
  paramMatrixProbeStatic(size(paramMatrixProbeStatic,1)+1:size(paramMatrixProbeStatic,1)+size(paramMatrixStatic{iLr},1),:) = paramMatrixStatic{iLr};
        end
clear tmp paramMatrixStatic    
    end
end
    end 
end
