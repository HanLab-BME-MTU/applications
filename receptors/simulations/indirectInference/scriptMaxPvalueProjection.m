%script to calculate the projection of the maximum value of p-value
%considering the different repetitions of the bootstrapping.

%Output
    % matrixMaxPvalueProjec: maximum p-value for each coordinate [rD,aP,lR]
    
    
currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/bootstrapping/results';
% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
rDtarget = {'rD14'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p7'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p4'};%,'lR0p3','lR0p4','lR0p5'};

   % values of rD, aP and lR of probe
    rDvals = [4;6;8;10;12;14;16];%20;40;60;80;100;140;160];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
   lRStr = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'};%'lR0p01lR0p02';'lR0p02lR0p04';'lR0p03lR0p06
    
%% maximum p-value  projection  
    
% allocate the matrix for max p-value projection
matrixMaxPvalueProj=zeros(length(rDvals),length(aPvals),length(lRStr));

% call p-matrix
temp= load([currDir,filesep,rDtarget{1},aPtarget{1},lRtarget{1},filesep,'pMatrix.mat']);
pMatrix=temp.pMatrix; 
      
 %iterate through the different receptor densities
     for rDTIndx = 1 : length(rDvals)
    
 %Iterate through association probability values per density
     for aPTIndx = 1 : length(aPvals)    
         
%iterate through the different labeling ratios         
     for lRindex = 1 : length(lRStr)
      
   % calculates the maximum p-value projection considering all
   % bootstrapping repetitions
   
   maxPValueproj= max(max(max(pMatrix(rDTIndx,aPTIndx,lRindex,:))));

   %update value in the matrixMaxPvalueProj
   
   if maxPValue >=0.05 
   matrixMaxPvalueProj(rDTIndx,aPTIndx,lRindex)=maxPValueproj;
   end
%% save the information 
   resultsDirBoot=[currDir,filesep,rDtarget{1},aPtarget{1},lRtarget{1}];
   save([resultsDirBoot,filesep,'matrixMaxPvalueProj'],'matrixMaxPvalueProj','-v7.3');  
   
     end 
     end
     end
    
   
   