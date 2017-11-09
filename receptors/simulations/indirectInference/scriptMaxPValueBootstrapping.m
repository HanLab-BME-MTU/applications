%script to determine the maximum value of p-value for each repetition of
%bootstraping 
%
%Output
    % matrixMaxPvalueCount: counts of the max p-value for the bootstrapping repetition
    % maxPvalueCoord: cell array with the first column the max p-values and
    % the second column their coordinates [rD,aP,lR]
    % maxPvalueCoordSig: same as maxPvalueCoord, but considering only
    % p-value >=0.05.
    
    
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
    lRvals = [0.1;0.2;0.3;0.4;0.5;0.6];%'lR0p01lR0p02';'lR0p02lR0p04';'lR0p03lR0p06
    BootNumber=1:101;

%% maximum p-value    
    
% allocate the cell array for the p-value and respective coordinates
maxPvalueCoord=cell(length(BootNumber),2);
maxPvalueCoordSig=cell(length(BootNumber),2);
    
    
    
    for lRindx=1:length(lRvals)
     
     for lRTIndx = 1 : length(lRtarget)

     %Iterate through association probability values per density
     for aPTIndx = 1 : length(aPtarget)
%         
%                
%         %iterate through the different labeling ratios
     for rDTIndx = 1 : length(rDtarget)

    %load pvalue matrix

temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx},filesep,'pMatrix.mat']);
pMatrix=temp.pMatrix;     
 

% Allocate space for the maximum p-value matrix

matrixMaxPValueCount=zeros(size(pMatrix(:,:,:,1)));

 for  bootIndx=1:length(BootNumber)    
        
   % calculates the maximum p-value for each repetition of the bootstrapping.
   % It will give the maximum value of pvalue considering the 3 dimensions:
   % (rD,aP,lR)
   maxPValue = max(max(max(pMatrix(:,:,:,bootIndx))));

    %find the index of the maximum p-value
    [rDindex,aPindex,lRindex]=ind2sub(size(pMatrix), find(pMatrix(:,:,:,bootIndx)==maxPValue));
 
    
% I am considering two formats for the cell with the max p-value and the respective
% coordinate [rDindex,aPindex,lRindex]:
% maxPvalueCoord will have all the p-values and maxPvalueCoordSig will have
% only the significant values.
%PS: the count considers only the significant values. 

   maxPvalueCoord(bootIndx,1)={maxPValue};
   maxPvalueCoord(bootIndx,2)={[rDindex aPindex lRindex]};
   
% for the significant values

  if maxPValue >=0.05 
      
   %count max p-value for in the matrix of max p-value
   matrixMaxPValueCount(rDindex,aPindex,lRindex)=matrixMaxPValueCount(rDindex,aPindex,lRindex)+1; 
   
   %pvalue and coordinates
   maxPvalueCoordSig(bootIndx,1)={maxPValue};
   maxPvalueCoordSig(bootIndx,2)={[rDindex aPindex lRindex]};
   
  end


  
%% save the information 
   resultsDirBoot=[currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx}];
   save([resultsDirBoot,filesep,'matrixMaxPvalueCount'],'matrixMaxPValueCount','-v7.3');  
   save([resultsDirBoot,filesep,'maxPvalueCoord'],'maxPvalueCoord','-v7.3');  
   save([resultsDirBoot,filesep,'maxPvalueCoordSig'],'maxPvalueCoordSig','-v7.3');  
   
 end 
     end
     end
     end
    end
   
   