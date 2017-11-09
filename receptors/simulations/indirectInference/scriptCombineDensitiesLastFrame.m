%Script to combine density results from different equivalent
%"movies". will be saved in same directory as cell
%array of individual movie results
%
%Luciana de Oliveira, February 2017.

sourceRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170206';

%Define strings for directory hierarchy as needed
 rDDir = {'rD100'};%,'rD60','rD80','rD120','rD140','rD160'};
 aPDir = {'aP0p5'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
 outDirNum =1:30;
 lRDir = {'lR0p1'};%;'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'
fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
           %iterate through the different runs
    
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx},filesep,'ind'];
          
            %allocate cell array 
            
          densityPerMovie=cell(length(outDirNum),1);
            
            for outDirIndx = 1 : length(outDirNum)
                
                %read results
                tmp = load(fullfile(currDir,['clusterStatsLastFrame_',int2str(outDirNum(outDirIndx)),'.mat']));
                
                %store in cell array
             densityPerMovie{outDirIndx}= tmp.clusterStatsLastFrame.clusterDensity;
                
            end %for each outDir
            
     %calculate the maximum cluster size
     maxSizeDensity=zeros(1,length(outDirNum));
         for iMovie = 1 : length(outDirNum)
     maxSizeDensity(iMovie) = length(densityPerMovie{iMovie});
   
         end
       %calculate the max value for the cluster size of all movies
       maxSizeDensityAll = max(maxSizeDensity);

       
       % Allocate space for the matrix with the densities for all 
       
       matrixDensity=zeros(maxSizeDensityAll,length(outDirNum));
       
       % replace values for each movie
       
        for iMovie = 1 : length(outDirNum)     
  matrixDensity(1:length(densityPerMovie{iMovie}),iMovie) = densityPerMovie{iMovie};
        end
        
             %save directory
            
           saveDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
            
            %save combined results
            % save in the same format as function combineClusterRatesAndDensity
            
            save([saveDir,filesep,'paramMatrixLastFrame_',lRDir{lRDirIndx}],'matrixDensity','-v7.3');
            

        end %for each labelRatio
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear

